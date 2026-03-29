"""
NeuralSDF training script.

Usage examples:

  # Train on SCONE binary output:
  python train.py --input sphere.bin --output sphere_weights.bin

  # Train on Python-native sphere SDF (no SCONE required):
  python train.py --sphere --radius 5.0 --output sphere_weights.bin

  # Common overrides:
  python train.py --input shape.bin --output shape.bin \\
      --hidden-dim 64 --num-layers 3 --epochs 200 --lr 1e-3

Outputs:
  <output>          Binary weight file for SCONE (mlpWeightIO_mod.f90)
  <output>.txt      Text weight file (for debugging)
  <output>.pt       PyTorch checkpoint (for resuming / inspection)
"""

import argparse
import sys
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from model   import NeuralSDF
from sampler import (load_scone_binary, generate_sphere_sdf,
                     SdfDataset, make_train_val_split)
from export  import export_weights, export_text


# ---------------------------------------------------------------------------
# Training loop
# ---------------------------------------------------------------------------

def train(model, train_loader, val_loader, epochs, lr, weight_decay, device, verbose):
    """
    Train model with Adam + L1 loss (DeepLS Sec. 4.2).

    Returns:
        history : list of dicts {'epoch', 'train_loss', 'val_loss'}
    """
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    criterion = nn.L1Loss()

    history = []
    for epoch in range(1, epochs + 1):
        # --- Training ---
        model.train()
        train_loss = 0.0
        for coords, sdfs in train_loader:
            coords = coords.to(device)
            sdfs   = sdfs.to(device)
            optimizer.zero_grad()
            pred = model(coords)
            loss = criterion(pred, sdfs)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * len(coords)
        train_loss /= len(train_loader.dataset)

        # --- Validation ---
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for coords, sdfs in val_loader:
                coords = coords.to(device)
                sdfs   = sdfs.to(device)
                pred = model(coords)
                val_loss += criterion(pred, sdfs).item() * len(coords)
        val_loss /= len(val_loader.dataset)

        history.append({'epoch': epoch, 'train_loss': train_loss, 'val_loss': val_loss})

        if verbose and (epoch % max(1, epochs // 20) == 0 or epoch == 1):
            print(f"  Epoch {epoch:4d}/{epochs}  "
                  f"train={train_loss:.6f}  val={val_loss:.6f}")

    return history


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Train a NeuralSDF MLP and export weights for SCONE.')

    # Data source (mutually exclusive)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument('--input', metavar='FILE',
                     help='SCONE sdfSampler binary output file')
    src.add_argument('--sphere', action='store_true',
                     help='Generate training data from analytic sphere SDF')

    # Sphere parameters (only used with --sphere)
    p.add_argument('--radius', type=float, default=1.0,
                   help='Sphere radius (default 1.0, only with --sphere)')
    p.add_argument('--center', type=float, nargs=3, default=[0.0, 0.0, 0.0],
                   metavar=('CX', 'CY', 'CZ'),
                   help='Sphere center (default 0 0 0, only with --sphere)')
    p.add_argument('--bbox-min', type=float, nargs=3,
                   metavar=('XMIN', 'YMIN', 'ZMIN'),
                   help='Bounding box minimum (default: center - 1.5*radius)')
    p.add_argument('--bbox-max', type=float, nargs=3,
                   metavar=('XMAX', 'YMAX', 'ZMAX'),
                   help='Bounding box maximum (default: center + 1.5*radius)')
    p.add_argument('--n-samples', type=int, default=100_000,
                   help='Number of training samples for --sphere (default 100000)')

    # Output
    p.add_argument('--output', required=True, metavar='FILE',
                   help='Output binary weight file path (e.g. sphere_weights.bin)')

    # Architecture
    p.add_argument('--hidden-dim', type=int, default=128,
                   help='Hidden layer width (default 128)')
    p.add_argument('--num-layers', type=int, default=4,
                   help='Total number of weight matrices including output (default 4)')
    p.add_argument('--activation', choices=['leakyrelu', 'relu', 'tanh'],
                   default='leakyrelu',
                   help='Hidden layer activation (default leakyrelu)')
    p.add_argument('--leaky-alpha', type=float, default=0.01,
                   help='Negative slope for LeakyReLU (default 0.01)')

    # Training
    p.add_argument('--epochs', type=int, default=100,
                   help='Number of training epochs (default 100)')
    p.add_argument('--lr', type=float, default=1e-3,
                   help='Adam learning rate (default 1e-3)')
    p.add_argument('--weight-decay', type=float, default=1e-6,
                   help='Adam weight decay (default 1e-6)')
    p.add_argument('--batch-size', type=int, default=4096,
                   help='Mini-batch size (default 4096)')
    p.add_argument('--val-fraction', type=float, default=0.1,
                   help='Validation set fraction (default 0.1)')
    p.add_argument('--seed', type=int, default=42,
                   help='Random seed (default 42)')
    p.add_argument('--float32', action='store_true',
                   help='Train in float32 (default: float64)')

    # Misc
    p.add_argument('--quiet', action='store_true',
                   help='Suppress per-epoch progress output')

    return p.parse_args()


def main():
    args = parse_args()
    verbose = not args.quiet
    dtype = torch.float32 if args.float32 else torch.float64

    torch.manual_seed(args.seed)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    if verbose:
        print(f"Device: {device}")

    # ------------------------------------------------------------------
    # Load / generate training data
    # ------------------------------------------------------------------
    if args.input:
        if verbose:
            print(f"Loading SCONE binary: {args.input}")
        points, sdfs, bbox = load_scone_binary(args.input)
        bbox_min = bbox['min']
        bbox_max = bbox['max']
    else:
        # Analytic sphere
        center = np.array(args.center, dtype=np.float64)
        r = args.radius
        bbox_min = np.array(args.bbox_min, dtype=np.float64) if args.bbox_min \
                   else center - 1.5 * r
        bbox_max = np.array(args.bbox_max, dtype=np.float64) if args.bbox_max \
                   else center + 1.5 * r
        if verbose:
            print(f"Generating sphere SDF: radius={r}, center={center}, "
                  f"n_samples={args.n_samples}")
        points, sdfs = generate_sphere_sdf(
            args.n_samples, r, center, bbox_min, bbox_max, seed=args.seed)

    if verbose:
        print(f"  Samples: {len(points):,}")
        print(f"  SDF range: [{sdfs.min():.4f}, {sdfs.max():.4f}]")
        print(f"  Bbox: {bbox_min} to {bbox_max}")

    # ------------------------------------------------------------------
    # Dataset and data loaders
    # ------------------------------------------------------------------
    dataset = SdfDataset(points, sdfs, bbox_min, bbox_max)
    train_ds, val_ds = make_train_val_split(dataset, args.val_fraction, args.seed)

    # Cast dataset tensors to target dtype
    # SdfDataset produces float64; if float32 requested, re-wrap
    if args.float32:
        train_ds.dataset.coords = train_ds.dataset.coords.float()
        train_ds.dataset.sdfs   = train_ds.dataset.sdfs.float()

    train_loader = DataLoader(train_ds, batch_size=args.batch_size,
                              shuffle=True,  num_workers=0)
    val_loader   = DataLoader(val_ds,   batch_size=args.batch_size,
                              shuffle=False, num_workers=0)

    if verbose:
        print(f"  Train: {len(train_ds):,}  Val: {len(val_ds):,}")

    # ------------------------------------------------------------------
    # Model
    # ------------------------------------------------------------------
    model = NeuralSDF(hidden_dim  = args.hidden_dim,
                      num_layers  = args.num_layers,
                      activation  = args.activation,
                      leaky_alpha = args.leaky_alpha)
    if dtype == torch.float64:
        model = model.double()

    # Set sdf_scale to cover the data range (1.05× margin)
    sdf_abs_max = float(np.abs(sdfs).max())
    sdf_scale   = sdf_abs_max * 1.05
    model.set_sdf_scale(sdf_scale)

    if verbose:
        print(f"\nModel: {model.num_layers} layers, {model.hidden_dim} hidden, "
              f"activation={model.activation}, sdf_scale={sdf_scale:.4f}")
        print(f"  Parameters: {model.parameter_count():,}")

    # ------------------------------------------------------------------
    # Train
    # ------------------------------------------------------------------
    if verbose:
        print(f"\nTraining ({args.epochs} epochs, lr={args.lr})...")
    history = train(model, train_loader, val_loader,
                    epochs       = args.epochs,
                    lr           = args.lr,
                    weight_decay = args.weight_decay,
                    device       = device,
                    verbose      = verbose)

    final = history[-1]
    print(f"\nFinal  train={final['train_loss']:.6f}  val={final['val_loss']:.6f}")

    # Move model back to CPU for export
    model = model.cpu()
    if dtype == torch.float64:
        model = model.double()

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------
    export_weights(model, args.output, bbox_min, bbox_max)
    export_text(model, args.output + '.txt', bbox_min, bbox_max)

    # Save PyTorch checkpoint
    pt_path = args.output + '.pt'
    torch.save({
        'model_state_dict': model.state_dict(),
        'hidden_dim'  : model.hidden_dim,
        'num_layers'  : model.num_layers,
        'activation'  : model.activation,
        'leaky_alpha' : model.leaky_alpha,
        'sdf_scale'   : model.sdf_scale,
        'bbox_min'    : bbox_min.tolist(),
        'bbox_max'    : bbox_max.tolist(),
        'history'     : history,
    }, pt_path)
    print(f"Saved checkpoint: {pt_path}")


if __name__ == '__main__':
    main()
