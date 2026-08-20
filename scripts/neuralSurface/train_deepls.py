"""
DeepLS training, paper-faithful architecture (Chabra et al., ECCV 2020,
arXiv:2003.10983) — a single shared decoder network conditioned on a
per-voxel latent code, trained jointly (auto-decoder style) across all active
voxels' data at once.

An earlier "independent small MLP per voxel" attempt, with no shared decoder
or latent codes at all, is not used: on the sphere validation it
underperformed the global MLP substantially (Δk −253 pcm, 0.60% vs 0.070%
misclass), and closer reading of the paper showed it was never really
DeepLS — the paper shares one full-size decoder across all voxels via
latent codes, it does not train independent networks.

Architecture, matched to the paper's stated choices (arXiv:2003.10983 App. B):
  - Shared decoder: 4 layers, 128 hidden (same size as our global MLP/DeepSDF)
  - Latent code: 125-dim (their own ablation-chosen size) -> decoder input
    dim = 125 + 3 = 128, matching their "128 input neurons" exactly
  - Adam, initial LR 0.01, decayed twice over training (step decay here)
  - Auto-decoder: latent codes are free parameters optimised jointly with the
    decoder weights via backprop, no separate encoder network

Deliberate deviation from the paper: loss stays binary BCE on sign(sdf), not
their truncated-SDF regression — this project's own Phase 1 result (see
Neural.md/project_phase1_complete) already established BCE beats SDF
regression for SCONE's actual halfspace-query use case, and changing the loss
at the same time as the architecture would confound which change caused any
result difference. Everything else follows the paper.

Per-voxel local coordinate normalisation uses an extended 1.5x receptive
field bbox, so neighbouring voxels' training data overlap slightly and avoid
a hard seam at voxel boundaries.

Usage:
  python train_deepls.py --input ../../data/rbsphere_unit_train.bin \\
      --output ../../weights/rbsphere_deepls_unit.bin \\
      --nvox 8 8 8 --latent-dim 125 --hidden-dim 128 --num-layers 4 \\
      --epochs 500 --device cpu
"""

import argparse
import copy
import json
import os
import time
import numpy as np
import torch
import torch.nn as nn

from model import NeuralSDF, _ACTIVATION_NAMES
from sampler import load_scone_binary
from deepls_common import classify_and_gather
from export_deepls import export_deepls

_SCONE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
_STATUS_FILE = os.path.join(_SCONE_ROOT, 'deepls_training_status.json')


def _write_status(data):
    tmp = _STATUS_FILE + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, _STATUS_FILE)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True)
    p.add_argument('--output', required=True)

    p.add_argument('--nvox', type=int, nargs=3, default=[8, 8, 8], metavar=('NX', 'NY', 'NZ'))
    p.add_argument('--receptive-field', type=float, default=1.5)
    p.add_argument('--min-points', type=int, default=200)

    p.add_argument('--latent-dim', type=int, default=125,
                   help='Per-voxel latent code size (default 125, the paper\'s chosen value — '
                        'combined with 3 xyz dims this gives a 128-wide decoder input, matching '
                        'the paper\'s "128 input neurons" exactly)')
    p.add_argument('--hidden-dim', type=int, default=128,
                   help='Shared decoder hidden width (default 128, matching the paper and our '
                        'own global MLP)')
    p.add_argument('--num-layers', type=int, default=4,
                   help='Shared decoder layer count (default 4, matching the paper/DeepSDF)')
    p.add_argument('--activation', choices=['leakyrelu', 'relu', 'tanh'], default='leakyrelu')
    p.add_argument('--leaky-alpha', type=float, default=0.01)

    p.add_argument('--epochs', type=int, default=500)
    p.add_argument('--lr', type=float, default=0.01,
                   help='Initial learning rate (default 0.01, matching the paper)')
    p.add_argument('--lr-decay-epochs', type=float, nargs=2, default=[0.5, 0.75],
                   help='Fractions of --epochs at which LR is multiplied by --lr-decay-factor '
                        '(default 0.5 0.75, matching the paper\'s "decreased twice")')
    p.add_argument('--lr-decay-factor', type=float, default=0.1)
    p.add_argument('--latent-reg', type=float, default=1e-4,
                   help='L2 regularisation weight on latent codes (standard auto-decoder '
                        'practice, keeps the latent space compact). 0 = disabled.')
    p.add_argument('--batch-size', type=int, default=16384)
    p.add_argument('--checkpoint-every', type=int, default=100, metavar='N',
                   help='Export a checkpoint every N epochs (0 = disabled, default 100). '
                        'Deep-copies decoder+latents to CPU before export so the live GPU '
                        'training state is never touched -- unlike the final export, this '
                        'must not call decoder.cpu() directly (nn.Module.cpu() mutates '
                        'parameters in place, which would break the training loop mid-run).')
    p.add_argument('--val-fraction', type=float, default=0.1)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--device', default='auto', choices=['auto', 'cpu', 'cuda'])
    p.add_argument('--float32', action='store_true',
                   help='Train in float32 instead of float64 (default float64, matching '
                        'train.py\'s convention). GPU FP64 throughput is often much lower than '
                        'FP32 on consumer/typical cards -- try this if a GPU run looks slow. '
                        'Weights are always exported as float64 regardless (export upcasts), '
                        'so this only affects training speed/precision, not the output format.')
    p.add_argument('--quiet', action='store_true')
    return p.parse_args()


def main():
    args = parse_args()
    verbose = not args.quiet
    nvox = np.array(args.nvox, dtype=np.int64)
    dtype = torch.float32 if args.float32 else torch.float64

    if args.device == 'auto':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device(args.device)
    torch.manual_seed(args.seed)

    if verbose:
        print(f"Device: {device}")
        print(f"Loading: {args.input}")
    points, sdfs, bbox = load_scone_binary(args.input)
    bbox_min = bbox['min']
    bbox_max = bbox['max']

    grid_origin = bbox_min.copy()
    voxel_size = (bbox_max - bbox_min) / nvox
    if verbose:
        print(f"  {len(points):,} samples, bbox {bbox_min} .. {bbox_max}")
        print(f"  Voxel grid: {nvox[0]}x{nvox[1]}x{nvox[2]}, voxel size {voxel_size}")

    # Fast, vectorised classification + extended-region gather (see
    # deepls_common.py -- the original per-voxel-loop version took minutes on
    # the teapot's 32^3 grid; this is seconds, validated identical output on
    # the sphere's 8^3 grid).
    status_grid, active_voxels, gathered, n_fallback_classify, n_fallback_gather = \
        classify_and_gather(points, sdfs, bbox_min, bbox_max, nvox,
                            receptive_field=args.receptive_field,
                            min_points=args.min_points, verbose=verbose)
    n_active = len(active_voxels)

    # ------------------------------------------------------------------
    # Phase 2: build ONE joint dataset spanning every active voxel's
    # extended-receptive-field points, tagged with a latent index per point.
    # ------------------------------------------------------------------
    voxel_bboxes = {}
    all_xyz_norm = []
    all_labels = []
    all_vidx = []

    for vidx, (ix, iy, iz) in enumerate(active_voxels):
        local_points, local_labels_raw, ext_lo, ext_hi = gathered[(ix, iy, iz)]

        voxel_bboxes[(ix, iy, iz)] = (ext_lo, ext_hi)
        xyz_norm = 2.0 * (local_points - ext_lo) / (ext_hi - ext_lo) - 1.0
        all_xyz_norm.append(xyz_norm)
        all_labels.append(local_labels_raw)
        all_vidx.append(np.full(len(local_points), vidx, dtype=np.int64))

    xyz_norm = np.concatenate(all_xyz_norm, axis=0)
    sdf_labels = np.concatenate(all_labels, axis=0)
    vidx_arr = np.concatenate(all_vidx, axis=0)
    n_points = len(xyz_norm)
    if verbose:
        print(f"Joint training set: {n_points:,} points across {n_active} active voxels "
              f"(avg {n_points/max(1,n_active):.0f} pts/voxel)")

    # ------------------------------------------------------------------
    # Model: shared decoder (NeuralSDF, in_dim = latent_dim + 3) + a free
    # latent-code table, one row per active voxel. Auto-decoder: both are
    # optimised jointly, no separate encoder.
    # ------------------------------------------------------------------
    decoder = NeuralSDF(hidden_dim=args.hidden_dim, num_layers=args.num_layers,
                        activation=args.activation, leaky_alpha=args.leaky_alpha,
                        in_dim=args.latent_dim + 3)
    decoder = (decoder.double() if dtype == torch.float64 else decoder.float()).to(device)
    decoder.set_sdf_scale(1.0)

    latents = nn.Parameter(torch.randn(n_active, args.latent_dim, dtype=dtype,
                                       device=device) * 0.01)

    # Pin the WHOLE joint dataset to `device` once, rather than the more
    # obvious DataLoader-over-CPU-tensors approach -- that would re-transfer
    # every batch, every epoch (500 epochs = 500x the transfer traffic for no
    # reason, since this dataset is easily small enough to fit in GPU memory
    # entirely). Manual index-based shuffling below replaces DataLoader.
    xyz_t = torch.tensor(xyz_norm, dtype=dtype, device=device)
    sdf_t = torch.tensor(sdf_labels, dtype=dtype, device=device)
    vidx_t = torch.tensor(vidx_arr, dtype=torch.long, device=device)

    n_val = max(1, int(n_points * args.val_fraction))
    n_train = n_points - n_val
    perm = torch.randperm(n_points, generator=torch.Generator().manual_seed(args.seed))
    train_idx = perm[:n_train].to(device)
    val_idx = perm[n_train:].to(device)

    optimizer = torch.optim.Adam(list(decoder.parameters()) + [latents], lr=args.lr)
    decay_epochs = {int(f * args.epochs) for f in args.lr_decay_epochs}

    if verbose:
        n_params = sum(p.numel() for p in decoder.parameters())
        print(f"Shared decoder: {args.num_layers} layers, {args.hidden_dim} hidden, "
              f"{n_params:,} params. Latent table: {n_active} x {args.latent_dim} "
              f"= {n_active*args.latent_dim:,} params. dtype={dtype}")

    def _batches(idx_tensor, batch_size, shuffle):
        if shuffle:
            idx_tensor = idx_tensor[torch.randperm(len(idx_tensor), device=idx_tensor.device)]
        for start in range(0, len(idx_tensor), batch_size):
            yield idx_tensor[start:start + batch_size]

    t_start = time.time()
    history = []
    for epoch in range(1, args.epochs + 1):
        if epoch in decay_epochs:
            for g in optimizer.param_groups:
                g['lr'] *= args.lr_decay_factor
            if verbose:
                print(f"  LR decayed to {optimizer.param_groups[0]['lr']:.2e} at epoch {epoch}")

        decoder.train()
        train_loss_sum = torch.zeros((), dtype=dtype, device=device)
        for b_idx in _batches(train_idx, args.batch_size, shuffle=True):
            xyz_b, label_b, vidx_b = xyz_t[b_idx], sdf_t[b_idx], vidx_t[b_idx]
            lat_b = latents[vidx_b]
            decoder_in = torch.cat([lat_b, xyz_b], dim=1)
            pred = decoder(decoder_in).squeeze(-1)

            p_prob = (pred + 1) / 2
            t_prob = (label_b.sign() + 1) / 2
            loss = torch.nn.functional.binary_cross_entropy(p_prob.clamp(1e-7, 1 - 1e-7), t_prob)
            if args.latent_reg > 0:
                loss = loss + args.latent_reg * lat_b.pow(2).sum(-1).mean()

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            # Accumulate as a GPU tensor -- one .item() sync per epoch below,
            # not one per batch (each .item() call forces a CPU-GPU sync that
            # stalls the pipeline instead of letting kernels queue ahead).
            train_loss_sum += loss.detach() * len(b_idx)
        train_loss = (train_loss_sum / n_train).item()

        decoder.eval()
        val_loss_sum = torch.zeros((), dtype=dtype, device=device)
        with torch.no_grad():
            for b_idx in _batches(val_idx, args.batch_size, shuffle=False):
                xyz_b, label_b, vidx_b = xyz_t[b_idx], sdf_t[b_idx], vidx_t[b_idx]
                lat_b = latents[vidx_b]
                decoder_in = torch.cat([lat_b, xyz_b], dim=1)
                pred = decoder(decoder_in).squeeze(-1)
                p_prob = (pred + 1) / 2
                t_prob = (label_b.sign() + 1) / 2
                vloss = torch.nn.functional.binary_cross_entropy(
                    p_prob.clamp(1e-7, 1 - 1e-7), t_prob)
                val_loss_sum += vloss.detach() * len(b_idx)
        val_loss = (val_loss_sum / n_val).item()

        history.append({'epoch': epoch, 'train_loss': train_loss, 'val_loss': val_loss})
        elapsed = time.time() - t_start
        eta = (elapsed / epoch) * (args.epochs - epoch)
        if verbose and (epoch % 10 == 0 or epoch == 1 or epoch == args.epochs):
            print(f"  epoch {epoch}/{args.epochs}  train={train_loss:.5f}  val={val_loss:.5f}  "
                  f"ETA {eta/60:.1f} min", flush=True)
        _write_status({
            'status': 'training', 'epoch': epoch, 'total_epochs': args.epochs,
            'train_loss': train_loss, 'val_loss': val_loss,
            'elapsed_seconds': round(elapsed, 1), 'eta_seconds': round(eta, 1),
        })

        if args.checkpoint_every > 0 and epoch % args.checkpoint_every == 0:
            ckpt_path = f"{args.output.rsplit('.', 1)[0]}_ep{epoch:04d}.bin"
            try:
                # Deep-copy BEFORE .cpu() -- decoder.cpu() on the live object
                # would mutate its parameters in place and pull the actual
                # training-loop decoder off the GPU mid-run. The copy is
                # what gets moved/exported; `decoder`/`latents` above are
                # never touched.
                decoder_ckpt = copy.deepcopy(decoder).cpu()
                latents_ckpt = latents.detach().cpu().numpy()
                export_deepls(ckpt_path, status_grid, grid_origin, voxel_size,
                                     decoder_ckpt, latents_ckpt, active_voxels, voxel_bboxes,
                                     latent_dim=args.latent_dim, hidden_dim=args.hidden_dim,
                                     num_layers=args.num_layers,
                                     activation_type=_ACTIVATION_NAMES[args.activation],
                                     leaky_alpha=args.leaky_alpha)
                if verbose:
                    print(f"  checkpoint saved: {ckpt_path}", flush=True)
            except Exception as ckpt_err:
                # A failed checkpoint should never abort the training run --
                # print and continue, don't let a disk/IO hiccup lose
                # everything.
                print(f"  WARNING: checkpoint save failed at epoch {epoch}: {ckpt_err}",
                      flush=True)

    _write_status({'status': 'exporting', 'epoch': args.epochs, 'total_epochs': args.epochs})

    export_deepls(args.output, status_grid, grid_origin, voxel_size,
                         decoder.cpu(), latents.detach().cpu().numpy(),
                         active_voxels, voxel_bboxes,
                         latent_dim=args.latent_dim, hidden_dim=args.hidden_dim,
                         num_layers=args.num_layers,
                         activation_type=_ACTIVATION_NAMES[args.activation],
                         leaky_alpha=args.leaky_alpha)

    _write_status({'status': 'completed', 'epoch': args.epochs, 'total_epochs': args.epochs,
                  'elapsed_seconds': round(time.time() - t_start, 1),
                  'output_file': args.output})


if __name__ == '__main__':
    main()
