"""
NeuralSDF training script.

Usage examples:

  # Train on Python-native sphere SDF (no SCONE required):
  python train.py --sphere --radius 1.0 --n-samples 100000 \\
      --output /home/jamie/SCONE/weights/sphere_neural.bin --epochs 1000

  # Train on SCONE binary output:
  python train.py --input sphere.bin --output sphere_weights.bin

  # Resume after a crash:
  python train.py --sphere --radius 1.0 --n-samples 100000 \\
      --output /home/jamie/SCONE/weights/sphere_neural.bin --epochs 1000 \\
      --resume /home/jamie/SCONE/weights/sphere_neural_ep050.pt

  # Common overrides:
  python train.py --input shape.bin --output shape.bin \\
      --hidden-dim 64 --num-layers 3 --epochs 200 --lr 1e-3

Outputs:
  <output>                              Final binary weights for SCONE
  <output>.txt                          Text weights (debugging)
  <output>.pt                           Final PyTorch checkpoint
  <SCONE_ROOT>/weights/<stem>_epNNN.bin Periodic SCONE binary checkpoints
  <SCONE_ROOT>/weights/<stem>_epNNN.pt  Resumable PyTorch checkpoints

Crash resilience:
  Status is written after EVERY epoch to:
    <SCONE_ROOT>/neural_training_status.json
  Every crash (exception or SIGTERM) is APPENDED to:
    <SCONE_ROOT>/neural_crash_log.jsonl
  SIGTERM saves an emergency checkpoint before exiting.
  SIGKILL (OOM killer) leaves the status file at the last completed epoch.
  Default checkpoint interval: 50 epochs (--checkpoint-every N).
"""

import argparse
import copy
import json
import os
import signal
import sys
import traceback
from datetime import datetime

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from model   import NeuralSDF
from sampler import (load_scone_binary, generate_sphere_sdf,
                     SdfDataset, make_train_val_split)
from export  import export_weights, export_text


_SCONE_ROOT  = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_STATUS_FILE = os.path.join(_SCONE_ROOT, 'neural_training_status.json')
_CRASH_LOG   = os.path.join(_SCONE_ROOT, 'neural_crash_log.jsonl')

# Set by SIGTERM handler; checked after every training epoch.
_SIGTERM_RECEIVED = False


class _SigtermInterrupt(Exception):
    """Raised by the per-epoch callback when SIGTERM is detected."""


def _sigterm_handler(signum, frame):
    global _SIGTERM_RECEIVED
    _SIGTERM_RECEIVED = True


def _get_mem_mb():
    """Current process RSS in MB. Reads /proc/self/status; falls back to resource module."""
    try:
        with open('/proc/self/status') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    return int(line.split()[1]) / 1024.0
    except OSError:
        pass
    try:
        import resource as _res
        ru = _res.getrusage(_res.RUSAGE_SELF)
        return ru.ru_maxrss / 1024.0  # Linux: kB -> MB
    except Exception:
        pass
    return None


def _write_status(data):
    """Atomically write training status to the project-root JSON file."""
    data['last_update'] = datetime.now().isoformat(timespec='seconds')
    tmp = _STATUS_FILE + '.tmp'
    try:
        with open(tmp, 'w') as f:
            json.dump(data, f, indent=2)
        os.replace(tmp, _STATUS_FILE)
    except OSError:
        pass  # Non-fatal


def _write_crash_log(record):
    """Append one crash record to the persistent JSONL crash log (never truncated)."""
    record.setdefault('timestamp', datetime.now().isoformat(timespec='seconds'))
    try:
        with open(_CRASH_LOG, 'a') as f:
            f.write(json.dumps(record) + '\n')
    except OSError:
        pass  # Non-fatal


# ---------------------------------------------------------------------------
# Training loop
# ---------------------------------------------------------------------------

def train(model, train_loader, val_loader, epochs, lr, weight_decay, device, verbose,
          checkpoint_every=50, on_checkpoint=None, on_epoch=None,
          history=None, start_epoch=1, eta_min=0.0, use_scheduler=True,
          clamp_delta=0.1, eikonal_weight=0.1, eikonal_scale=1.0):
    """
    Train model with Adam + clamped L1 loss + optional Eikonal regularisation.

    Args:
        checkpoint_every : save a checkpoint every this many epochs (0 = disabled)
        on_checkpoint    : callable(epoch, model, history_so_far) invoked after each
                           checkpoint epoch; history_so_far covers start_epoch..epoch
        on_epoch         : callable(epoch, model, train_loss, val_loss) called after
                           every epoch. May raise _SigtermInterrupt to abort.
        history          : list to append per-epoch records into. If None a fresh list
                           is used. Pass in a list from the caller to read it after
                           an exception (the reference stays valid through the raise).
        start_epoch      : first epoch number (default 1; higher when resuming)
        clamp_delta      : SDF clamping threshold (physical units). Loss is computed on
                           clamp(pred, -δ, δ) vs clamp(gt, -δ, δ), concentrating gradient
                           signal on the near-surface region. 0 = disabled.
        eikonal_weight   : Weight λ for Eikonal regularisation term
                           λ * mean(|∇f| - eikonal_scale)². 0 = disabled.
        eikonal_scale    : Target gradient norm in normalised coordinates.
                           For isotropic bbox of half-width h: eikonal_scale = h
                           (physical Eikonal |∇f_phys|=1 → |∇f_norm|=h).

    Returns:
        history : the same list that was appended to
    """
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    criterion = nn.L1Loss()

    # Cosine annealing: decays lr from initial value down to eta_min over the
    # full training run, eliminating the oscillation caused by constant lr Adam.
    # For resumed runs, last_epoch re-positions the scheduler correctly.
    scheduler = None
    if use_scheduler:
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, T_max=epochs, eta_min=eta_min,
            last_epoch=start_epoch - 2 if start_epoch > 1 else -1
        )

    use_eikonal = eikonal_weight > 0.0
    use_clamp   = clamp_delta > 0.0

    if history is None:
        history = []

    for epoch in range(start_epoch, epochs + 1):
        # --- Training ---
        model.train()
        train_loss = 0.0
        train_sdf_loss = 0.0
        train_eik_loss = 0.0
        for coords, sdfs in train_loader:
            sdfs = sdfs.to(device)

            # Eikonal requires input gradients; enable before forward pass
            if use_eikonal:
                coords = coords.to(device).requires_grad_(True)
            else:
                coords = coords.to(device)

            optimizer.zero_grad()
            pred = model(coords)

            # Clamped SDF loss: focus gradient signal near the surface
            if use_clamp:
                sdf_loss = criterion(torch.clamp(pred, -clamp_delta, clamp_delta),
                                     torch.clamp(sdfs,  -clamp_delta, clamp_delta))
            else:
                sdf_loss = criterion(pred, sdfs)

            # Eikonal regularisation: penalise |∇f| deviating from eikonal_scale
            if use_eikonal:
                grad = torch.autograd.grad(
                    pred.sum(), coords, create_graph=True
                )[0]
                eik_loss = (grad.norm(dim=1) - eikonal_scale).pow(2).mean()
                loss = sdf_loss + eikonal_weight * eik_loss
                train_eik_loss += eik_loss.item() * len(coords)
            else:
                loss = sdf_loss

            loss.backward()
            optimizer.step()
            train_sdf_loss += sdf_loss.item() * len(coords)
            train_loss     += loss.item()     * len(coords)

        n_train = len(train_loader.dataset)
        train_loss     /= n_train
        train_sdf_loss /= n_train
        train_eik_loss /= n_train

        # --- Validation: SDF loss only (no eikonal, no grad tracking) ---
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for coords, sdfs in val_loader:
                coords = coords.to(device)
                sdfs   = sdfs.to(device)
                pred   = model(coords)
                if use_clamp:
                    val_loss += criterion(
                        torch.clamp(pred, -clamp_delta, clamp_delta),
                        torch.clamp(sdfs,  -clamp_delta, clamp_delta)
                    ).item() * len(coords)
                else:
                    val_loss += criterion(pred, sdfs).item() * len(coords)
        val_loss /= len(val_loader.dataset)

        if scheduler:
            scheduler.step()

        history.append({'epoch': epoch, 'train_loss': train_loss,
                        'train_sdf_loss': train_sdf_loss,
                        'train_eik_loss': train_eik_loss,
                        'val_loss': val_loss})

        if verbose:
            current_lr = optimizer.param_groups[0]['lr']
            eik_str = f'  eik={train_eik_loss:.4f}' if use_eikonal else ''
            print(f"  Epoch {epoch:4d}/{epochs}  "
                  f"sdf={train_sdf_loss:.6f}  val={val_loss:.6f}"
                  f"{eik_str}  lr={current_lr:.2e}")

        # Per-epoch callback: status update, SIGTERM check.
        # history already has this epoch appended, so emergency checkpoints are complete.
        if on_epoch:
            on_epoch(epoch, model, train_loss, val_loss)

        if checkpoint_every > 0 and on_checkpoint and epoch % checkpoint_every == 0:
            on_checkpoint(epoch, model, history[:])

    return history


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Train a NeuralSDF MLP and export weights for SCONE.')

    src = p.add_mutually_exclusive_group(required=False)
    src.add_argument('--input', metavar='FILE',
                     help='SCONE sdfSampler binary output file')
    src.add_argument('--sphere', action='store_true',
                     help='Generate training data from analytic sphere SDF')

    p.add_argument('--radius', type=float, default=1.0)
    p.add_argument('--center', type=float, nargs=3, default=[0.0, 0.0, 0.0],
                   metavar=('CX', 'CY', 'CZ'))
    p.add_argument('--bbox-min', type=float, nargs=3,
                   metavar=('XMIN', 'YMIN', 'ZMIN'))
    p.add_argument('--bbox-max', type=float, nargs=3,
                   metavar=('XMAX', 'YMAX', 'ZMAX'))
    p.add_argument('--n-samples', type=int, default=100_000)

    p.add_argument('--output', required=True, metavar='FILE',
                   help='Output binary weight file path')
    p.add_argument('--resume', metavar='CHECKPOINT_PT',
                   help='Resume from a .pt checkpoint. --sphere or --input still required.')

    p.add_argument('--hidden-dim',   type=int,   default=128)
    p.add_argument('--num-layers',   type=int,   default=4)
    p.add_argument('--activation',   choices=['leakyrelu', 'relu', 'tanh'], default='leakyrelu')
    p.add_argument('--leaky-alpha',  type=float, default=0.01)

    p.add_argument('--epochs',       type=int,   default=100)
    p.add_argument('--lr',           type=float, default=1e-3)
    p.add_argument('--weight-decay', type=float, default=1e-6)
    p.add_argument('--batch-size',   type=int,   default=16384)
    p.add_argument('--no-compile',   action='store_true',
                   help='Disable torch.compile (enabled by default on PyTorch >= 2.0)')
    p.add_argument('--no-scheduler', action='store_true',
                   help='Disable cosine LR annealing (use constant lr throughout)')
    p.add_argument('--eta-min',      type=float, default=1e-5,
                   help='Minimum LR for cosine annealing (default: 1e-5)')
    p.add_argument('--clamp-delta',  type=float, default=0.1,
                   help='SDF clamping threshold in physical units (default: 0.1). '
                        'Loss is L1(clamp(pred,±δ), clamp(gt,±δ)), concentrating '
                        'gradients near the surface. Set 0 to disable.')
    p.add_argument('--eikonal-weight', type=float, default=0.1,
                   help='Weight λ for Eikonal regularisation term '
                        'λ·mean(|∇f|−target)² (default: 0.1). Set 0 to disable.')
    p.add_argument('--no-eikonal',   action='store_true',
                   help='Disable Eikonal regularisation entirely (equivalent to --eikonal-weight 0)')
    p.add_argument('--val-fraction', type=float, default=0.1)
    p.add_argument('--seed',         type=int,   default=42)
    p.add_argument('--float32',      action='store_true',
                   help='Train in float32 (default: float64)')

    p.add_argument('--checkpoint-every', type=int, default=50, metavar='N',
                   help='Save .pt + .bin checkpoint every N epochs (default 50; 0=off)')
    p.add_argument('--checkpoint-dir', metavar='DIR', default=None,
                   help='Checkpoint directory (default: <SCONE_ROOT>/weights/)')

    p.add_argument('--quiet', action='store_true')

    return p.parse_args()


def main():
    args = parse_args()
    verbose = not args.quiet
    dtype = torch.float64 if not args.float32 else torch.float32

    # Install SIGTERM handler (WSL2/systemd may send this before SIGKILL)
    signal.signal(signal.SIGTERM, _sigterm_handler)

    if not args.sphere and not args.input:
        print("error: one of --sphere or --input is required", file=sys.stderr)
        sys.exit(1)

    torch.manual_seed(args.seed)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    if verbose:
        print(f"Device: {device}", flush=True)

    # ------------------------------------------------------------------
    # Resume checkpoint
    # ------------------------------------------------------------------
    resumed_history = []
    start_epoch = 1
    resume_ckpt = None

    # Write an early status record immediately so any SIGKILL during data
    # loading or model setup is visible in neural_training_status.json.
    _write_status({
        'status':               'running',
        'phase':                'initialising',
        'command':              ' '.join(sys.argv),
        'output_file':          args.output,
        'total_epochs':         args.epochs,
        'start_epoch':          start_epoch,
        'last_epoch_completed': start_epoch - 1,
        'last_checkpoint':      args.resume,
        'last_val_loss':        None,
        'last_train_loss':      None,
        'current_ram_mb':       round(_get_mem_mb() or 0, 1),
        'peak_ram_mb':          round(_get_mem_mb() or 0, 1),
        'start_time':           datetime.now().isoformat(timespec='seconds'),
        'error':                None,
        'traceback':            None,
    })

    if args.resume:
        if verbose:
            print(f"\nLoading checkpoint: {args.resume}")
        resume_ckpt = torch.load(args.resume, map_location='cpu', weights_only=False)
        resumed_history = resume_ckpt.get('history', [])
        start_epoch = len(resumed_history) + 1

        ckpt_arch = {k: resume_ckpt[k]
                     for k in ('hidden_dim', 'num_layers', 'activation', 'leaky_alpha')}
        if (args.hidden_dim != ckpt_arch['hidden_dim'] or
                args.num_layers != ckpt_arch['num_layers'] or
                args.activation != ckpt_arch['activation']):
            print("  WARNING: architecture args differ from checkpoint — "
                  "using checkpoint architecture.", file=sys.stderr)
        args.hidden_dim  = ckpt_arch['hidden_dim']
        args.num_layers  = ckpt_arch['num_layers']
        args.activation  = ckpt_arch['activation']
        args.leaky_alpha = ckpt_arch['leaky_alpha']

        if resumed_history:
            last = resumed_history[-1]
            print(f"  Resuming: start epoch {start_epoch}/{args.epochs}, "
                  f"last val_loss={last['val_loss']:.6f}")

        if start_epoch > args.epochs:
            print(f"  Already at epoch {start_epoch - 1} (>= --epochs {args.epochs}). Done.")
            sys.exit(0)

    # ------------------------------------------------------------------
    # Training data
    # ------------------------------------------------------------------
    if args.input:
        if verbose:
            print(f"\nLoading SCONE binary: {args.input}", flush=True)
        points, sdfs, bbox = load_scone_binary(args.input)
        bbox_min = bbox['min']
        bbox_max = bbox['max']
    else:
        center = np.array(args.center, dtype=np.float64)
        r = args.radius
        bbox_min = np.array(args.bbox_min, dtype=np.float64) if args.bbox_min \
                   else center - 1.5 * r
        bbox_max = np.array(args.bbox_max, dtype=np.float64) if args.bbox_max \
                   else center + 1.5 * r
        if verbose:
            print(f"\nGenerating sphere SDF: radius={r}, center={center}, "
                  f"n_samples={args.n_samples}", flush=True)
        points, sdfs = generate_sphere_sdf(
            args.n_samples, r, center, bbox_min, bbox_max, seed=args.seed)

    # Eikonal target: physical |∇f|=1 becomes |∇f_norm|=h in normalised coords,
    # where h = mean half-width of the bounding box.
    eikonal_scale = float(np.mean((bbox_max - bbox_min) / 2.0))
    eikonal_weight = 0.0 if args.no_eikonal else args.eikonal_weight

    if verbose:
        print(f"  Samples: {len(points):,}")
        print(f"  SDF range: [{sdfs.min():.4f}, {sdfs.max():.4f}]")
        print(f"  Bbox: {bbox_min} to {bbox_max}")
        print(f"  Eikonal scale (normalised target |∇f|): {eikonal_scale:.4f}")

    # ------------------------------------------------------------------
    # Dataset and loaders
    # ------------------------------------------------------------------
    dataset = SdfDataset(points, sdfs, bbox_min, bbox_max)
    train_ds, val_ds = make_train_val_split(dataset, args.val_fraction, args.seed)

    if args.float32:
        train_ds.dataset.coords = train_ds.dataset.coords.float()
        train_ds.dataset.sdfs   = train_ds.dataset.sdfs.float()

    train_loader = DataLoader(train_ds, batch_size=args.batch_size,
                              shuffle=True,  num_workers=0, pin_memory=False)
    val_loader   = DataLoader(val_ds,   batch_size=args.batch_size,
                              shuffle=False, num_workers=0, pin_memory=False)

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

    sdf_abs_max = float(np.abs(sdfs).max())
    model.set_sdf_scale(sdf_abs_max * 1.05)

    if resume_ckpt is not None:
        model.load_state_dict(resume_ckpt['model_state_dict'])
        model.set_sdf_scale(resume_ckpt['sdf_scale'])
        if verbose:
            print(f"  Loaded weights (sdf_scale={resume_ckpt['sdf_scale']:.4f})")

    if verbose:
        print(f"\nModel: {model.num_layers} layers, {model.hidden_dim} hidden, "
              f"activation={model.activation}, sdf_scale={model.sdf_scale:.4f}")
        print(f"  Parameters: {model.parameter_count():,}")

    # ------------------------------------------------------------------
    # torch.compile (PyTorch 2.0+): compiles the model graph for faster
    # CPU execution. The compiled wrapper shares parameters with `model`,
    # so checkpointing always reads from the original `model` object.
    # ------------------------------------------------------------------
    train_model = model  # what gets passed to the training loop
    if not args.no_compile and hasattr(torch, 'compile'):
        try:
            train_model = torch.compile(model)
            if verbose:
                print("  torch.compile: enabled (use --no-compile to disable)")
        except Exception as e:
            if verbose:
                print(f"  torch.compile: skipped ({e})")

    # ------------------------------------------------------------------
    # Checkpoint directory
    # ------------------------------------------------------------------
    if args.checkpoint_dir is None:
        args.checkpoint_dir = os.path.join(_SCONE_ROOT, 'weights')
    if args.checkpoint_every > 0:
        os.makedirs(args.checkpoint_dir, exist_ok=True)

    ckpt_stem = os.path.join(
        args.checkpoint_dir,
        os.path.splitext(os.path.basename(args.output))[0]
    )

    # ------------------------------------------------------------------
    # Status dict (written after every epoch)
    # ------------------------------------------------------------------
    init_mem = _get_mem_mb()
    status = {
        'status':               'running',
        'phase':                'training',
        'command':              ' '.join(sys.argv),
        'output_file':          args.output,
        'total_epochs':         args.epochs,
        'start_epoch':          start_epoch,
        'last_epoch_completed': start_epoch - 1,
        'last_checkpoint':      args.resume,
        'last_val_loss':        resumed_history[-1]['val_loss']   if resumed_history else None,
        'last_train_loss':      resumed_history[-1]['train_loss'] if resumed_history else None,
        'current_ram_mb':       round(init_mem, 1) if init_mem else None,
        'peak_ram_mb':          round(init_mem, 1) if init_mem else None,
        'start_time':           datetime.now().isoformat(timespec='seconds'),
        'error':                None,
        'traceback':            None,
    }
    _write_status(status)

    # ------------------------------------------------------------------
    # _save_checkpoint helper (used by on_checkpoint and on_epoch/SIGTERM)
    # ------------------------------------------------------------------
    def _save_checkpoint(epoch, model_dev, hist_new, suffix=''):
        """
        Deepcopy model to CPU, write .bin + .pt checkpoint.
        hist_new is the new history only (resumed_history is prepended here).
        suffix is appended to the epoch stamp (e.g. '_emergency').
        Returns the .pt path, or None on failure.
        """
        try:
            m_cpu = copy.deepcopy(model_dev).cpu()
            full_hist = resumed_history + hist_new

            ckpt_bin = f"{ckpt_stem}_ep{epoch:03d}{suffix}.bin"
            export_weights(m_cpu, ckpt_bin, bbox_min, bbox_max)

            ckpt_pt = f"{ckpt_stem}_ep{epoch:03d}{suffix}.pt"
            torch.save({
                'model_state_dict': m_cpu.state_dict(),
                'hidden_dim':  m_cpu.hidden_dim,
                'num_layers':  m_cpu.num_layers,
                'activation':  m_cpu.activation,
                'leaky_alpha': m_cpu.leaky_alpha,
                'sdf_scale':   m_cpu.sdf_scale,
                'bbox_min':    bbox_min.tolist(),
                'bbox_max':    bbox_max.tolist(),
                'history':     full_hist,
            }, ckpt_pt)
            return ckpt_pt
        except Exception as save_err:
            print(f"  WARNING: checkpoint save failed: {save_err}", file=sys.stderr)
            return None

    # ------------------------------------------------------------------
    # Per-epoch callback: update status + handle SIGTERM
    # ------------------------------------------------------------------
    # new_history is the list passed into train(); it is populated in-place
    # by the training loop so on_epoch always sees completed epochs.
    new_history = []

    def on_epoch(epoch, model_dev, train_loss, val_loss):
        mem = _get_mem_mb()
        status['last_epoch_completed'] = epoch
        status['last_train_loss']      = train_loss
        status['last_val_loss']        = val_loss
        if mem is not None:
            status['current_ram_mb'] = round(mem, 1)
            status['peak_ram_mb']    = round(max(status.get('peak_ram_mb') or 0, mem), 1)
        _write_status(status)

        if not _SIGTERM_RECEIVED:
            return

        print(f"\nSIGTERM received after epoch {epoch}. Saving emergency checkpoint...",
              flush=True)
        ckpt_pt = _save_checkpoint(epoch, model_dev, new_history, suffix='_emergency')

        status['status'] = 'crashed'
        status['error']  = f'SIGTERM received at epoch {epoch}'
        if ckpt_pt:
            status['last_checkpoint'] = ckpt_pt
            print(f"  Saved: {ckpt_pt}", flush=True)
        else:
            print("  WARNING: emergency checkpoint save failed.", flush=True)
        _write_status(status)

        _write_crash_log({
            'crash_type':           'sigterm',
            'last_epoch_completed': epoch,
            'total_epochs':         args.epochs,
            'error':                f'SIGTERM received at epoch {epoch}',
            'traceback':            None,
            'last_val_loss':        val_loss,
            'last_train_loss':      train_loss,
            'peak_ram_mb':          status.get('peak_ram_mb'),
            'last_checkpoint':      status.get('last_checkpoint'),
            'command':              ' '.join(sys.argv),
        })
        raise _SigtermInterrupt(f'SIGTERM at epoch {epoch}')

    # ------------------------------------------------------------------
    # Periodic checkpoint callback
    # ------------------------------------------------------------------
    def on_checkpoint(epoch, model_dev, history_so_far):
        ckpt_pt = _save_checkpoint(epoch, model_dev, history_so_far)
        if ckpt_pt:
            status['last_checkpoint'] = ckpt_pt
        _write_status(status)
        if verbose:
            print(f"  Checkpoint ep{epoch:03d}: "
                  f"{ckpt_pt if ckpt_pt else '(save failed)'}")

    if args.checkpoint_every == 0:
        on_checkpoint_cb = None
    else:
        on_checkpoint_cb = on_checkpoint

    # ------------------------------------------------------------------
    # Train
    # ------------------------------------------------------------------
    if verbose:
        print(f"\nTraining epochs {start_epoch}–{args.epochs}  "
              f"(lr={args.lr}, wd={args.weight_decay}, "
              f"checkpoint every {args.checkpoint_every} epochs)...")

    try:
        history = train(train_model, train_loader, val_loader,
                        epochs           = args.epochs,
                        lr               = args.lr,
                        weight_decay     = args.weight_decay,
                        device           = device,
                        verbose          = verbose,
                        checkpoint_every = args.checkpoint_every,
                        on_checkpoint    = on_checkpoint_cb,
                        on_epoch         = on_epoch,
                        history          = new_history,
                        start_epoch      = start_epoch,
                        eta_min          = args.eta_min,
                        use_scheduler    = not args.no_scheduler,
                        clamp_delta      = args.clamp_delta,
                        eikonal_weight   = eikonal_weight,
                        eikonal_scale    = eikonal_scale)

    except _SigtermInterrupt:
        # Emergency checkpoint already saved inside on_epoch.
        sys.exit(0)

    except KeyboardInterrupt:
        status['status'] = 'crashed'
        status['error']  = 'KeyboardInterrupt'
        _write_status(status)
        _write_crash_log({
            'crash_type':           'keyboard_interrupt',
            'last_epoch_completed': status['last_epoch_completed'],
            'total_epochs':         args.epochs,
            'error':                'KeyboardInterrupt',
            'traceback':            None,
            'last_val_loss':        status.get('last_val_loss'),
            'last_train_loss':      status.get('last_train_loss'),
            'peak_ram_mb':          status.get('peak_ram_mb'),
            'last_checkpoint':      status.get('last_checkpoint'),
            'command':              ' '.join(sys.argv),
        })
        print('\nInterrupted. Status and crash log updated.', flush=True)
        raise

    except Exception as e:
        tb = traceback.format_exc()
        status['status']    = 'crashed'
        status['error']     = f"{type(e).__name__}: {e}"
        status['traceback'] = tb
        _write_status(status)
        _write_crash_log({
            'crash_type':           'exception',
            'last_epoch_completed': status['last_epoch_completed'],
            'total_epochs':         args.epochs,
            'error':                f"{type(e).__name__}: {e}",
            'traceback':            tb,
            'last_val_loss':        status.get('last_val_loss'),
            'last_train_loss':      status.get('last_train_loss'),
            'peak_ram_mb':          status.get('peak_ram_mb'),
            'last_checkpoint':      status.get('last_checkpoint'),
            'command':              ' '.join(sys.argv),
        })
        print(f"\nCrashed: {status['error']}", flush=True)
        print(f"Crash appended to {_CRASH_LOG}", flush=True)
        raise

    # ------------------------------------------------------------------
    # Export final weights
    # ------------------------------------------------------------------
    final = history[-1]
    print(f"\nFinal  train={final['train_loss']:.6f}  val={final['val_loss']:.6f}")

    model = model.cpu()
    if dtype == torch.float64:
        model = model.double()

    export_weights(model, args.output, bbox_min, bbox_max)
    export_text(model, args.output + '.txt', bbox_min, bbox_max)

    full_history = resumed_history + history
    pt_path = args.output + '.pt'
    torch.save({
        'model_state_dict': model.state_dict(),
        'hidden_dim':  model.hidden_dim,
        'num_layers':  model.num_layers,
        'activation':  model.activation,
        'leaky_alpha': model.leaky_alpha,
        'sdf_scale':   model.sdf_scale,
        'bbox_min':    bbox_min.tolist(),
        'bbox_max':    bbox_max.tolist(),
        'history':     full_history,
    }, pt_path)
    print(f"Saved final checkpoint: {pt_path}")

    status['status']               = 'completed'
    status['last_epoch_completed'] = args.epochs
    status['last_val_loss']        = final['val_loss']
    status['last_checkpoint']      = pt_path
    _write_status(status)


if __name__ == '__main__':
    main()
