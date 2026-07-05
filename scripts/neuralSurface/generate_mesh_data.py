"""
Generate binary inside/outside training data from a mesh file and save to disk.

Output is written in the SCONE binary format so it can be loaded with
  train.py --input <output_file> --mode binary --loss bce ...

Unit convention: SCONE works in cm. PLY files are often in metres.
Use --scale 100 to convert a metre-scale mesh to cm before sampling so that
the bbox stored in the weights file matches SCONE coordinate units directly.

Usage:
  python generate_mesh_data.py --mesh meshes/bunny_watertight.ply \
      --n-samples 200000 --scale 100 --output bunny_train.bin
"""

import argparse
import struct
import sys
import time
import numpy as np

# Allow running from the repo root or from this directory
import os
sys.path.insert(0, os.path.dirname(__file__))
from sampler import generate_mesh_binary


def save_scone_binary(filename, points, labels, bbox_min, bbox_max):
    """
    Write points + labels (used as the SDF field) in SCONE binary format.

    Format:
      [int32]        count
      [float64 × 6]  bbox: xmin ymin zmin xmax ymax zmax
      [float64 × 4]  × count: x y z label
    """
    count = len(points)
    data = np.column_stack([points, labels]).astype(np.float64)
    bbox_flat = np.concatenate([bbox_min, bbox_max]).astype(np.float64)

    with open(filename, 'wb') as f:
        f.write(struct.pack('<i', count))
        f.write(bbox_flat.tobytes())
        f.write(data.tobytes())


def main():
    p = argparse.ArgumentParser(description='Generate mesh training data for NeuralSDF.')
    p.add_argument('--mesh',     required=True, metavar='FILE',
                   help='Mesh file (PLY/OBJ/STL); must be watertight for reliable labels')
    p.add_argument('--n-samples', type=int, default=200_000, metavar='N')
    p.add_argument('--output',   required=True, metavar='FILE',
                   help='Output binary file (SCONE format, loadable with --input)')
    p.add_argument('--near-fraction',     type=float, default=0.65,
                   help='Fraction of samples placed near the surface (default 0.65)')
    p.add_argument('--near-distance-rel', type=float, default=0.05,
                   help='Near-surface band as fraction of bbox diagonal (default 0.05)')
    p.add_argument('--bbox-padding',      type=float, default=0.1,
                   help='Fractional padding added to mesh bounds (default 0.1)')
    p.add_argument('--scale', type=float, default=1.0, metavar='FACTOR',
                   help='Multiply all mesh coordinates by this factor before sampling '
                        '(e.g. --scale 100 converts metres to cm). Default 1.0 (no scaling).')
    p.add_argument('--seed', type=int, default=42)
    args = p.parse_args()

    print(f"Mesh:      {args.mesh}", flush=True)
    print(f"Samples:   {args.n_samples:,}", flush=True)
    print(f"Scale:     {args.scale}x", flush=True)
    print(f"Output:    {args.output}", flush=True)

    t0 = time.time()
    print("\nLoading mesh and sampling points...", flush=True)
    points, labels, bbox_min, bbox_max = generate_mesh_binary(
        args.mesh,
        args.n_samples,
        near_fraction=args.near_fraction,
        near_distance_rel=args.near_distance_rel,
        bbox_padding=args.bbox_padding,
        seed=args.seed,
        scale=args.scale,
    )
    t1 = time.time()

    n_inside  = int((labels < 0).sum())
    n_outside = int((labels > 0).sum())
    print(f"Done in {t1 - t0:.1f}s  "
          f"({n_inside:,} inside / {n_outside:,} outside, "
          f"ratio {n_inside / len(labels):.3f})", flush=True)
    print(f"BBox min: {bbox_min}", flush=True)
    print(f"BBox max: {bbox_max}", flush=True)

    save_scone_binary(args.output, points, labels, bbox_min, bbox_max)
    size_mb = os.path.getsize(args.output) / 1e6
    print(f"\nWritten {args.output}  ({size_mb:.1f} MB)", flush=True)


if __name__ == '__main__':
    main()
