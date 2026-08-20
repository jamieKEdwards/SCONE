"""
Shared, vectorised voxel classification + extended-region point gathering,
used by train_deepls.py.

Added 2026-08-18 after the original per-voxel-loop implementation (nested
Python loop over every voxel, each doing a full boolean mask over the WHOLE
point array) proved fine at the sphere's scale (512 voxels x 200k points) but
took minutes just to CLASSIFY (no training) the teapot's proposed 32x32x32
grid (32,768 voxels x 2M points) -- 640x more total work than the sphere run
that took seconds. This is a straight algorithmic fix, not a data problem:
data/teapot_train.bin's format/content is fine (verified separately).

Approach: single-pass spatial binning instead of O(n_voxels * n_points):
  1. Compute each point's voxel index once (vectorised).
  2. Group point indices by voxel via one sort (classic numpy "groupby").
  3. Classify every voxel via np.bincount of positive/negative label counts
     per voxel -- O(N) total, not O(N) PER voxel.
  4. For each active voxel's extended (1.5x) receptive field: since the
     extension is only 0.25*voxel_size past each face, it can only ever reach
     into the immediately adjacent voxel -- so gather candidates from the
     3x3x3 neighbour block's pre-built point lists (small) and apply the
     exact box filter to just that candidate set, instead of the full 2M
     points.

Correctness: validated against the original brute-force implementation on
the sphere data (data/rbsphere_unit_train.bin, 8x8x8 grid) -- identical
status_grid and identical (order-independent) per-voxel point sets.
"""

import numpy as np

STATUS_CONST_INSIDE  = 0
STATUS_CONST_OUTSIDE = 1
STATUS_HAS_MLP       = 2


def classify_and_gather(points, sdfs, bbox_min, bbox_max, nvox,
                        receptive_field=1.5, min_points=200, verbose=True):
    """
    Args:
        points, sdfs : full training set, ndarray (N,3) / (N,)
        bbox_min, bbox_max : world-space bounding box, ndarray (3,)
        nvox : voxel grid dimensions, ndarray (3,) int
        receptive_field : extended receptive field multiple of voxel size
        min_points : minimum points required in a voxel's extended region;
                     falls back to the N nearest points overall if short

    Returns:
        status_grid  : ndarray (nx,ny,nz) int32, STATUS_* codes
        active_voxels: list of (ix,iy,iz) 0-based, in Fortran column-major
                       scan order (ix fastest) -- matches the binary export
                       format's required traversal order
        gathered     : dict {(ix,iy,iz): (local_points, local_sdfs, bbox_lo, bbox_hi)}
                       for every active voxel
        n_fallback_classify : voxels with an empty core (grid finer than
                       sample density), classified via nearest single point
        n_fallback_gather   : active voxels whose 3x3x3-neighbour candidate
                       set was still short of min_points, widened via a
                       global nearest-N fallback (rare unless the grid is
                       very fine relative to sample density)
    """
    nvox = np.asarray(nvox, dtype=np.int64)
    bbox_min = np.asarray(bbox_min, dtype=np.float64)
    bbox_max = np.asarray(bbox_max, dtype=np.float64)
    voxel_size = (bbox_max - bbox_min) / nvox
    labels = np.sign(sdfs)
    labels[labels == 0] = 1.0

    nx, ny, nz = (int(v) for v in nvox)
    n_total = nx * ny * nz

    # ------------------------------------------------------------------
    # Step 1-2: bin every point into its voxel, group by voxel via one sort.
    # ------------------------------------------------------------------
    idx3 = np.clip(np.floor((points - bbox_min) / voxel_size).astype(np.int64),
                   0, nvox - 1)
    lin_idx = idx3[:, 0] + idx3[:, 1] * nx + idx3[:, 2] * nx * ny  # ix fastest

    order = np.argsort(lin_idx, kind='stable')
    sorted_lin = lin_idx[order]
    # Boundaries of each voxel's block within `order`
    boundaries = np.searchsorted(sorted_lin, np.arange(n_total + 1))

    def point_indices_in_voxel(v):
        lo, hi = boundaries[v], boundaries[v + 1]
        return order[lo:hi]

    # ------------------------------------------------------------------
    # Step 3: classify via bincount (O(N) total).
    # ------------------------------------------------------------------
    pos_mask = (labels > 0).astype(np.int64)
    neg_mask = (labels < 0).astype(np.int64)
    n_pos = np.bincount(lin_idx, weights=pos_mask, minlength=n_total).astype(np.int64)
    n_neg = np.bincount(lin_idx, weights=neg_mask, minlength=n_total).astype(np.int64)

    status_flat = np.full(n_total, STATUS_HAS_MLP, dtype=np.int32)
    status_flat[(n_pos > 0) & (n_neg == 0)] = STATUS_CONST_OUTSIDE
    status_flat[(n_neg > 0) & (n_pos == 0)] = STATUS_CONST_INSIDE
    empty_core = (n_pos == 0) & (n_neg == 0)

    n_fallback_classify = int(empty_core.sum())
    if n_fallback_classify > 0:
        # Grid finer than sample density in a few spots: fall back to the
        # single nearest point's sign (rare, done unvectorised since it's
        # a handful of voxels at most).
        empty_lin = np.nonzero(empty_core)[0]
        for v in empty_lin:
            ix = v % nx
            iy = (v // nx) % ny
            iz = v // (nx * ny)
            center = bbox_min + (np.array([ix, iy, iz]) + 0.5) * voxel_size
            nearest = np.argmin(np.sum((points - center) ** 2, axis=1))
            status_flat[v] = STATUS_CONST_INSIDE if labels[nearest] < 0 else STATUS_CONST_OUTSIDE

    status_grid = status_flat.reshape(nz, ny, nx).transpose(2, 1, 0)  # -> (nx,ny,nz)

    active_voxels = [(ix, iy, iz)
                     for iz in range(nz) for iy in range(ny) for ix in range(nx)
                     if status_grid[ix, iy, iz] == STATUS_HAS_MLP]

    if verbose:
        n_in = int((status_flat == STATUS_CONST_INSIDE).sum())
        n_out = int((status_flat == STATUS_CONST_OUTSIDE).sum())
        print(f"Classification: {n_in} constant-inside, {n_out} constant-outside, "
              f"{len(active_voxels)} active (of {n_total} total; "
              f"{n_fallback_classify} used nearest-point fallback)")

    # ------------------------------------------------------------------
    # Step 4: extended-region gather via 3x3x3 neighbour blocks.
    # ------------------------------------------------------------------
    gathered = {}
    n_fallback_gather = 0
    for (ix, iy, iz) in active_voxels:
        center = bbox_min + (np.array([ix, iy, iz]) + 0.5) * voxel_size
        ext_half = 0.5 * receptive_field * voxel_size
        lo, hi = center - ext_half, center + ext_half

        candidate_lists = []
        for dz in (-1, 0, 1):
            jz = iz + dz
            if jz < 0 or jz >= nz:
                continue
            for dy in (-1, 0, 1):
                jy = iy + dy
                if jy < 0 or jy >= ny:
                    continue
                for dx in (-1, 0, 1):
                    jx = ix + dx
                    if jx < 0 or jx >= nx:
                        continue
                    v = jx + jy * nx + jz * nx * ny
                    candidate_lists.append(point_indices_in_voxel(v))
        candidates = np.concatenate(candidate_lists) if candidate_lists else np.array([], dtype=np.int64)
        cand_points = points[candidates]
        box_mask = np.all((cand_points >= lo) & (cand_points <= hi), axis=1)
        sel = candidates[box_mask]

        if len(sel) < min_points:
            n_fallback_gather += 1
            d2 = np.sum((points - center) ** 2, axis=1)
            sel = np.argsort(d2)[:min_points]

        gathered[(ix, iy, iz)] = (points[sel], sdfs[sel], lo, hi)

    if verbose and n_fallback_gather > 0:
        print(f"  ({n_fallback_gather}/{len(active_voxels)} active voxels used the "
              f"global nearest-N fallback for their extended region -- grid may be "
              f"fine relative to sample density)")

    return status_grid, active_voxels, gathered, n_fallback_classify, n_fallback_gather
