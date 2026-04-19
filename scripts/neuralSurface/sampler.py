"""
SDF training data utilities.

Supports two data sources:
  1. SCONE binary output from the --sample-sdf command (primary workflow)
  2. Python-native analytic SDF for simple primitives (proof-of-concept / testing)

The mixed sampling strategy follows DeepLS (Chabra et al., ECCV 2020, Sec. 4.3):
  32.5% tight near-surface — uniform samples filtered to |sdf| < 0.01R
  32.5% loose near-surface — uniform samples filtered to |sdf| < 0.1R
  25%   volumetric points  — uniform samples throughout bbox
  10%   on-surface points  — surface samples perturbed along surface normal

SCONE binary format (written by sdfSampler Fortran utility):
  [Header]
    count   int32         number of sample records
    bbox    float64 × 6   xmin, ymin, zmin, xmax, ymax, zmax
  [Records, count times]
    x, y, z, sdf   float64 × 4
"""

import numpy as np
import torch
from torch.utils.data import Dataset, random_split


# ---------------------------------------------------------------------------
# SCONE binary file loader
# ---------------------------------------------------------------------------

def load_scone_binary(filename):
    """
    Load SDF samples from a SCONE sdfSampler binary output file.

    Returns:
        points : ndarray of shape (N, 3), world-space coordinates
        sdfs   : ndarray of shape (N,),   signed distance values
        bbox   : dict with keys 'min' (shape 3) and 'max' (shape 3)
    """
    import struct

    with open(filename, 'rb') as f:
        count = struct.unpack('<i', f.read(4))[0]
        bbox_flat = np.frombuffer(f.read(48), dtype=np.float64)  # 6 × 8 bytes
        data = np.frombuffer(f.read(count * 4 * 8), dtype=np.float64)

    data = data.reshape(count, 4)
    points = data[:, :3]
    sdfs   = data[:,  3]
    bbox   = {'min': bbox_flat[:3], 'max': bbox_flat[3:]}
    return points, sdfs, bbox


# ---------------------------------------------------------------------------
# Python-native analytic SDF generators (for testing without SCONE)
# ---------------------------------------------------------------------------

def _rejection_sample_near_surface(rng, n, bbox_min, bbox_max, center, radius, distance):
    """Return n points uniformly from bbox with |sdf| < distance (rejection sampling)."""
    collected = []
    n_collected = 0
    oversample = max(4, int(np.ceil((np.prod(bbox_max - bbox_min)) /
                                    (4 * np.pi * radius**2 * 2 * distance))))
    batch_size = max(n * oversample, 4096)
    while n_collected < n:
        batch = rng.uniform(bbox_min, bbox_max, size=(batch_size, 3))
        sdf_b = np.linalg.norm(batch - center, axis=1) - radius
        mask = np.abs(sdf_b) < distance
        accepted = np.column_stack([batch[mask], sdf_b[mask]])
        collected.append(accepted)
        n_collected += len(accepted)
    arr = np.vstack(collected)[:n]
    return arr[:, :3], arr[:, 3]


def generate_sphere_sdf(n_samples, radius, center, bbox_min, bbox_max,
                        near_fraction=0.65, near_distance=None,
                        near_distance_tight=None, seed=42):
    """
    Generate SDF training samples for a sphere using analytic formula.
    SDF(r) = |r - center| - radius   (negative inside, positive outside)

    Near-surface samples are split equally between a tight band (near_distance_tight)
    and a loose band (near_distance), giving the network good coverage at both
    fine and coarse scales near the boundary.

    Args:
        n_samples          : Total number of samples to generate
        radius             : Sphere radius
        center             : Array-like (3,) — sphere centre
        bbox_min           : Array-like (3,) — sampling bounding box minimum
        bbox_max           : Array-like (3,) — sampling bounding box maximum
        near_fraction      : Fraction of samples concentrated near the surface
                             (split equally between tight and loose bands).
                             Default 0.65.
        near_distance      : Loose near-surface band: keep |sdf| < near_distance.
                             Defaults to 0.1 * radius if None.
        near_distance_tight: Tight near-surface band: keep |sdf| < near_distance_tight.
                             Defaults to 0.01 * radius if None.
        seed               : Random seed

    Returns:
        points : ndarray (N, 3)
        sdfs   : ndarray (N,)
    """
    rng = np.random.default_rng(seed)
    center   = np.asarray(center,   dtype=np.float64)
    bbox_min = np.asarray(bbox_min, dtype=np.float64)
    bbox_max = np.asarray(bbox_max, dtype=np.float64)

    if near_distance is None:
        near_distance = 0.1 * radius
    if near_distance_tight is None:
        near_distance_tight = 0.01 * radius

    n_near_each = int(n_samples * near_fraction / 2)   # per band
    n_surface   = int(n_samples * 0.10)
    n_volume    = n_samples - 2 * n_near_each - n_surface

    all_points = []
    all_sdfs   = []

    # --- Volumetric samples (25% of total) ---
    pts_vol = rng.uniform(bbox_min, bbox_max, size=(n_volume, 3))
    sdf_vol = np.linalg.norm(pts_vol - center, axis=1) - radius
    all_points.append(pts_vol)
    all_sdfs.append(sdf_vol)

    # --- Tight near-surface (32.5%): |sdf| < 0.01R ---
    pts_t, sdf_t = _rejection_sample_near_surface(
        rng, n_near_each, bbox_min, bbox_max, center, radius, near_distance_tight)
    all_points.append(pts_t)
    all_sdfs.append(sdf_t)

    # --- Loose near-surface (32.5%): |sdf| < 0.1R ---
    pts_l, sdf_l = _rejection_sample_near_surface(
        rng, n_near_each, bbox_min, bbox_max, center, radius, near_distance)
    all_points.append(pts_l)
    all_sdfs.append(sdf_l)

    # --- On-surface samples: random point on sphere, perturb along normal ---
    # Random unit vectors give uniform points on the sphere surface
    dirs = rng.standard_normal(size=(n_surface, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    surface_pts = center + radius * dirs
    # Perturb along normal (= unit direction from centre) by Gaussian noise
    perturb = rng.normal(scale=near_distance * 0.3, size=n_surface)
    pts_surf = surface_pts + dirs * perturb[:, np.newaxis]
    sdf_surf = np.linalg.norm(pts_surf - center, axis=1) - radius
    all_points.append(pts_surf)
    all_sdfs.append(sdf_surf)

    points = np.vstack(all_points)
    sdfs   = np.concatenate(all_sdfs)

    # Shuffle
    idx = rng.permutation(len(points))
    return points[idx], sdfs[idx]


# ---------------------------------------------------------------------------
# PyTorch Dataset
# ---------------------------------------------------------------------------

class SdfDataset(Dataset):
    """
    PyTorch Dataset for SDF training.

    Normalises input coordinates per-axis to [-1, 1] using the bounding box.
    SDF values are kept in world units (not normalised).

    Args:
        points  : ndarray (N, 3) world-space coordinates
        sdfs    : ndarray (N,)   signed distance values
        bbox_min: ndarray (3,)   bounding box minimum
        bbox_max: ndarray (3,)   bounding box maximum
    """

    def __init__(self, points, sdfs, bbox_min, bbox_max):
        self.bbox_min = np.asarray(bbox_min, dtype=np.float64)
        self.bbox_max = np.asarray(bbox_max, dtype=np.float64)

        # Normalise coordinates per axis to [-1, 1]
        points_norm = 2.0 * (points - self.bbox_min) / (self.bbox_max - self.bbox_min) - 1.0

        self.coords = torch.tensor(points_norm, dtype=torch.float64)
        self.sdfs   = torch.tensor(sdfs, dtype=torch.float64).unsqueeze(1)

    def __len__(self):
        return len(self.coords)

    def __getitem__(self, idx):
        return self.coords[idx], self.sdfs[idx]


def make_train_val_split(dataset, val_fraction=0.1, seed=42):
    """
    Split a dataset into training and validation subsets.

    Returns:
        train_dataset, val_dataset
    """
    n_val   = max(1, int(len(dataset) * val_fraction))
    n_train = len(dataset) - n_val
    generator = torch.Generator().manual_seed(seed)
    return random_split(dataset, [n_train, n_val], generator=generator)
