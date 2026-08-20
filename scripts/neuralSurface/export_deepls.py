"""
Export a paper-faithful DeepLS model (shared decoder + per-voxel latent
codes) to the binary format read by SCONE's deepLSWeightIO_mod.f90.

Binary file layout (all little-endian):
  [Header]
    magic_number     int32   = 0x4E534446 ("NSDF")
    version          int32   = 3
    latent_dim       int32
    hidden_dim       int32   (shared decoder architecture)
    num_layers       int32
    activation_type  int32
    leaky_alpha      float64
    nvox_x/y/z       int32 x3
    grid_origin      float64 x3
    voxel_size       float64 x3
    n_active_voxels  int32

  [Voxel Map]  (nvox_x*nvox_y*nvox_z int32, Fortran column-major, ix fastest)
    status: 0 = constant inside, 1 = constant outside, 2 = has latent code

  [Shared Decoder]  (written once — same layer serialisation as export.py's
                     single-MLP format: Fortran column-major weight matrices)
    [Layers 1..num_layers] weight (out,in) col-major, bias (out)

  [Active Voxel Latent Codes]  (n_active_voxels blocks, same traversal order
                                as status==2 entries in the Voxel Map)
    voxel_index(3)  : int32 x3   (1-based, self-validation)
    bbox_min(3)     : float64 x3 (local normalisation bbox — extended 1.5x
                                   receptive field, same convention as v2)
    bbox_max(3)     : float64 x3
    latent(latent_dim) : float64 x latent_dim

  [Validation test vector]  (first active voxel, decoder+latent combined)
    test_input(3)   : float64 x3  (world-space point inside the first active voxel)
    test_output     : float64     (expected combined evaluate() output there)

This is dramatically more compact than the v2 (independent-MLP) format: one
shared decoder (same size as a global MLP) instead of N independent decoders,
plus a small latent vector per voxel instead of a full weight set per voxel.
"""

import struct
import numpy as np
import torch

MAGIC_NUMBER   = 0x4E534446
FORMAT_VERSION = 3

STATUS_CONST_INSIDE  = 0
STATUS_CONST_OUTSIDE = 1
STATUS_HAS_MLP       = 2


def _to_float64(model):
    import copy
    m = copy.deepcopy(model)
    return m.double()


def export_deepls(filename, status_grid, grid_origin, voxel_size,
                  decoder, latents, active_voxels, voxel_bboxes,
                  latent_dim, hidden_dim, num_layers, activation_type, leaky_alpha):
    """
    Args:
        decoder       : trained NeuralSDF (in_dim = latent_dim+3), CPU, any dtype
        latents       : ndarray (n_active, latent_dim) — row i corresponds to
                        active_voxels[i]
        active_voxels : list of (ix,iy,iz) 0-based, in the SAME order as `latents`'
                        rows and as encountered scanning status_grid in Fortran
                        column-major order (ix fastest) — caller must ensure this,
                        train_deepls.py's classification loop already does.
        voxel_bboxes  : dict {(ix,iy,iz) -> (bbox_min, bbox_max)}
    """
    nvox_x, nvox_y, nvox_z = status_grid.shape
    grid_origin = np.asarray(grid_origin, dtype=np.float64)
    voxel_size  = np.asarray(voxel_size, dtype=np.float64)
    n_active = len(active_voxels)

    # Sanity: active_voxels must appear in Fortran column-major scan order,
    # matching the Voxel Map's traversal — verify against status_grid directly
    # rather than trusting the caller silently.
    expected_order = [(ix, iy, iz)
                       for iz in range(nvox_z)
                       for iy in range(nvox_y)
                       for ix in range(nvox_x)
                       if status_grid[ix, iy, iz] == STATUS_HAS_MLP]
    if expected_order != list(active_voxels):
        raise ValueError("active_voxels order does not match status_grid's Fortran "
                         "column-major scan order — Fortran reader would desync.")

    decoder64 = _to_float64(decoder)

    with open(filename, 'wb') as f:
        f.write(struct.pack('<i', MAGIC_NUMBER))
        f.write(struct.pack('<i', FORMAT_VERSION))
        f.write(struct.pack('<i', latent_dim))
        f.write(struct.pack('<i', hidden_dim))
        f.write(struct.pack('<i', num_layers))
        f.write(struct.pack('<i', activation_type))
        f.write(struct.pack('<d', leaky_alpha))
        f.write(struct.pack('<i', nvox_x))
        f.write(struct.pack('<i', nvox_y))
        f.write(struct.pack('<i', nvox_z))
        f.write(grid_origin.tobytes())
        f.write(voxel_size.tobytes())
        f.write(struct.pack('<i', n_active))

        f.write(status_grid.astype(np.int32).tobytes(order='F'))

        for layer in decoder64.linear_layers:
            W = layer.weight.detach().cpu().numpy()
            b = layer.bias.detach().cpu().numpy()
            f.write(W.flatten(order='F').astype(np.float64).tobytes())
            f.write(b.astype(np.float64).tobytes())

        first_test_input = None
        first_test_output = None
        for i, (ix, iy, iz) in enumerate(active_voxels):
            bbox_min, bbox_max = voxel_bboxes[(ix, iy, iz)]
            bbox_min = np.asarray(bbox_min, dtype=np.float64)
            bbox_max = np.asarray(bbox_max, dtype=np.float64)
            latent = latents[i].astype(np.float64)

            f.write(struct.pack('<iii', ix + 1, iy + 1, iz + 1))
            f.write(bbox_min.tobytes())
            f.write(bbox_max.tobytes())
            f.write(latent.tobytes())

            if first_test_input is None:
                first_test_input = 0.5 * (bbox_min + bbox_max)
                xyz_norm = 2.0 * (first_test_input - bbox_min) / (bbox_max - bbox_min) - 1.0
                decoder_in = np.concatenate([latent, xyz_norm])
                with torch.no_grad():
                    x = torch.tensor(decoder_in, dtype=torch.float64).unsqueeze(0)
                    first_test_output = float(decoder64(x).item())

        if n_active > 0:
            f.write(np.asarray(first_test_input, dtype=np.float64).tobytes())
            f.write(struct.pack('<d', first_test_output))

    print(f"Exported DeepLS (shared decoder) to: {filename}")
    print(f"  Voxel grid: {nvox_x}x{nvox_y}x{nvox_z} ({n_active} active)")
    print(f"  Shared decoder: {num_layers} layers, {hidden_dim} hidden")
    print(f"  Latent dim: {latent_dim}")
