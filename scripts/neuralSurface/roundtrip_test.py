"""
Numpy-only round-trip test: Python binary writer → Fortran reader.

Creates a 2-layer MLP (3→2→1) with known weights, writes the NSDF
binary file using the same format as export.py, then prints the
expected output so it can be verified against Fortran.

Run:
    python3 roundtrip_test.py

The script writes /tmp/roundtrip_test.bin.
Run the Fortran companion test with:
    ./Build/mlpRoundtrip_sandbox.out
"""

import struct
import math
import numpy as np

MAGIC_NUMBER   = 0x4E534446
FORMAT_VERSION = 1

# ---------------------------------------------------------------------------
# Network definition (matches mlpInference_sandbox.f90 Test setup)
# Layer 1: W = [[1,0,0],[0,1,0]], b = [0,0]
# Layer 2: W = [[1,1]],           b = [0]
# Activation: LeakyReLU(0.01), sdfScale=1.0, bbox=[-1,1]^3
# ---------------------------------------------------------------------------

W1 = np.array([[1.0, 0.0, 0.0],
               [0.0, 1.0, 0.0]], dtype=np.float64)   # shape (2,3)
b1 = np.array([0.0, 0.0],       dtype=np.float64)

W2 = np.array([[1.0, 1.0]],     dtype=np.float64)   # shape (1,2)
b2 = np.array([0.0],            dtype=np.float64)

hidden_dim      = 2
num_layers      = 2
activation_type = 1          # 1 = LeakyReLU
leaky_alpha     = 0.01
sdf_scale       = 1.0
bbox_min        = np.array([-1.0, -1.0, -1.0], dtype=np.float64)
bbox_max        = np.array([ 1.0,  1.0,  1.0], dtype=np.float64)

# ---------------------------------------------------------------------------
# Forward pass in numpy (reference implementation)
# ---------------------------------------------------------------------------

def leaky_relu(x, alpha=0.01):
    return np.where(x >= 0, x, alpha * x)

def normalise(point, bmin, bmax):
    return 2.0 * (point - bmin) / (bmax - bmin) - 1.0

def forward(point):
    x = normalise(point, bbox_min, bbox_max)
    h = leaky_relu(W1 @ x + b1)      # (2,)
    z = W2 @ h + b2                   # (1,)
    return math.tanh(z[0]) * sdf_scale

# Test point: [0.5, 0.5, 0.5]  — same as sandbox Test 2
test_point  = np.array([0.5, 0.5, 0.5], dtype=np.float64)
test_output = forward(test_point)
print(f"Test point : {test_point}")
print(f"Test output: {test_output:.17g}   (expected tanh(1.0) = {math.tanh(1.0):.17g})")

# ---------------------------------------------------------------------------
# Write binary file (same format as export.py)
# ---------------------------------------------------------------------------

out_path = '/tmp/roundtrip_test.bin'

with open(out_path, 'wb') as f:
    # Header
    f.write(struct.pack('<i', MAGIC_NUMBER))
    f.write(struct.pack('<i', FORMAT_VERSION))
    f.write(struct.pack('<i', 3))              # input_dim
    f.write(struct.pack('<i', hidden_dim))
    f.write(struct.pack('<i', num_layers))
    f.write(struct.pack('<i', activation_type))
    f.write(struct.pack('<d', leaky_alpha))
    f.write(struct.pack('<d', sdf_scale))

    # Normalisation bbox
    f.write(bbox_min.tobytes())
    f.write(bbox_max.tobytes())

    # Embedded test vector
    f.write(test_point.tobytes())
    f.write(struct.pack('<d', test_output))

    # Layer 1: W1 in Fortran column-major order, then b1
    # W1 shape (2,3): flatten(order='F') goes column-by-column
    f.write(W1.flatten(order='F').astype(np.float64).tobytes())
    f.write(b1.astype(np.float64).tobytes())

    # Layer 2: W2 in Fortran column-major order, then b2
    f.write(W2.flatten(order='F').astype(np.float64).tobytes())
    f.write(b2.astype(np.float64).tobytes())

print(f"\nWrote: {out_path}")
print(f"  Header:      magic=0x{MAGIC_NUMBER:08X}, version={FORMAT_VERSION}")
print(f"  Arch:        input=3, hidden={hidden_dim}, layers={num_layers}, "
      f"act={activation_type}, alpha={leaky_alpha}, scale={sdf_scale}")
print(f"  Bbox:        {bbox_min} to {bbox_max}")
print(f"\nNow run: ./Build/mlpRoundtrip_sandbox.out")
