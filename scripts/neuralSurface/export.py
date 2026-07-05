"""
Export trained NeuralSDF weights to the binary format read by SCONE's mlpWeightIO_mod.f90.

Binary file layout (all little-endian, matching native x86 byte order):
  [Header]
    magic_number     int32   = 0x4E534446 ("NSDF")
    version          int32   = 1
    input_dim        int32   (always 3)
    hidden_dim       int32
    num_layers       int32
    activation_type  int32   (1=LeakyReLU, 2=ReLU, 3=tanh)
    leaky_alpha      float64
    sdf_scale        float64

  [Normalisation]
    bbox_min         float64 × 3
    bbox_max         float64 × 3

  [Validation test vector]
    test_input       float64 × 3   (one test point in world coordinates)
    test_output      float64       (expected model output at test_input)

  [Layers 1 .. num_layers]
    weight           float64 matrix in Fortran column-major order
    bias             float64 vector

CRITICAL: Weight matrices are written in Fortran column-major (column-by-column) order.
PyTorch stores weights as (out_features, in_features) in C row-major order.
We write W.numpy().flatten(order='F') which reorders to column-major so that
Fortran's matmul(W_fort(out,in), h(in)) gives the correct result.
The embedded test vector catches any transpose/byte-order errors at load time.
"""

import struct
import numpy as np
import torch

MAGIC_NUMBER    = 0x4E534446
FORMAT_VERSION  = 1


def export_weights(model, filename, bbox_min, bbox_max, test_point=None):
    """
    Write model weights to a binary file in the NSDF format.

    Args:
        model      : Trained NeuralSDF instance
        filename   : Output file path (e.g. 'sphere_weights.bin')
        bbox_min   : Array-like of shape (3,) — training bounding box minimum
        bbox_max   : Array-like of shape (3,) — training bounding box maximum
        test_point : Optional array-like of shape (3,) in world coordinates.
                     If None, uses the centre of the bounding box.
                     Used to embed a validation test vector in the file.
    """
    bbox_min = np.asarray(bbox_min, dtype=np.float64)
    bbox_max = np.asarray(bbox_max, dtype=np.float64)

    # Choose test point in world coordinates
    if test_point is None:
        test_point = 0.5 * (bbox_min + bbox_max)
    test_point = np.asarray(test_point, dtype=np.float64)

    in_dim = getattr(model, 'in_dim', 3)

    # Pad bbox and test_point to 3 elements (Fortran reader always reads 3)
    bbox_min_3 = np.zeros(3); bbox_min_3[:in_dim] = bbox_min[:in_dim]
    bbox_max_3 = np.ones(3);  bbox_max_3[:in_dim] = bbox_max[:in_dim]
    test_point_3 = np.zeros(3); test_point_3[:in_dim] = test_point[:in_dim]

    # Evaluate model at test point using only in_dim coordinates
    test_point_norm = _normalise(test_point[:in_dim], bbox_min[:in_dim], bbox_max[:in_dim])
    model.eval()
    with torch.no_grad():
        x = torch.tensor(test_point_norm, dtype=torch.float64).unsqueeze(0)
        model_f64 = _to_float64(model)
        test_output = float(model_f64(x).item())

    with open(filename, 'wb') as f:
        # Header
        f.write(struct.pack('<i', MAGIC_NUMBER))
        f.write(struct.pack('<i', FORMAT_VERSION))
        f.write(struct.pack('<i', in_dim))                     # input_dim
        f.write(struct.pack('<i', model.hidden_dim))
        f.write(struct.pack('<i', model.num_layers))
        f.write(struct.pack('<i', model.activation_type_code))
        f.write(struct.pack('<d', model.leaky_alpha))
        f.write(struct.pack('<d', model.sdf_scale))

        # Normalisation bounding box (always 3 × float64 for Fortran reader)
        f.write(bbox_min_3.tobytes())
        f.write(bbox_max_3.tobytes())

        # Embedded test vector (always 3 × float64 for Fortran reader)
        f.write(test_point_3.tobytes())
        f.write(struct.pack('<d', test_output))

        # Layer weights and biases
        for i, layer in enumerate(model.linear_layers):
            W = layer.weight.detach().cpu().to(torch.float64).numpy()  # shape (out, in)
            b = layer.bias.detach().cpu().to(torch.float64).numpy()    # shape (out,)

            # Write W in Fortran column-major order so Fortran's
            # matmul(W_fort(out,in), h(in)) gives the correct result.
            # W.flatten(order='F') = column-by-column = Fortran native.
            f.write(W.flatten(order='F').astype(np.float64).tobytes())
            f.write(b.astype(np.float64).tobytes())

    print(f"Exported weights to: {filename}")
    print(f"  Architecture: {model.num_layers} layers, {model.hidden_dim} hidden, "
          f"activation={model.activation}, sdf_scale={model.sdf_scale:.4f}")
    print(f"  Parameters:   {model.parameter_count():,}")
    print(f"  Test point:   {test_point}")
    print(f"  Test output:  {test_output:.8f}")
    print(f"  Bbox:         {bbox_min} to {bbox_max}")


def export_text(model, filename, bbox_min, bbox_max, test_point=None):
    """
    Write model weights to a plain-text file (one value per line).
    Same logical order as the binary format. Useful for debugging.

    The magic number is written as its decimal integer value so the
    Fortran text reader can verify format consistency.
    """
    bbox_min = np.asarray(bbox_min, dtype=np.float64)
    bbox_max = np.asarray(bbox_max, dtype=np.float64)

    if test_point is None:
        test_point = 0.5 * (bbox_min + bbox_max)
    test_point = np.asarray(test_point, dtype=np.float64)

    test_point_norm = _normalise(test_point, bbox_min, bbox_max)
    model.eval()
    with torch.no_grad():
        x = torch.tensor(test_point_norm, dtype=torch.float64).unsqueeze(0)
        model_f64 = _to_float64(model)
        test_output = float(model_f64(x).item())

    with open(filename, 'w') as f:
        f.write(f"# NSDF text weight file — generated by export.py\n")
        f.write(f"{MAGIC_NUMBER}\n")
        f.write(f"{FORMAT_VERSION}\n")
        f.write(f"3\n")                          # input_dim
        f.write(f"{model.hidden_dim}\n")
        f.write(f"{model.num_layers}\n")
        f.write(f"{model.activation_type_code}\n")
        f.write(f"{model.leaky_alpha:.17g}\n")
        f.write(f"{model.sdf_scale:.17g}\n")
        # bbox
        for v in bbox_min:
            f.write(f"{v:.17g}\n")
        for v in bbox_max:
            f.write(f"{v:.17g}\n")
        # test vector
        for v in test_point:
            f.write(f"{v:.17g}\n")
        f.write(f"{test_output:.17g}\n")
        # Layer weights and biases (Fortran column-major order)
        for layer in model.linear_layers:
            W = layer.weight.detach().cpu().to(torch.float64).numpy()
            b = layer.bias.detach().cpu().to(torch.float64).numpy()
            for val in W.flatten(order='F'):
                f.write(f"{val:.17g}\n")
            for val in b:
                f.write(f"{val:.17g}\n")

    print(f"Exported text weights to: {filename}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _normalise(point, bbox_min, bbox_max):
    """Normalise world-space point to [-1, 1] per axis."""
    return 2.0 * (point - bbox_min) / (bbox_max - bbox_min) - 1.0


def _to_float64(model):
    """Return a float64 copy of the model for precise test-vector evaluation."""
    import copy
    m = copy.deepcopy(model)
    m = m.double()
    return m
