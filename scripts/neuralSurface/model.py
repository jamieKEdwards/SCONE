"""
NeuralSDF MLP architecture.

Maps 3D coordinates to a signed distance value.
Negative output = inside surface (negative halfspace).
Positive output = outside surface (positive halfspace).

Architecture follows Chabra et al., "Deep Local Shapes" (ECCV 2020):
  - num_layers fully-connected layers (includes output layer)
  - hidden_dim neurons per hidden layer
  - LeakyReLU (or ReLU / tanh) between hidden layers
  - tanh output activation, scaled by sdf_scale

The model operates on per-axis normalised coordinates in [-1, 1].
Normalisation is applied externally (see sampler.py / train.py).
sdf_scale is set after inspecting the training data range.

Activation type codes (must match mlpInference_mod.f90 constants):
  1 = LeakyReLU
  2 = ReLU
  3 = tanh
"""

import torch
import torch.nn as nn

ACTIVATION_LEAKYRELU = 1
ACTIVATION_RELU      = 2
ACTIVATION_TANH      = 3

_ACTIVATION_NAMES = {
    'leakyrelu': ACTIVATION_LEAKYRELU,
    'relu':      ACTIVATION_RELU,
    'tanh':      ACTIVATION_TANH,
}


class NeuralSDF(nn.Module):
    """
    MLP signed distance function.

    Args:
        hidden_dim   : Width of each hidden layer (default 128)
        num_layers   : Total number of weight matrices, including output layer (default 4)
        activation   : Hidden layer activation: 'leakyrelu', 'relu', or 'tanh' (default 'leakyrelu')
        leaky_alpha  : Negative slope for LeakyReLU (default 0.01)

    Architecture for num_layers=4, hidden_dim=128:
        Linear(3, 128)   + LeakyReLU
        Linear(128, 128) + LeakyReLU
        Linear(128, 128) + LeakyReLU
        Linear(128, 1)   + tanh * sdf_scale
    """

    def __init__(self, hidden_dim=128, num_layers=4, activation='leakyrelu', leaky_alpha=0.01,
                 in_dim=3):
        super().__init__()

        if num_layers < 2:
            raise ValueError("num_layers must be >= 2 (at least one hidden + output layer)")
        if activation not in _ACTIVATION_NAMES:
            raise ValueError(f"activation must be one of {list(_ACTIVATION_NAMES)}")
        if in_dim not in (2, 3):
            raise ValueError("in_dim must be 2 or 3")

        self.in_dim        = in_dim
        self.hidden_dim    = hidden_dim
        self.num_layers    = num_layers
        self.activation    = activation
        self.leaky_alpha   = leaky_alpha
        self.sdf_scale     = 1.0          # set via set_sdf_scale() before export

        # Build layer list
        layers = []
        _first = in_dim
        for i in range(num_layers):
            out_dim = hidden_dim if i < num_layers - 1 else 1
            layers.append(nn.Linear(_first, out_dim))
            _first = hidden_dim

        self.linear_layers = nn.ModuleList(layers)

        # Hidden activation function
        if activation == 'leakyrelu':
            self.hidden_act = nn.LeakyReLU(negative_slope=leaky_alpha)
        elif activation == 'relu':
            self.hidden_act = nn.ReLU()
        else:  # tanh
            self.hidden_act = nn.Tanh()

    def forward(self, x):
        """
        Forward pass.

        Args:
            x : Tensor of shape (batch, 3), coordinates normalised to [-1, 1]

        Returns:
            Tensor of shape (batch, 1) with signed distance values
        """
        h = x
        for i, layer in enumerate(self.linear_layers):
            h = layer(h)
            if i < self.num_layers - 1:
                h = self.hidden_act(h)
        # Output activation: tanh scaled by sdf_scale
        return torch.tanh(h) * self.sdf_scale

    def set_sdf_scale(self, scale):
        """Set the SDF output scale (multiplier on tanh output)."""
        if scale <= 0:
            raise ValueError("sdf_scale must be positive")
        self.sdf_scale = float(scale)

    @property
    def activation_type_code(self):
        """Integer code for Fortran weight file header."""
        return _ACTIVATION_NAMES[self.activation]

    def parameter_count(self):
        """Total number of trainable parameters."""
        return sum(p.numel() for p in self.parameters())
