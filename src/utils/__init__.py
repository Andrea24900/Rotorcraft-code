"""Utility modules for plotting and matrix operations."""

from .plotting import PlotProperties, MyPlot
from .matrix_repositories import build_mass_matrix, build_damping_matrix

__all__ = [
    'PlotProperties',
    'MyPlot',
    'build_mass_matrix',
    'build_damping_matrix',
]
