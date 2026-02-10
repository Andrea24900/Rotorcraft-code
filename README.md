# Rotor Dynamics Analysis - Python

[![Version](https://img.shields.io/github/v/tag/Andrea24900/Rotorcraft-code)](https://github.com/Andrea24900/Rotorcraft-code/tags)

Python implementation of rotorcraft dynamics stability analysis.

## Overview

Stability analysis tools for rotorcraft systems:
- **LTI (Linear Time-Invariant)**: Eigenvalue analysis
- **LTP (Linear Time-Periodic)**: Floquet theory with monodromy matrix
- **HD (Harmonic Decomposition)**: Frequency domain method (Lopez-Prasad)

## Project Structure

```
rotor_dynamics_python/
├── src/                     # Source code
│   ├── config/              # Configuration
    |── analysis/            # Modal participation and parameteric analysis
│   ├── stability_analysis/  # LTI, LTP, HD, sorting-free Hill analysis
│   ├── analysis/            # Continuation methods
│   └── utils/               # Plotting utilities
├── tests/                   # No tests (code already validated in known cases)
├── pyproject.toml           # Package configuration
└── main_*.py                # Example scripts
```

## Installation

### From Release (Recommended)
```bash
# Install specific version
pip install git+https://github.com/Andrea24900/rotor_dynamics_python.git@v1.3.0

# Or clone specific version
git clone --branch v1.3.0 https://github.com/Andrea24900/rotor_dynamics_python.git
cd rotor_dynamics_python
pip install -e ".[dev]"
```

### From Source (Development)
```bash
git clone https://github.com/Andrea24900/rotor_dynamics_python.git
cd rotor_dynamics_python
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Development install (recommended)
pip install -e ".[dev]"

# Or without dev dependencies
pip install -e .
```

## Quick Start

```bash
# Run example analysis
python main_lti_test.py
```

Or use the API:
```python
from src.rotor_build import RotorBuild
from src.stability_analysis.base import StabilityAnalysis

rotor = RotorBuild.build_all()
modes = StabilityAnalysis.run_stability(rotor)

for sol in modes.modal_solution:
    print(f"RPM: {sol.OMEGA_RPM:.1f}, Dampings: {sol.damping}")
```

## Continuation Analysis

Robust eigenvalue tracking using predictor-corrector algorithm.

**Example:**
```bash
python main_lti_continuation_test.py
```

```python
from src.stability_analysis.continuation_analysis import ContinuationAnalysis

rotor = RotorBuild.build_all()
rotor.problem.continuation = "YES"

cont = ContinuationAnalysis(rotor).continuation()

for sol in cont.modal_solution:
    print(f"RPM: {sol.OMEGA_RPM:.1f}, Damping: {sol.damping}")
```

### Parametric Continuation

Continuation with respect to any parameter (not just RPM):

```python
# Sensitivity w.r.t. custom parameter
sensitivity = ContinuationAnalysis.numerical_sensitivity(
    matrix_type='SS',        # 'SS', 'MON', or 'HD'
    A_handle_all=state_matrix_func,
    OMEGA=100.0,
    wrt_omega=False,         # False = use 'par' parameter
    par=5.0                  # Parameter value
)
```

**Configuration** (`problem_definition.py`):
```python
continuation: Literal["YES", "NO"] = "NO"
continuation_tolerance: float = 0.01
continuation_max_iter: int = 100
```

## Configuration

Edit `src/config/problem_definition.py`. The configuration is split into two classes:

- **`RotorCharacteristics`**: physical properties of the rotor (blades, dampers, geometry)
- **`ProblemDefinition`**: analysis settings (solver, RPM range, tolerances)

```python
@dataclass
class RotorCharacteristics:
    number_blades: int = 4
    damper_connection: Literal["H2B", "B2B"] = "H2B"
    damper_activation: Literal["ALL", "ODI", "CUSTOM"] = "ALL"
    damper_activation_vector: np.ndarray  # For CUSTOM mode, e.g., np.array([0,1,1,1])
    # ... blade, hub, and damper physical properties

@dataclass
class ProblemDefinition:
    lower_rotor_RPM: float = 20.0
    higher_rotor_RPM: float = 400.0
    number_points: int = 100

    required_solver: Literal["LTI", "LTP", "HD"] = "LTI"
    number_harmonics: int = 1          # HD analysis
    hd_use_complex: bool = False       # Complex HD formulation

    rotor_characteristics: RotorCharacteristics = ...
```

### Damper Connection Types
- **H2B (Hub-to-Blade)**: Each damper connects hub to individual blade (diagonal damping matrix)
- **B2B (Blade-to-Blade)**: Each damper connects adjacent blades cyclically (damping matrix with diagonal and extra-diagonal coupling terms)

### Damper Activation Modes
- **ALL**: All dampers active (LTI system)
- **ODI**: One Damper Inoperative - first damper disabled (LTP system)
- **CUSTOM**: User-defined via `damper_activation_vector` (e.g., `np.array([0,1,1,1])` disables first damper)

### B2B Configuration Example

```python
from src.config.problem_definition import ProblemDefinition, RotorCharacteristics
from src.rotor_build import RotorBuild
from src.stability_analysis.continuation_analysis import ContinuationAnalysis

problem = ProblemDefinition(
    rotor_characteristics=RotorCharacteristics(
        number_blades=4, damper_connection="B2B", damper_activation="ODI"),
    required_solver="LTP",
    continuation="YES"
)

rotor = RotorBuild(problem)
analysis = ContinuationAnalysis(rotor).continuation()
```

**Note**: B2B geometric parameters (`l0`, `lE`, `gamma_E`, `phi`) are computed automatically by `RotorCharacteristics.__post_init__` when `damper_connection="B2B"`.

## Features

### Implemented ✅
- LTI, LTP, HD stability analysis
- Continuation methods (OMEGA and parametric)
- Mass, damping, stiffness matrices for arbitrary number of blades
- H2B and B2B damper configurations
- Complex HD formulation
- Parametric analysis (OMEGA and one parameter)
- Visualization tools
- Sorting-free Hill's method

### In Progress 🚧
- Modal participation analysis

## Recent Improvements (February 2026)
### v1.3.0 - Sorting-free stability method
- **Sorting-free stability**: Added the 'SortingFreeStability class, which implements the method by Bayer (2023) to compute the Floquet modes from Hill's matrix, aka the complex HD matrix. 
- **Improvement of plots**: Added some features to the plot_***_generic method. 

### v1.2.0 - Code generalisation
- **Separation of variables**: Separated the solver options from the rotor characteristics in 'problem_definition'. Adapted code accordingly. Updated and clarified some comments. 
- **Cleanup of old code**: Deleted old code with outdated tests.

## Results

Comprehensive stability analysis results and visualizations for all methods (LTI, LTP, HD) are available in [RESULTS.md](RESULTS.md), including:
- Damping and frequency plots across RPM ranges
- Continuation analysis demonstrations
- Method comparisons and validation
- Complete bibliography

## Dependencies

- NumPy (>=1.24.0), SciPy (>=1.10.0), Matplotlib (>=3.7.0)
- Dev: pytest, pytest-cov

Automatically installed with `pip install -e ".[dev]"`

### Building

```bash
python -m build  # Creates dist/ packages
```

## License

GNU

## Citation

If you use this software in your research, please cite:

```bibtex
@software{bassi2026rotor,
  author = {Bassi, Andrea},
  title = {Rotor Dynamics Analysis - Python},
  year = {2026},
  version = {1.1.0},
  url = {https://github.com/Andrea24900/rotor_dynamics_python}
}

@mastersthesis{bassi2024thesis,
  author = {Bassi, Andrea},
  title = {Periodic ground resonance analysis with non-linear elements},
  school = {Politecnico di Milano},
  year = {2024}
}
```

## Contact

Andrea Bassi
GitHub: [@Andrea24900](https://github.com/Andrea24900)
Politecnico di Milano

## Acknowledgments

Based on the Master's thesis *Periodic ground resonance analysis with non-linear elements* (2024) and the original MATLAB implementation developed for that work.
