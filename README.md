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
│   ├── stability_analysis/  # LTI, LTP, HD analysis
│   ├── analysis/            # Continuation methods
│   └── utils/               # Plotting utilities
├── tests/                   # Unit tests (78 total)
├── pyproject.toml           # Package configuration
└── main_*.py                # Example scripts
```

## Installation

### From Release (Recommended)
```bash
# Install specific version
pip install git+https://github.com/Andrea24900/rotor_dynamics_python.git@v1.0.1

# Or clone specific version
git clone --branch v1.0.1 https://github.com/Andrea24900/rotor_dynamics_python.git
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

Edit `src/config/problem_definition.py`:

```python
@dataclass
class ProblemDefinition:
    number_blades: int = 4
    damper_connection: Literal["H2B", "B2B"] = "H2B"
    damper_activation: Literal["ALL", "ODI", "CUSTOM"] = "ALL"
    damper_activation_vector: np.ndarray  # For CUSTOM mode, e.g., np.array([0,1,1,1])

    lower_rotor_RPM: float = 20.0
    higher_rotor_RPM: float = 400.0
    number_points: int = 100

    required_solver: Literal["LTI", "LTP", "HD"] = "LTI"
    number_harmonics: int = 1          # HD analysis
    hd_use_complex: bool = False       # Complex HD formulation
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
from src.config.problem_definition import ProblemDefinition
from src.rotor_build import RotorBuild
from src.stability_analysis.continuation_analysis import ContinuationAnalysis

# B2B requires using ProblemDefinition constructor for proper geometric initialization
problem = ProblemDefinition(
    number_blades=4,
    damper_connection="B2B",
    damper_activation="ODI",
    required_solver="LTP",
    continuation="YES"
)

rotor = RotorBuild(problem)
analysis = ContinuationAnalysis(rotor).continuation()
```

**Note**: For B2B configuration, always use the `ProblemDefinition` constructor with `damper_connection="B2B"` to ensure proper initialization of B2B geometric parameters (`l0`, `lE`, `gamma_E`, `phi`).

## Features

### Implemented ✅
- LTI, LTP, HD stability analysis
- Continuation methods (OMEGA and parametric)
- Mass, damping, stiffness matrices for 3/4/5/7 blades
- H2B and B2B damper configurations
- Complex HD formulation
- Visualization tools
- 78 comprehensive tests

### In Progress 🚧
- Modal participation analysis

## Conversion Status

| Module | Status | Notes |
|--------|--------|-------|
| Core Infrastructure | ✅ | Configuration, matrices |
| LTI/LTP/HD Analysis | ✅ | All methods complete |
| Continuation | ✅ | OMEGA and parametric |
| Plotting | ✅ | Damping/frequency plots |
| Modal Participation | 🚧 | Structure exists |

## Recent Improvements (January 2026)

### v1.0.2 - B2B Damper Configuration
- **B2B Matrix Generation**: Added `damping_activation_B2B()` function in `matrix_generation_GR.py` implementing blade-to-blade damper matrices with:
  - Diagonal terms: contributions from dampers on both sides of each blade (`C_before`, `C_after`)
  - Extra-diagonal terms: coupling between adjacent blades (`C_zetaed`)
  - Stiffness coefficients: `K_before`, `K_after`, `K_zetaed` with proper geometric formulation
- **B2B Test Suite**: Added continuation tests for B2B configuration (`tests_to_run_2.py`):
  - LTI continuation with ALL dampers
  - LTP continuation with ODI
  - HD continuation with ODI (1 harmonic)
- **Separate Output Directories**: B2B figures saved to `extra_figures/B2B/`, H2B to `extra_figures/H2B/`
- **ProblemDefinition Constructor**: All test files updated to use constructor approach for proper B2B geometric initialization

### v1.0.1 - API Simplification & Code Cleanup
- **Simplified API Signatures**:
  - `ltp_single_point(OMEGA)` and `hd_single_point(OMEGA)` now only require OMEGA; period T computed internally as `T = 2π/OMEGA`
  - `monodromy_computer(A_handle, OMEGA, rtol, atol)` signature with OMEGA instead of T
  - `hd_computer(A_handle, OMEGA, number_harmonics, time_samples, use_complex)` signature
- **Unified Matrix Generation**: New `build_MCK_matrices()` function in `matrix_generation_GR.py` replaces separate `build_mass_matrix` and `build_damping_matrix` functions
- **Configurable Integration Tolerances**: Added `ode_rtol` and `ode_atol` parameters to `ProblemDefinition` for ODE solver control
- **Docstring Cleanup**: Reduced verbose docstrings to concise documentation while preserving essential information
- **Removed Legacy Code**: Deleted MATLAB `OOP_build/` directory and obsolete `matrix_repositories.py`
- **HDStability Inheritance**: Now inherits directly from `StabilityAnalysis` instead of `LTPStability`

### v1.0.0 - System Matrices Refactor
- **Matrix Generation**: Refactored `matrix_generation_GR.py` with proper MBC transformation
- **LTP Performance**: Vectorized monodromy integration (single n² ODE instead of n separate integrations)
- **Damper Activation**: Fixed CUSTOM mode with proper `damper_activation_vector` support
- **Mass Matrix Reuse**: `build_damping_matrix` accepts pre-computed `mass_matrix_rot` to avoid redundant computation
- **Integration Tolerances**: Adjusted to (rtol=1e-4, atol=1e-6) for better performance

### Key Enhancements
- **Continuation**: Fixed conjugation bug, parametric continuation support, code deduplication
- **HD Analysis**: Complex formulation option, vectorized computation (2-5x faster)
- **Nomenclature**: Unified to "HD" (Harmonic Decomposition), snake_case methods
- **Validation**: Comprehensive input checks with clear error messages
- **Performance**: Configurable tolerances, automatic warnings for large systems
- **Testing**: 78 tests with 100% pass rate

## Results

Comprehensive stability analysis results and visualizations for all methods (LTI, LTP, HD) are available in [RESULTS.md](RESULTS.md), including:
- Damping and frequency plots across RPM ranges
- Continuation analysis demonstrations
- Method comparisons and validation
- Complete bibliography

## Verification

All code verified against MATLAB implementation:
- Physical parameters and matrix formulations
- HD computer with FFT coefficients
- 8 H-matrix assignment functions

See [CONVERSION_GUIDE.md](CONVERSION_GUIDE.md) for details.

## Dependencies

- NumPy (>=1.24.0), SciPy (>=1.10.0), Matplotlib (>=3.7.0)
- Dev: pytest, pytest-cov

Automatically installed with `pip install -e ".[dev]"`

## Development

### Running Tests (78 total)

```bash
# All tests
python -m pytest tests/ -v

# With coverage
python -m pytest tests/ --cov=src --cov-report=term-missing

# Specific test file
python -m pytest tests/test_lti_stability.py -v
```

**Test breakdown:**
- LTI stability: 27 tests
- LTP stability: 32 tests
- HD stability: 19 tests
- Core modules: 63 tests (matrices, rotor build, config)

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
  version = {1.0.1},
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
