# Rotor Dynamics Analysis - Python

Python implementation of rotorcraft dynamics stability analysis, converted from MATLAB.

## Overview

This project provides tools for analyzing the stability of rotorcraft systems using various methods:
- **LTI (Linear Time-Invariant)** stability analysis
- **LTP (Linear Time-Periodic)** stability analysis using Floquet theory
- **HD (Harmonic Decomposition)** stability analysis using the method from Lopez and Prasad

## Project Structure

```
rotor_dynamics_python/
├── src/                          # Source code
│   ├── config/                   # Configuration and problem definition
│   ├── stability_analysis/       # Stability analysis modules
│   ├── analysis/                 # Modal participation and continuation
│   └── utils/                    # Utilities and plotting
├── scripts/                      # Executable scripts
├── tests/                        # Unit tests
│   ├── test_problem_definition.py
│   └── test_package_config.py    # Tests for pyproject.toml validation
├── pyproject.toml                # Modern package configuration (PEP 621)
├── setup.py                      # Minimal setup file for compatibility
└── requirements.txt              # Legacy dependencies (use pyproject.toml)
```

## Installation

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd rotor_dynamics_python
   ```

2. **Create a virtual environment** (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install the package**:

   **For development (editable install - recommended):**
   ```bash
   pip install -e .
   ```
   This allows you to make changes to the code and see them reflected immediately.

   **For normal use:**
   ```bash
   pip install .
   ```

   **With development dependencies (includes pytest, pytest-cov):**
   ```bash
   pip install -e ".[dev]"
   ```

The package uses modern Python packaging standards with `pyproject.toml` (PEP 621 compliant).

## Quick Start

### Basic Stability Analysis

Run the LTI stability analysis test:
```bash
python main_lti_test.py
```

Or use the API directly:
```python
from src.rotor_build import RotorBuild
from src.stability_analysis.base import StabilityAnalysis

# Build rotor system
rotor = RotorBuild.build_all()

# Run stability analysis (automatically selects solver based on config)
modes = StabilityAnalysis.run_stability(rotor)

# Access results
for sol in modes.modal_solution:
    print(f"RPM: {sol.OMEGA_RPM:.1f}, Dampings: {sol.damping}")
```

### Running the Test Script

The `main_lti_test.py` script performs a complete LTI stability analysis:
- Builds a 4-blade rotor system with H2B dampers
- Computes eigenvalues across the RPM range (20-400 RPM)
- Generates damping and frequency plots
- Shows all 12 eigenvalues (dampings) and positive frequencies only

## Configuration

The default configuration can be modified in `src/config/problem_definition.py`:

```python
@dataclass
class ProblemDefinition:
    # Rotor configuration
    number_blades: Literal[3, 4, 5, 7] = 4
    damper_connection: Literal["H2B", "B2B"] = "H2B"
    damper_activation: Literal["ALL", "ODI", "2DI ADJ", "2DI OPP", "3DI"] = "ODI"

    # Solution parameters
    lower_rotor_RPM: float = 20.0
    higher_rotor_RPM: float = 400.0
    number_points: int = 50

    # Solver configuration
    required_solver: Literal["LTI", "LTP", "HD"] = "LTI"
    number_harmonics: int = 1  # for HD analysis
```

### Configuration Options

- **number_blades**: Blade count (3, 4, 5, or 7)
- **damper_connection**:
  - `H2B`: Hub-to-Blade dampers
  - `B2B`: Blade-to-Blade dampers
- **damper_activation**:
  - `ALL`: All dampers active (LTI system)
  - `ODI`: One Damper Inactive per mode (LTP system)
- **required_solver**: Analysis method (LTI, LTP, or HD)
- **RPM range**: Operating speed range for analysis
- **number_harmonics**: Harmonics to include in HD analysis

## Features

### Implemented ✅
- **Problem Definition**: Complete dataclass-based configuration system
- **Mass Matrix Repository**: All blade configurations (3, 4, 5, 7) including time-varying matrices
- **Damping & Stiffness Matrices**: H2B and B2B configurations with ALL/ODI damper activation
- **LTI Stability Analysis**: Eigenvalue analysis across RPM range
- **LTP Stability Analysis**: Monodromy matrix computation via Floquet theory
- **HD Stability Analysis**: Complete Harmonic Decomposition method with FFT-based harmonic decomposition
- **Continuation Analysis**: Predictor-corrector method for eigenvalue tracking across parameter ranges
- **State Matrix Construction**: Full state-space formulation for all configurations
- **Visualization**: Damping and frequency plots vs RPM (generic scatter plots and mode-tracked line plots)

### In Progress 🚧
- Modal participation analysis

### Planned ⏳
- Documentation and examples

## Conversion Status

| Module | Status | Files | Notes |
|--------|--------|-------|-------|
| **Core Infrastructure** | ✅ Complete | `src/config/`, `src/utils/matrix_repositories.py` | All data structures and matrix formulations |
| **Rotor Build** | ✅ Complete | `src/rotor_build.py` | State-space construction |
| **Stability Base** | ✅ Complete | `src/stability_analysis/base.py` | HD_computer with 8 H-matrix helpers |
| **LTI Stability** | ✅ Complete | `src/stability_analysis/lti_stability.py` | Eigenvalue analysis |
| **LTP Stability** | ✅ Complete | `src/stability_analysis/ltp_stability.py` | Monodromy matrices |
| **Continuation** | ✅ Complete | `src/stability_analysis/continuation_analysis.py` | Predictor-corrector tracking |
| **Plotting** | ✅ Complete | `src/utils/plotting.py` | Generic and mode-tracked plots |
| **Modal Participation** | 🚧 Partial | `src/analysis/modal_participation.py` | Structure exists |

Legend: ✅ Complete | 🚧 In Progress | ⏳ Not Started

## Verification Against MATLAB

All converted code has been verified against the original MATLAB implementation:

### Data Verification ✅
- **Rotor Characteristics**: All physical parameters match exactly
  - Blade properties (mass, inertia, static moment)
  - Hub properties (mass, damping, stiffness)
  - Damper coefficients (H2B and B2B configurations)
  - Geometric parameters (B2B configuration)

### Matrix Formula Verification ✅
- **Mass Matrices**: All coefficients verified for 3, 4, 5, and 7-blade configurations
- **Damping Matrices**: All time-varying terms match MATLAB expressions
- **Stiffness Matrices**: Gyroscopic and centrifugal terms verified
- **HD Computer**: All 8 H-matrix assignment functions match MATLAB logic

### Implementation Details
See `CONVERSION_GUIDE.md` for detailed documentation of:
- Matrix repository implementation
- Harmonic Decomposition computer with FFT coefficients
- H-matrix block assembly algorithms
- Harmonic interaction formulas

## Dependencies

### Runtime Dependencies
- **NumPy** (>=1.24.0): Array operations and linear algebra
- **SciPy** (>=1.10.0): ODE solvers, interpolation, FFT
- **Matplotlib** (>=3.7.0): Plotting and visualization

### Development Dependencies
- **pytest** (>=7.0.0): Testing framework
- **pytest-cov** (>=4.0.0): Test coverage reporting
- **tomli**: TOML file parsing (Python <3.11)

All dependencies are automatically installed when you run `pip install -e .` or `pip install -e ".[dev]"`

## Development

### Running Tests

**Run all tests:**
```bash
python -m pytest tests/ -v
```

**Run specific test file:**
```bash
python -m pytest tests/test_package_config.py -v
```

**Run with coverage report:**
```bash
python -m pytest tests/ --cov=src --cov-report=term-missing
```

**Run tests without pytest:**
```bash
python tests/test_package_config.py
```

### Package Configuration Tests

The `test_package_config.py` file validates your `pyproject.toml` configuration:
- ✅ TOML file exists and is valid
- ✅ Build system properly configured
- ✅ Project metadata (name, version, description)
- ✅ Dependencies correctly specified
- ✅ Optional dev dependencies defined
- ✅ Authors and license information
- ✅ Package discovery configuration
- ✅ Python version requirements

### Building the Package

**Build distribution packages (wheel and source):**
```bash
python -m build
```

This creates installable packages in the `dist/` directory.

### Package Structure

The package follows modern Python packaging standards:
- `pyproject.toml`: Modern PEP 621 compliant configuration
- `setup.py`: Minimal backward-compatibility shim
- All package metadata, dependencies, and build configuration in `pyproject.toml`

### Code Style
This project follows PEP 8 style guidelines.

## License

GNU

## Contact

Andrea Bassi
GitHub: Andrea24900

## Acknowledgments

From my (Andrea Bassi) thesis work.
