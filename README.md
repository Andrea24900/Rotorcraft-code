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
│   ├── test_problem_definition.py  # Configuration tests (13 tests)
│   ├── test_rotor_build.py         # Rotor build tests (16 tests)
│   ├── test_matrix_repositories.py # Matrix construction tests (25 tests)
│   └── test_package_config.py      # Package config validation (9 tests)
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

### Continuation Analysis

The continuation method provides **robust eigenvalue tracking** across parameter ranges by using a predictor-corrector algorithm. This avoids mode crossing issues common in standard eigenvalue analysis.

#### Running Continuation Analysis

**Test script:**
```bash
python main_lti_continuation_test.py
```

**API usage (continuation w.r.t. OMEGA):**
```python
from src.rotor_build import RotorBuild
from src.stability_analysis.continuation_analysis import ContinuationAnalysis

# Build rotor system
rotor = RotorBuild.build_all()
rotor.problem.continuation = "YES"

# Run continuation analysis
cont_analysis = ContinuationAnalysis(rotor)
cont_analysis = cont_analysis.continuation()

# Access tracked eigenvalue branches
for sol in cont_analysis.modal_solution:
    print(f"RPM: {sol.OMEGA_RPM:.1f}")
    print(f"Dampings: {sol.damping}")
    print(f"Frequencies: {sol.frequency}")
```

#### Parametric Continuation (New Feature!)

The continuation method now supports **parametric continuation** with respect to any system parameter, not just rotor speed (OMEGA).

**Example - Continuation w.r.t. damping coefficient:**
```python
from src.stability_analysis.continuation_analysis import ContinuationAnalysis

# Define state matrix function with parameter
def state_matrix_with_damping(t, damping_coeff, omega):
    """State matrix A(t, c, Ω) where c is the damping parameter."""
    # Your state matrix construction with variable damping_coeff
    return A

# Create continuation analysis
cont_analysis = ContinuationAnalysis(rotor)

# Compute sensitivity w.r.t. damping coefficient (not OMEGA)
sensitivity_matrix = ContinuationAnalysis.numerical_sensitivity(
    matrix_type='SS',           # or 'MON', 'HD'
    A_handle_all=state_matrix_with_damping,
    step_size_h=0.01,
    OMEGA=100.0,                # Fixed rotor speed
    wrt_omega=False,            # Key parameter!
    par=5.0                     # Current damping coefficient value
)
```

**Key Parameters:**
- `wrt_omega=True` (default): Compute derivative w.r.t. rotor speed OMEGA
- `wrt_omega=False`: Compute derivative w.r.t. parameter `par`
  - When using `wrt_omega=False`, you **must** provide `par` value
  - Your `A_handle_all` function signature changes:
    - For OMEGA: `A(omega)` or `A(t, omega)`
    - For parameter: `A(par, omega)` or `A(t, par, omega)`

**Supported Matrix Types:**
- `'SS'`: State-space (LTI systems)
- `'MON'`: Monodromy matrix (LTP systems)
- `'HD'`: Harmonic Decomposition (time-periodic systems)

**Use Cases for Parametric Continuation:**
- Damping coefficient variation studies
- Stiffness parameter sensitivity
- Geometric parameter sweeps
- Control parameter optimization
- Bifurcation analysis with respect to design parameters

#### Continuation Method Details

The continuation algorithm uses:
- **Predictor step**: First-order sensitivity-based prediction of next eigenvalue
- **Corrector step**: Newton-Raphson iterations to refine the eigenvalue
- **Normalization constraint**: Ensures unique eigenvector determination

**Configuration parameters** (in `problem_definition.py`):
```python
continuation: Literal["YES", "NO"] = "NO"
continuation_tolerance: float = 0.01    # Convergence tolerance
continuation_max_iter: int = 100        # Max corrector iterations
step_h: float = 0.1                     # Finite difference step size
```

**Benefits of continuation method:**
- ✅ Smooth eigenvalue branch tracking (no mode jumping)
- ✅ Works with degenerate/clustered eigenvalues
- ✅ Handles bifurcations naturally
- ✅ More robust than standard eigenvalue sorting
- ✅ Supports parametric studies beyond just RPM variation

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
    hd_use_complex: bool = False  # Use complex HD formulation
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
- **hd_use_complex**: Use complex formulation for HD solver (faster for many harmonics)

## Features

### Implemented ✅
- **Problem Definition**: Complete dataclass-based configuration system
- **Mass Matrix Repository**: All blade configurations (3, 4, 5, 7) including time-varying matrices
- **Damping & Stiffness Matrices**: H2B and B2B configurations with ALL/ODI damper activation
- **LTI Stability Analysis**: Eigenvalue analysis across RPM range
- **LTP Stability Analysis**: Monodromy matrix computation via Floquet theory with configurable tolerances
- **HD Stability Analysis**: Complete Harmonic Decomposition method with:
  - FFT-based harmonic decomposition
  - Real and complex formulation options
  - Vectorized performance optimizations
  - Comprehensive input validation
- **Continuation Analysis**: Predictor-corrector method for eigenvalue tracking across parameter ranges
- **State Matrix Construction**: Full state-space formulation for all configurations
- **Visualization**: Damping and frequency plots vs RPM (generic scatter plots and mode-tracked line plots)
- **Performance Features**:
  - Automatic warnings for large system computations
  - Configurable integration tolerances for monodromy matrix computation
  - Optimized array operations for HD method

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

## Recent Improvements (January 2026)

### Continuation Analysis Refactoring (January 2026)

The continuation analysis module has been significantly refactored for improved clarity, performance, and flexibility:

#### 1. Fixed Critical Mathematical Bug
- **Conjugation inconsistency**: Corrected normalization constraint in sensitivity systems
- All three augmented systems now use mathematically correct `2 * v^H` (conjugate transpose)
- Improved predictor accuracy and convergence speed
- The corrector was already correct, which is why code worked but with suboptimal prediction

#### 2. Nomenclature Consistency
- **HD vs HB**: Unified all references to use "HD" (Harmonic Decomposition) consistently
- Updated variable names: `A_HB*` → `A_HD*`
- Updated parameter names: `HB_parameters` → `HD_parameters`
- Clearer and more consistent codebase

#### 3. Removed Dead Code
- Eliminated unused forward difference (FW) derivative code
- Kept only centered difference (CEN) implementation
- Reduced `numerical_sensitivity` from ~105 to ~78 lines
- Code is now simpler and more maintainable

#### 4. Parametric Continuation Support
- **New feature**: Continuation with respect to any parameter, not just OMEGA
- Replaced confusing `par=1000` magic number with explicit `wrt_omega` boolean
- Clear API: `wrt_omega=True` (OMEGA derivative) or `wrt_omega=False` (parameter derivative)
- Enables parametric studies: damping coefficients, stiffness, geometric parameters, etc.

#### 5. Code Deduplication and Refactoring
- **Extracted helper methods**: Eliminated code duplication in continuation algorithm
  - `_build_augmented_matrix()`: Constructs augmented system matrix (used in 3 places)
  - `_compute_residual()`: Computes residual for corrector (used in 2 places)
  - `_solve_augmented_system()`: Solves with conditioning check (used in 3 places)
- **Reduced code size**: Main continuation method is now ~60% cleaner
- **Improved maintainability**: Changes to augmented system only need to be made once
- **Better documentation**: Each helper method has clear docstrings explaining purpose

#### 6. Bug Fixes
- Fixed `state_matrix_A_handles` → `state_matrix_function` attribute name
- Tests now run successfully

### Stability Analysis Enhancements

The stability analysis module has been significantly improved with better performance, robustness, and user experience:

#### 1. Comprehensive Input Validation
- All critical methods (`monodromy_computer`, `hd_computer`, `hd_computer_complex`) now validate inputs
- Clear error messages for invalid parameters (negative periods, non-square matrices, invalid harmonics)
- Protection against None returns and type mismatches

#### 2. Configurable Tolerances
- Integration tolerances for monodromy computation are now configurable
- Class constants: `DEFAULT_RTOL = 1e-6`, `DEFAULT_ATOL = 1e-8`
- Can be overridden per-call: `monodromy_computer(A, T, rtol=1e-8, atol=1e-10)`

#### 3. Performance Optimizations
- **Vectorized HD computation**: 2-5x faster for large time arrays
- Replaced loop-based array filling with vectorized operations
- Applied to both real and complex HD formulations

#### 4. Complex HD Formulation
- New configuration option: `hd_use_complex = True` in problem definition
- Uses complex exponentials instead of sine/cosine terms
- Automatically selected when `use_complex=True` in `hd_computer()`

#### 5. Performance Warnings
- Automatic warnings for large system computations (n > 50)
- Informs users about expected computation time
- Suggests using reduced-order models for n > 100

#### 6. Code Quality Improvements
- DRY principle applied to RPM-to-rad/s conversions
- Better code organization and maintainability
- 35 new comprehensive tests for stability analysis base class

#### 7. Method Naming Conventions
- **Breaking Change**: All methods now follow PEP 8 snake_case convention
- Renamed: `HD_computer` → `hd_computer`
- Renamed: `HD_computer_complex` → `hd_computer_complex`
- Renamed: `LTP_single_point` → `ltp_single_point`
- Renamed: `LTP_full_range` → `ltp_full_range`
- Renamed: `HD_single_point` → `hd_single_point`
- Renamed: `HD_full_range` → `hd_full_range`
- Improved code consistency across the entire codebase

#### 8. Test Coverage
- Total test suite: 98 tests (35 for stability analysis alone)
- All tests passing with 100% success rate
- Covers edge cases, validation, and integration scenarios

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

The project includes **98 comprehensive tests** covering all core functionality.

**Run all tests:**
```bash
python -m pytest tests/ -v
```

**Run specific test file:**
```bash
python -m pytest tests/test_matrix_repositories.py -v
python -m pytest tests/test_rotor_build.py -v
python -m pytest tests/test_problem_definition.py -v
python -m pytest tests/test_package_config.py -v
python -m pytest tests/test_stability_base.py -v
```

**Run with coverage report:**
```bash
python -m pytest tests/ --cov=src --cov-report=term-missing
```

**Run tests without pytest (standalone):**
```bash
python tests/test_matrix_repositories.py
python tests/test_rotor_build.py
python tests/test_problem_definition.py
python tests/test_package_config.py
python tests/test_stability_base.py
```

### Test Suite Overview

**test_stability_base.py** (35 tests) - Stability analysis core methods:
- ✅ Class constants and initialization
- ✅ Monodromy matrix computation (simple, periodic, custom tolerances)
- ✅ HD computer (real and complex formulations)
- ✅ Input validation for all methods
- ✅ Integration with rotor systems
- ✅ Performance warnings for large systems
- ✅ Edge cases (small periods, many harmonics)
- ✅ ModalSolution dataclass
- ✅ OMEGA range assignment (L2R and R2L)

**test_package_config.py** (9 tests) - Package configuration validation:
- ✅ TOML file exists and is valid
- ✅ Build system properly configured
- ✅ Project metadata (name, version, description)
- ✅ Dependencies correctly specified
- ✅ Optional dev dependencies defined
- ✅ Authors and license information
- ✅ Package discovery configuration
- ✅ Python version requirements

**test_problem_definition.py** (13 tests) - Configuration and rotor characteristics:
- ✅ Default values and custom initialization
- ✅ H2B and B2B damper configurations
- ✅ Geometric calculations for different blade counts
- ✅ Solver types (LTI, LTP, HD)
- ✅ Damper activation modes
- ✅ RPM range and continuation parameters

**test_rotor_build.py** (16 tests) - State-space matrix construction:
- ✅ Initialization with/without auto_build
- ✅ Mass matrix properties (symmetry, positive definiteness)
- ✅ Mass matrix inverse validation
- ✅ State matrix dimensions and structure
- ✅ Time-varying behavior
- ✅ Different blade counts and damper modes
- ✅ Error handling and caching

**test_matrix_repositories.py** (25 tests) - Matrix repository functions:
- ✅ Mass matrix construction for 3, 4, 5, 7-blade rotors
- ✅ Matrix symmetry and positive definiteness
- ✅ Damping/stiffness matrices for H2B and B2B configurations
- ✅ ALL vs ODI damper activation modes
- ✅ Speed and time dependencies
- ✅ Gyroscopic coupling effects
- ✅ 7-blade coupling coefficient validation
- ✅ Numerical stability (no NaN/Inf values)

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

From my (Andrea Bassi) thesis work and the MATLAB implementation
