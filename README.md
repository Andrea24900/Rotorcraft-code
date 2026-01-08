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
└── requirements.txt              # Python dependencies
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

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

## Quick Start

### Basic Stability Analysis

```python
from src.rotor_build import RotorBuild
from src.stability_analysis.base import StabilityAnalysis

# Build rotor system
rotor = RotorBuild.build_all()

# Run stability analysis
modes = StabilityAnalysis.run_stability(rotor)

# Plot results
# (plotting functionality to be added)
```

## Configuration

Edit `src/config/problem_definition.py` to configure:
- Number of blades (3, 4, 5, 7)
- Damper configuration (H2B or B2B)
- Analysis type (LTI, LTP, HD)
- RPM range and number of points
- Number of harmonics (for HD analysis)

## Features

- ✅ Multiple stability analysis methods
- ✅ Modal participation analysis
- ✅ Support for various rotor configurations
- ✅ Continuation methods for bifurcation analysis
- 🚧 Visualization tools (in progress)
- 🚧 Unit tests (in progress)

## Conversion Status

| Module | Status | Notes |
|--------|--------|-------|
| Problem Definition | 🚧 In Progress | Core data structures |
| Rotor Build | 🚧 In Progress | Main rotor class |
| LTI Stability | ⏳ Pending | - |
| LTP Stability | ⏳ Pending | - |
| HD Stability | ⏳ Pending | - |
| Modal Participation | ⏳ Pending | - |
| Plotting | ⏳ Pending | - |

Legend: ✅ Complete | 🚧 In Progress | ⏳ Pending

## Dependencies

- **NumPy**: Array operations and linear algebra
- **SciPy**: ODE solvers, interpolation, FFT
- **Matplotlib**: Plotting and visualization

## Development

### Running Tests
```bash
pytest tests/
```

### Code Style
This project follows PEP 8 style guidelines.

## License

GNU

## Contact

Andrea Bassi
GitHub: Andrea24900

## Acknowledgments

From my (Andrea Bassi) thesis work.
