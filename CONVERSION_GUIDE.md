# MATLAB to Python Conversion Guide

MATLAB rotorcraft dynamics analysis converted to Python/NumPy.

## Status: Complete ✅

All core functionality converted and tested:
- LTI, LTP, HD stability analysis
- Continuation methods (predictor-corrector)
- Matrix repositories (3/4/5/7 blades)
- Plotting system (matplotlib)
- 6 test scripts

## Test Scripts

```bash
# Standard methods
python main_lti_test.py  # LTI (eigenvalue)
python main_ltp_test.py  # LTP (Floquet)
python main_hd_test.py   # HD (harmonic decomposition)

# Continuation methods (branch tracking)
python main_lti_continuation_test.py
python main_ltp_continuation_test.py
python main_hd_continuation_test.py

# Comparison
python main_ltp_comparison.py
```

**Use continuation for:** bifurcations, mode crossings, parameter studies with degeneracies

## Conversion Status

| Component | Status | Python Files |
|-----------|--------|--------------|
| Core Infrastructure | ✅ | problem_definition.py, matrix_repositories.py |
| Rotor Build | ✅ | rotor_build.py |
| Stability Analysis | ✅ | lti_stability.py, ltp_stability.py, hd_stability.py |
| Continuation | ✅ | continuation_analysis.py |
| Plotting | ✅ | plotting.py |
| Test Scripts | ✅ | main_*_test.py (6 scripts) |
| Modal Participation | 🚧 | modal_participation.py (placeholder) |

## Key Components

### Matrix Repositories
- Mass, damping, stiffness for 3/4/5/7 blades
- H2B and B2B damper configurations
- ALL activation (all dampers active) or ODI activation (One Damper Inoperative)
- Time-varying terms and gyroscopic coupling
- Verified against MATLAB

### Stability Analysis
- **LTI**: Eigenvalue analysis
- **LTP**: Floquet theory with monodromy matrix
- **HD**: Harmonic decomposition with FFT (real and complex formulations)
- All 8 H-matrix helper functions for HD

### Continuation Algorithm
- **Predictor**: Sensitivity-based prediction `λ_next = λ + (dλ/dΩ)·ΔΩ`
- **Corrector**: Newton-Raphson iterations
- Tracks eigenvalue branches across parameter space
- Handles near-degeneracies and avoided crossings

### Plotting
- Generic scatter plots (all eigenvalues)
- Order-specific line plots (tracked modes)
- MATLAB-compatible color scheme

## Next Steps

1. Modal participation analysis
2. Campbell diagram plotting
3. Advanced continuation features (bifurcation detection, branch switching)

## Technical Details

### Mass Matrix
- 3-blade (5×5), 4-blade (6×6), 5-blade (7×7), 7-blade (9×9)
- 7-blade has time-varying trigonometric coupling
- Hub-blade coupling via static moment terms

### Damping & Stiffness
- H2B (Hub-to-Blade) and B2B (Blade-to-Blade) configurations
- ALL activation (LTI) or ODI activation (One Damper Inoperative - LTP)
- Gyroscopic coupling: `±n·Ib·ω`
- Centrifugal stiffening: `Sb·e·ω²`

### HD Computer
- FFT-based coefficient extraction
- HD matrix size: `(1+2N)×n` where N=harmonics, n=state_dim
- 8 helper functions for H-matrix block assembly
- Handles harmonic sum and difference interactions

### Continuation Algorithm
- **Predictor**: Sensitivity-based prediction `λ_next = λ + (dλ/dΩ)·ΔΩ`
- **Corrector**: Newton-Raphson iterations
- Tracks eigenvalue branches across parameter space
- Handles near-degeneracies and avoided crossings

### Plotting
- Generic scatter plots (all eigenvalues)
- Order-specific line plots (tracked modes)
- MATLAB-compatible color scheme

## MATLAB → Python Conversion Patterns

- Classes: `classdef` → `class`, `properties` → `__init__`
- Functions: `@(x)` → `lambda x:` or `def`
- Linear algebra: `A\b` → `np.linalg.solve(A, b)`, `eig(A)` → `np.linalg.eig(A)`
- ODE: `ode45` → `scipy.integrate.solve_ivp`

## Common Issues

- **Indexing**: MATLAB 1-based → Python 0-based
- **Imports**: Use relative imports within `src/`, absolute from scripts
- **Vectors**: Be explicit about array shapes (use `.reshape(-1, 1)` for column vectors)
