# MATLAB to Python Conversion Guide

## Project Overview

This document tracks the conversion of a MATLAB rotorcraft dynamics analysis codebase to Python/NumPy.

## Conversion Status

### Phase 1: Core Infrastructure ✅ (Structure Complete)

| File | Status | Notes |
|------|--------|-------|
| `problem_definition.m` | ✅ Complete | Converted to dataclasses |
| `rotor_definition_script.m` | ✅ Complete | Integrated into ProblemDefinition |
| `plot_properties.m` | ✅ Complete | Converted to dataclasses |
| `mass_matrix_repository.m` | 🚧 Placeholder | Needs implementation |
| `damping_matrix_repository.m` | 🚧 Placeholder | Needs implementation |

### Phase 2: Core Classes ✅ (Structure Complete)

| File | Status | Notes |
|------|--------|-------|
| `rotor_build.m` | ✅ Structure | Core logic needs matrix implementations |
| `stability_analysis.m` | ✅ Structure | Base class with static methods |

### Phase 3: Stability Subclasses ✅ (Structure Complete)

| File | Status | Notes |
|------|--------|-------|
| `LTI_stability.m` | ✅ Complete | Basic implementation done |
| `LTP_stability.m` | ✅ Complete | Monodromy computation done |
| `HD_stability.m` | 🚧 Partial | HD_computer needs H-matrix functions |

### Phase 4: Analysis Tools 🚧 (Partial)

| File | Status | Notes |
|------|--------|-------|
| `modal_participation_analysis.m` | 🚧 Placeholder | Complex FFT logic needed |
| `continuation_analysis.m` | ⏳ Not Started | - |

### Phase 5: Plotting 🚧 (Partial)

| File | Status | Notes |
|------|--------|-------|
| `my_plot.m` | 🚧 Placeholder | Needs matplotlib implementation |

### Phase 6: Scripts 🚧 (One Example)

| File | Status | Notes |
|------|--------|-------|
| `CALL.m` | ✅ Structure | Basic workflow implemented |
| `CALL_DMD.m` | ⏳ Not Started | DMD analysis |
| `CALL_DMD_LTP.m` | ⏳ Not Started | DMD for LTP |
| `CALL_DMD_NL.m` | ⏳ Not Started | DMD for nonlinear |

## Next Steps Priority

### High Priority (Required for Basic Functionality)

1. **Implement mass_matrix_repository.py**
   - Convert switch/case logic for blades 3, 4, 5, 7
   - Matrix assembly for each configuration
   
2. **Implement damping_matrix_repository.py**
   - Time-dependent damping matrices
   - Stiffness matrices
   - Damper activation logic
   
3. **Complete HD_computer H-matrix functions**
   - H0M, HiCMjC, HiCMjS, HiSMjC, HiSMjS assignment
   - Matrix block assembly

### Medium Priority (For Full Features)

4. **Modal participation analysis**
   - LTP: FFT of V(t), coefficient extraction
   - HD: Harmonic component extraction
   
5. **Plotting implementation**
   - plot_damping_generic, plot_damping_order
   - plot_mod_part
   - Campbell diagrams

6. **DMD scripts**
   - Basic DMD implementation
   - LTP and nonlinear variants

### Low Priority (Advanced Features)

7. **Continuation analysis**
   - Bifurcation tracking
   - Predictor-corrector methods

## Key MATLAB → Python Patterns

### Class Conversion
```matlab
% MATLAB
classdef MyClass
    properties
        prop1
    end
    methods (Static)
        function obj = build()
        end
    end
    methods
        function obj = MyClass(arg)
        end
    end
end
```

```python
# Python
class MyClass:
    def __init__(self, arg):
        self.prop1 = None
    
    @staticmethod
    def build():
        pass
```

### Function Handles
```matlab
% MATLAB
A = @(t, omega) compute_matrix(t, omega);
result = A(0, 25);
```

```python
# Python
def A(t, omega):
    return compute_matrix(t, omega)

result = A(0, 25)
```

### Matrix Operations
```matlab
% MATLAB
x = A\b;              % → np.linalg.solve(A, b)
[V, D] = eig(A);      % → eigenvalues, eigenvectors = np.linalg.eig(A)
M_inv = inv(M);       % → M_inv = np.linalg.inv(M)
```

### ODE Integration
```matlab
% MATLAB
[t, x] = ode45(@odefun, [0 T], x0);
```

```python
# Python
from scipy.integrate import solve_ivp
sol = solve_ivp(odefun, [0, T], x0, method='DOP853')
t, x = sol.t, sol.y
```

## Testing Strategy

Each module should have unit tests in `tests/`:

1. **test_problem_definition.py** - Test configuration creation
2. **test_matrix_repositories.py** - Test matrix dimensions and symmetry
3. **test_rotor_build.py** - Test state matrix construction
4. **test_lti_stability.py** - Test eigenvalue computation
5. **test_ltp_stability.py** - Test monodromy matrix
6. **test_hd_stability.py** - Test Hill matrix assembly

## Common Issues & Solutions

### Issue: Import errors
**Solution**: Use relative imports within package, absolute from scripts
```python
# Within src/
from .config.problem_definition import ProblemDefinition

# From scripts/
from src.config.problem_definition import ProblemDefinition
```

### Issue: MATLAB's 1-based indexing
**Solution**: Carefully convert loops and array indexing
```matlab
% MATLAB (1-based)
for i = 1:n
    A(i, :) = ...
end
```

```python
# Python (0-based)
for i in range(n):
    A[i, :] = ...
```

### Issue: Row vs column vectors
**Solution**: Be explicit about array shapes
```python
# Column vector
x = np.array([[1], [2], [3]])  # shape (3, 1)
# or
x = np.array([1, 2, 3]).reshape(-1, 1)
```

## Development Workflow

1. **Create feature branch**: `git checkout -b feature/matrix-repository`
2. **Implement module**: Convert MATLAB logic to Python
3. **Add unit tests**: Create corresponding test file
4. **Test locally**: `pytest tests/test_module.py`
5. **Commit changes**: `git commit -m "Implement mass matrix repository"`
6. **Push and review**: `git push origin feature/matrix-repository`

## Resources

- NumPy documentation: https://numpy.org/doc/
- SciPy documentation: https://docs.scipy.org/doc/
- Matplotlib gallery: https://matplotlib.org/stable/gallery/
- Python style guide (PEP 8): https://pep8.org/

## Contact & Support

For questions about:
- **MATLAB implementation**: [Original author]
- **Python conversion**: [Your contact]
- **Mathematical theory**: [Reference papers/books]
