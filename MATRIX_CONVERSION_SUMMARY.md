# Matrix Repository Conversion Summary

## Overview
Successfully converted MATLAB matrix repository functions to Python, implementing mass, damping, and stiffness matrix construction for rotorcraft dynamics analysis.

## Converted Functions

### 1. Mass Matrix (`build_mass_matrix`)
**Source:** `mass_matrix_repository.m`

**Implemented Configurations:**
- **3 blades:** 5×5 static mass matrix
- **4 blades:** 6×6 static mass matrix
- **5 blades:** 7×7 static mass matrix
- **7 blades:** 9×9 **time-varying** mass matrix with periodic coupling terms

**Key Features:**
- Time-varying coupling coefficients for 7-blade configuration
- Periodic terms at 2ω, 4ω, and 6ω frequencies
- Clear variable naming: `blade_inertia`, `blade_static_moment`, `total_hub_mass`

### 2. Damping & Stiffness Matrices (`build_damping_matrix`)
**Source:** `damping_matrix_repository.m`

**Implemented Configurations:**

#### 3-Blade Rotor (5×5 matrices)
- **H2B/ALL:** Static matrices with gyroscopic coupling
- **H2B/ODI:** Time-varying matrices with 1ω and 2ω harmonics

#### 4-Blade Rotor (6×6 matrices)
- **H2B/ALL:** Static matrices with gyroscopic coupling
- **H2B/ODI:** Time-varying matrices with 1ω and 2ω harmonics

#### 5-Blade Rotor (7×7 matrices)
- **ALL:** Static matrices with dual gyroscopic modes (1ω and 2ω)
- **ODI:** Complex time-varying matrices with 1ω, 2ω, and 4ω harmonics

**Not Implemented:**
- **7-blade configuration:** Skipped due to extreme complexity (740+ lines with symbolic expressions)
- **Additional 4-blade modes:** "2DI ADJ", "2DI OPP", "3DI"
- **B2B configuration:** Blade-to-Blade damper connection

## Code Structure

```
src/utils/matrix_repositories.py
├── build_mass_matrix()              # Main entry point for mass matrices
├── build_damping_matrix()           # Main entry point for C and K matrices
├── _build_matrices_3_blade()        # 3-blade specific implementation
├── _build_matrices_4_blade()        # 4-blade specific implementation
└── _build_matrices_5_blade()        # 5-blade specific implementation
```

## Variable Naming Improvements

| MATLAB | Python | Description |
|--------|--------|-------------|
| `rotor_OMEGA` | `rotor_omega` | Rotor angular velocity |
| `blade_inertia_Ib` | `blade_inertia` | Blade moment of inertia |
| `blade_static_Sb` | `blade_static_moment` | Blade static moment |
| `nominal_damping_Cd` | `nominal_damping` | Nominal damping coefficient |
| `nominal_stiffness_Kd` | `nominal_stiffness` | Nominal stiffness coefficient |
| `lag_hinge_offset_e` | `lag_hinge_offset` | Lag hinge offset |
| N/A | `centrifugal_stiffness` | Computed centrifugal stiffening term |
| N/A | `gyroscopic_term` | Computed gyroscopic coupling term |
| N/A | `dynamic_stiffness` | Computed frequency-dependent stiffness |

## Test Results

### Test Coverage
✅ **3-blade:** H2B/ALL and H2B/ODI modes
✅ **4-blade:** H2B/ALL and H2B/ODI modes
✅ **5-blade:** ALL and ODI modes
✅ **7-blade:** Time-varying mass matrix
✅ **Rotor speed variation:** Testing from 10-40 rad/s
✅ **Matrix properties:** Symmetry, positive definiteness

### Physical Property Verification

**Mass Matrix (4-blade example):**
- Shape: 6×6
- Symmetric: ✅ True
- Positive definite: ✅ True
- Eigenvalues: All positive

**Damping Matrix (4-blade example):**
- Shape: 6×6
- Diagonal elements: All positive ✅
- Gyroscopic coupling: Proportional to rotor speed ✅

**Stiffness Matrix (4-blade example):**
- Shape: 6×6
- Speed-dependent: ✅ Verified
- Centrifugal stiffening: Increases with ω² ✅

### Sample Matrix Values

#### 4-Blade, H2B/ALL at 30 rad/s (286.5 RPM)

**Mass Matrix M (symmetric):**
```
M[0,0] = 4340.0 kg·m²       (collective inertia)
M[1,1] = 2170.0 kg·m²       (cyclic inertia)
M[4,4] = 8.41e+03 kg        (hub mass x-direction)
```

**Damping Matrix C:**
```
C[0,0] = 16,268 Nms/rad     (collective damping)
C[1,1] = 8,134 Nms/rad      (cyclic damping)
C[1,2] = 130,200 Nms/rad    (gyroscopic coupling)
```

**Stiffness Matrix K:**
```
K[0,0] = 317,224 Nm/rad     (collective stiffness + centrifugal)
K[1,1] = -1,794,388 Nm/rad  (cyclic stiffness - inertial softening)
K[1,2] = 244,020 Nm/rad     (gyroscopic coupling)
```

## Key Implementation Details

### Time-Varying Matrices
For ODI (On-Damper Individually) activation, matrices vary periodically with time:

```python
wt = rotor_omega * time
cos_wt = np.cos(wt)
sin_wt = np.sin(wt)
cos_2wt = np.cos(2 * wt)
sin_2wt = np.sin(2 * wt)
```

### Gyroscopic Coupling
Gyroscopic terms appear in both damping and stiffness matrices:

```python
# Damping matrix - off-diagonal coupling
C[1,2] = blade_inertia * rotor_omega
C[2,1] = -blade_inertia * rotor_omega

# Stiffness matrix - proportional to damping
K[1,2] = nominal_damping * rotor_omega
K[2,1] = -nominal_damping * rotor_omega
```

### Centrifugal Stiffening
Stiffness increases with rotor speed due to centrifugal forces:

```python
centrifugal_stiffness = blade_static_moment * lag_hinge_offset * rotor_omega**2
K[0,0] = nominal_stiffness + centrifugal_stiffness
```

### Dynamic Stiffness
Blade modes show inertial softening at higher speeds:

```python
dynamic_stiffness = (
    nominal_stiffness
    - blade_inertia * rotor_omega**2  # inertial softening
    + centrifugal_stiffness            # centrifugal stiffening
)
```

## Usage Example

```python
from src.config.problem_definition import ProblemDefinition
from src.utils.matrix_repositories import build_mass_matrix, build_damping_matrix

# Create problem configuration
problem = ProblemDefinition(
    number_blades=4,
    damper_connection="H2B",
    damper_activation="ODI"
)

# Build mass matrix
time = 0.0  # seconds
rotor_omega = 25.0  # rad/s (about 240 RPM)
M = build_mass_matrix(problem, time, rotor_omega)

# Build damping and stiffness function handles
C_func, K_func = build_damping_matrix(problem)

# Evaluate at specific time and speed
C = C_func(time, rotor_omega)
K = K_func(time, rotor_omega)

print(f"Mass matrix shape: {M.shape}")
print(f"Damping matrix shape: {C.shape}")
print(f"Stiffness matrix shape: {K.shape}")
```

## Future Work

### High Priority
1. Implement 7-blade damping/stiffness matrices (very complex)
2. Add B2B (Blade-to-Blade) damper configuration
3. Implement remaining 4-blade modes: "2DI ADJ", "2DI OPP", "3DI"

### Medium Priority
1. Add matrix validation functions
2. Implement symbolic differentiation for time-varying terms
3. Optimize matrix construction for real-time applications

### Low Priority
1. Add visualization functions for matrix patterns
2. Implement matrix export to MATLAB format
3. Add unit tests with known analytical solutions

## Notes

- The negative stiffness values in K[1,1] are physically correct - they represent the inertial softening effect at higher rotor speeds
- Time-varying matrices require careful numerical integration
- The 7-blade case has extremely complex coupling terms that may require symbolic computation
- All matrices are returned as NumPy arrays with float64 precision

## Files Created/Modified

1. `src/utils/matrix_repositories.py` - Main implementation (468 lines)
2. `test_matrix_repositories.py` - Comprehensive test suite (231 lines)
3. `MATRIX_CONVERSION_SUMMARY.md` - This documentation

## Testing

Run comprehensive tests:
```bash
python test_matrix_repositories.py
```

Expected output: All tests pass with detailed matrix information for each configuration.
