# Modal Participation Analysis

## ⚠️ Important Notice

**WARNING: Evaluation of modal participation computational logic still in progress.**

This implementation is under active development and validation. Results should be carefully reviewed and cross-checked against known references.

## Overview

Modal participation analysis shows how each harmonic component contributes to each mode, providing insight into the frequency content of the system's dynamic behavior.

## Quick Start

**HD Method:**
```bash
python main_hd_modal_participation.py
```

**LTP Method:**
```bash
python main_ltp_modal_participation.py
```

⚠️ **LTP Warning**: LTP modal participation is dependent on the integer multiple of Omega added to the complex logarithm. Caution is advised when interpreting results.

Both scripts:
1. Build rotor system and run continuation analysis
2. Compute modal participation factors
3. Display overview damping plot
4. Prompt for state and mode pair selection
5. Generate detailed plots

## Output Dimensions

**φ_jkn** - Participation factor for state j, mode k, harmonic n

- **j (states)**: 1-12 → [1-4: MBC velocity | 5-6: Hub velocity | 7-10: MBC | 11-12: Hub]
- **k (mode pairs)**: Complex conjugate pairs (mode pair k → eigenvalues 2k-1 and 2k)
- **n (harmonics)**: Frequency components (negative, zero, positive)

**Properties**: Range [0,1], Σ_n φ_jkn = 1.0 for each (j,k) pair

## Plots

1. **Overview Damping**: All eigenvalues vs RPM (guides mode selection)
2. **Individual Harmonics**: Each n separately (solid: n<0, dash-dot: n=0, dashed: n>0)
3. **Summed Harmonics**: Combined ±n contributions
4. **Mode Damping**: Selected mode pair stability (negative = stable)

## Example Usage

Run script → View system info → Review overview damping plot → Select state (1-12) and mode pair → View four detailed plots

## Programmatic Usage

```python
from src.rotor_build import RotorBuild
from src.stability_analysis.continuation_analysis import ContinuationAnalysis
from src.analysis.modal_participation import ModalParticipationAnalysis
from src.utils.plotting import MyPlot

rotor = RotorBuild.build_all()
rotor.problem.required_solver = "HD"
rotor.problem.continuation = "YES"

cont_analysis = ContinuationAnalysis(rotor).continuation()
modal_part = ModalParticipationAnalysis(cont_analysis)

# Access results: modal_part.modal_participation[i].phi_jkn (shape: states × modes × harmonics)
# Plot: MyPlot().plot_mod_part(..., state_index=6, mode_index=4, ...)  # 0-indexed
```

## Configuration

Adjust harmonics (line 78 in main_hd_modal_participation.py):
```python
rotor.problem.number_harmonics = 4  # Typical: 1-8. Higher = more detail, more computation
```

Adjust RPM range:
```python
rotor.problem.lower_rotor_RPM = 10
rotor.problem.higher_rotor_RPM = 500
rotor.problem.number_points = 50
```

## Files

- `src/analysis/modal_participation.py` - Core implementation (LTP: Floquet+FFT, HD: harmonic extraction)
- `src/utils/plotting.py` - Plotting functions
- `main_hd_modal_participation.py` - Interactive HD application
- `main_ltp_modal_participation.py` - Interactive LTP application

## Implementation Notes

**HD Method**: Eigenvectors already expanded in harmonics [X_0, X_1c, X_1s, ...]. Converts to complex exponential form: C^(±n) = (X_nc ± i·X_ns)/2

**LTP Method**:
⚠️ Warning: Dependent on integer multiple of Omega in complex logarithm for Floquet exponents.
1. Integrate transition matrix φ(t) over one period
2. Compute V(t) = φ(t)·eigenvectors·exp(-η·t)
3. Apply FFT to extract 2N harmonics on each side

## Troubleshooting

- **"No module named 'src'"**: Run from project root or add to PYTHONPATH
- **Plots not displaying**: Configure matplotlib backend (TkAgg or Qt5Agg)
- **Out of memory**: Reduce `number_harmonics` or `number_points`

## References

MATLAB implementation: `OOP_build/modal_participation_analysis.m`

---

**Last Updated**: 2026-01-13
