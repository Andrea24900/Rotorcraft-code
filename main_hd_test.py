"""Main file to test HD (Harmonic Decomposition) stability analysis.

This script performs HD stability analysis using the Lopez-Prasad method,
which expands time-periodic systems into an infinite-dimensional time-invariant
system in the frequency domain.

Harmonic Decomposition Method:
------------------------------
The HD method transforms the time-periodic system:
  dx/dt = A(t)x  where A(t+T) = A(t)

into an infinite-dimensional time-invariant system by expanding:
  A(t) = A₀ + Σ[Aₖᶜcos(kωt) + Aₖˢsin(kωt)]

The state vector is expanded as:
  x(t) = x₀ + Σ[xₖᶜcos(kωt) + xₖˢsin(kωt)]

This results in the Hill matrix (harmonic balance matrix):
  H = [H₀₀  H₀₁ᶜ  H₀₁ˢ  H₀₂ᶜ  ...]
      [H₁ᶜ₀ H₁ᶜ₁ᶜ H₁ᶜ₁ˢ ...]
      [H₁ˢ₀ H₁ˢ₁ᶜ H₁ˢ₁ˢ ...]
      [...]

The eigenvalues of H directly give the stability of the periodic system.

Advantages over LTP:
- No numerical integration (monodromy matrix)
- Direct eigenvalue computation
- Can capture high-frequency effects with sufficient harmonics

Disadvantages:
- Large matrices for high harmonic count: (1+2N)×n dimension
- Truncation error from finite harmonics
- More complex implementation

Usage:
------
python main_hd_test.py

Output:
-------
- Console: Analysis progress and stability results
- Plots: Damping and frequency vs RPM

Requirements:
-------------
- HD solver enabled
- Sufficient harmonics (typically N=1-3)
- Time samples for FFT (typically 200)
"""

import numpy as np
import matplotlib.pyplot as plt
from src.rotor_build import RotorBuild
from src.stability_analysis.hd_stability import HDStability
from src.utils.plotting import MyPlot


def main():
    """Run HD stability analysis and plot results.

    This function performs a complete HD stability analysis workflow:
    1. Build rotor system with HD solver
    2. Initialize HD stability analysis object
    3. Compute eigenvalues over RPM range using Harmonic Decomposition
    4. Display results and stability assessment
    5. Generate and display damping and frequency plots

    The HD method expands the periodic system into frequency domain,
    resulting in larger but time-invariant eigenvalue problems.

    Outputs:
    --------
    - Console output with progress information and key results
    - Damping plot: Real parts of eigenvalues vs RPM
    - Frequency plot: Imaginary parts of eigenvalues vs RPM

    Notes:
    ------
    - Matrix size: (1+2N)×n where N is number of harmonics, n is state dimension
    - For 4-blade rotor (n=12) with N=1: Hill matrix is 36×36
    - Computation time scales with (1+2N)³
    """

    print("="*60)
    print("HD (Harmonic Decomposition) Stability Analysis")
    print("="*60)

    # Step 1: Build rotor system
    print("\n[1/4] Building rotor system...")
    rotor = RotorBuild.build_all()

    # Override to ensure HD solver is used
    rotor.problem.required_solver = "HD"
    rotor.problem.number_harmonics = 1  # Start with 1 harmonic

    print(f"  [OK] Rotor system built successfully")
    print(f"    - Number of blades: {rotor.problem.number_blades}")
    print(f"    - Solver type: {rotor.problem.required_solver}")
    print(f"    - Damper connection: {rotor.problem.damper_connection}")
    print(f"    - RPM range: {rotor.problem.lower_rotor_RPM} - {rotor.problem.higher_rotor_RPM}")
    print(f"    - Number of points: {rotor.problem.number_points}")
    print(f"    - Number of harmonics (N): {rotor.problem.number_harmonics}")

    # Calculate expected matrix size
    state_size = 2 * rotor.problem.number_blades  # For 4 blades: 12
    hill_matrix_size = (1 + 2 * rotor.problem.number_harmonics) * state_size
    print(f"    - State space dimension: {state_size}")
    print(f"    - Hill matrix dimension: {hill_matrix_size}×{hill_matrix_size}")

    # Step 2: Create HD analysis object
    print("\n[2/4] Initializing HD stability analysis...")
    hd = HDStability(rotor)
    print(f"  [OK] HD analysis object created")

    # Step 3: Run HD analysis over full RPM range
    print("\n[3/4] Running HD analysis...")
    print(f"  Computing Hill matrices using FFT and harmonic expansion...")
    hd = hd.HD_full_range()

    print(f"\n  [OK] HD analysis completed")
    print(f"    - Number of solution points: {len(hd.modal_solution)}")
    print(f"    - Number of eigenvalues per point: {len(hd.modal_solution[0].damping)}")

    # Print sample results
    print(f"\n  Sample results at first operating point:")
    print(f"    - RPM: {hd.modal_solution[0].OMEGA_RPM:.2f}")
    print(f"    - Period T: {hd.modal_solution[0].T:.6f} s")
    print(f"    - First 6 damping values: {hd.modal_solution[0].damping[:6]}")
    print(f"    - First 6 frequency values: {hd.modal_solution[0].frequency[:6]}")

    # Check stability
    max_damping_overall = max([np.max(sol.damping) for sol in hd.modal_solution])
    if max_damping_overall < 0:
        print(f"\n  [OK] System is STABLE (max damping: {max_damping_overall:.4f})")
    else:
        print(f"\n  [WARNING] System is UNSTABLE (max damping: {max_damping_overall:.4f})")

    # Step 4: Generate plots
    print("\n[4/4] Generating plots...")

    # Create plotting object
    plotter = MyPlot()

    # Determine axis limits based on data
    rpm_min = rotor.problem.lower_rotor_RPM
    rpm_max = rotor.problem.higher_rotor_RPM

    all_damping = np.concatenate([sol.damping for sol in hd.modal_solution])
    all_frequency = np.concatenate([sol.frequency for sol in hd.modal_solution])

    # Damping limits
    damping_min = np.min(all_damping) - 0.5
    damping_max = max(1.0, np.max(all_damping) + 0.5)

    # Frequency limits - only positive frequencies
    positive_frequencies = all_frequency[all_frequency > 0]
    if len(positive_frequencies) > 0:
        freq_min = 0
        freq_max = np.max(positive_frequencies) + 5
    else:
        freq_min = 0
        freq_max = 50

    print(f"\n  Total number of eigenvalues per point: {len(hd.modal_solution[0].damping)}")
    print(f"  Note: HD method produces (1+2N)×state_size eigenvalues")

    # Plot 1: Damping (generic - all points)
    print(f"\n  - Creating damping plot (all eigenvalues, dots)...")
    fig1 = plotter.plot_damping_generic(
        hd.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(damping_min, damping_max)
    )
    plt.figure(fig1.number)
    plt.title(f"HD Damping vs RPM (N={rotor.problem.number_harmonics} harmonics)")
    plt.tight_layout()

    # Plot 2: Frequency (generic - all points)
    print(f"  - Creating frequency plot (all eigenvalues, dots)...")
    fig2 = plotter.plot_frequency_generic(
        hd.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(freq_min, freq_max)
    )
    plt.figure(fig2.number)
    plt.title(f"HD Frequency vs RPM (N={rotor.problem.number_harmonics} harmonics)")
    plt.tight_layout()

    print("\n" + "="*60)
    print("HD Analysis Complete!")
    print("="*60)
    print("\nThe Harmonic Decomposition method expands the periodic system")
    print("into frequency domain, solving a larger but time-invariant")
    print("eigenvalue problem. This avoids numerical integration but")
    print("requires sufficient harmonics for accuracy.")

    # Show plots
    print("\nDisplaying plots...")
    plt.show()


if __name__ == "__main__":
    main()
