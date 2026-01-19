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

This results in the Harmonic Decomposition (HD) matrix A_HD, see Lopez.

The eigenvalues of A_HD directly give the stability of the periodic system.

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
    - HD matrix size: (1+2N)×n where N is number of harmonics (user-defined), n is state dimension
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
    rotor.problem.damper_activation = "ODI"
    rotor.problem.number_harmonics = 1  # Start with 1 harmonic

    print(f"  [OK] Rotor system built successfully")
    print(f"    - Number of blades: {rotor.problem.number_blades}")
    print(f"    - Solver type: {rotor.problem.required_solver}")
    print(f"    - Damper connection: {rotor.problem.damper_connection}")
    print(f"    - RPM range: {rotor.problem.lower_rotor_RPM} - {rotor.problem.higher_rotor_RPM}")
    print(f"    - Number of points: {rotor.problem.number_points}")
    print(f"    - Number of harmonics (N): {rotor.problem.number_harmonics}")

    # Calculate HD matrix size
    state_size = 2 * rotor.problem.number_blades+4
    hd_matrix_size = (1 + 2 * rotor.problem.number_harmonics) * state_size
    print(f"    - State space dimension: {state_size}")
    print(f"    - HD matrix dimension: {hd_matrix_size}×{hd_matrix_size}")

    # Step 2: Create HD analysis object
    print("\n[2/4] Initializing HD stability analysis...")
    hd = HDStability(rotor)
    print(f"  [OK] HD analysis object created")

    # Step 3: Run HD analysis over full RPM range
    print("\n[3/4] Running HD analysis...")
    hd = hd.hd_full_range()

    print(f"\n  [OK] HD analysis completed")
    print(f"    - Number of solution points: {len(hd.modal_solution)}")

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

    # Plot 1: Damping (generic - all points)
    fig1 = plotter.plot_damping_generic(
        hd.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(damping_min, damping_max)
    )
    plt.figure(fig1.number)
    plt.title(f"HD Analysis: Damping vs RPM (N={rotor.problem.number_harmonics} harmonics)")
    plt.tight_layout()

    # Plot 2: Frequency (generic - all points)
    fig2 = plotter.plot_frequency_generic(
        hd.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(freq_min, freq_max)
    )
    plt.figure(fig2.number)
    plt.title(f"HD Analysis: Frequency vs RPM (N={rotor.problem.number_harmonics} harmonics)")
    plt.tight_layout()

    print("\n" + "="*60)
    print("HD Analysis Complete!")
    print("="*60)
    # Save figures at 600 DPI

    fig1.savefig("results_figures/hd_damping.png", dpi=600, bbox_inches='tight')

    fig2.savefig("results_figures/hd_frequency.png", dpi=600, bbox_inches='tight')

    print("\n  Figures saved to results_figures/")


    plt.show()


if __name__ == "__main__":
    main()
