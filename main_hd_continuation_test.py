"""Main file to test HD stability analysis with continuation method.

This script performs HD (Harmonic Decomposition) stability analysis using
the continuation method for tracking eigenvalue branches across parameter space.

HD Continuation Method:
-----------------------
Combines the Harmonic Decomposition approach with predictor-corrector algorithms:

1. HD Method: Expands periodic system A(t) into Hill matrix H
   - A(t) = A₀ + Σ[Aₖᶜcos(kωt) + Aₖˢsin(kωt)]
   - Results in (1+2N)×n dimensional time-invariant system
   - Direct eigenvalue computation

2. Continuation: Tracks eigenvalue branches smoothly
   - Predictor: Uses sensitivity dH/dΩ to extrapolate
   - Corrector: Newton-Raphson to find actual eigenvalue
   - Smooth mode tracking across RPM range

Benefits:
---------
- No numerical integration (unlike LTP continuation)
- Direct sensitivity computation from Hill matrix
- Larger but structured eigenvalue problems
- Smooth tracking of complex eigenvalue branches
- Essential for bifurcation analysis in frequency domain

Computational Cost:
-------------------
- Standard HD: O((1+2N)³·n³) per point
- HD Continuation: O((1+2N)³·n³ + (1+2N)²·n²·iter) per point
- Trade-off: ~2-3x time for robust branch tracking

Matrix Dimensions:
------------------
- HD matrix size: (1+2N)×n where N is number of harmonics (user-defined), n is state dimension
- Computational cost scales with ((1+2N)×n)³ for eigenvalue problem
- Sensitivity matrix has same dimensions as HD matrix

Usage:
------
python main_hd_continuation_test.py

Output:
-------
- Console: Detailed progress and convergence information
- Plots: Damping and frequency with smooth mode tracking
"""

import numpy as np
import matplotlib.pyplot as plt
from src.rotor_build import RotorBuild
from src.stability_analysis.continuation_analysis import ContinuationAnalysis
from src.utils.plotting import MyPlot


def main():
    """Run HD continuation stability analysis and plot results.

    This function performs a complete HD continuation analysis workflow:
    1. Build rotor system with HD solver and continuation enabled
    2. Initialize continuation analysis object
    3. Compute eigenvalues over RPM range using HD continuation
    4. Display results and stability assessment
    5. Generate and display plots with smooth mode tracking

    The HD continuation method combines frequency-domain expansion with
    predictor-corrector algorithms for robust eigenvalue tracking.

    Outputs:
    --------
    - Console output with progress and convergence information
    - Damping plot: Real parts with smooth line tracking
    - Frequency plot: Imaginary parts with smooth line tracking

    Notes:
    ------
    - Continuation tolerance can be adjusted in problem definition
    - Higher harmonics increase accuracy but also matrix size
    - Smooth tracking essential for identifying bifurcations
    """

    print("="*60)
    print("HD Continuation Stability Analysis")
    print("="*60)

    # Step 1: Build rotor system
    print("\n[1/4] Building rotor system...")
    rotor = RotorBuild.build_all()

    # Override to ensure HD solver with continuation is used
    rotor.problem.required_solver = "HD"
    rotor.problem.continuation = "YES"
    rotor.problem.number_harmonics = 1  # Start with 1 harmonic

    print(f"  [OK] Rotor system built successfully")
    print(f"    - Number of blades: {rotor.problem.number_blades}")
    print(f"    - Solver type: {rotor.problem.required_solver}")
    print(f"    - Continuation: {rotor.problem.continuation}")
    print(f"    - Damper connection: {rotor.problem.damper_connection}")
    print(f"    - RPM range: {rotor.problem.lower_rotor_RPM} - {rotor.problem.higher_rotor_RPM}")
    print(f"    - Number of points: {rotor.problem.number_points}")
    print(f"    - Number of harmonics (N): {rotor.problem.number_harmonics}")
    print(f"    - Continuation tolerance: {rotor.problem.continuation_tolerance}")
    print(f"    - Max iterations: {rotor.problem.continuation_max_iter}")

    # Calculate HD matrix size
    state_size = 2 * rotor.problem.number_blades
    hd_matrix_size = (1 + 2 * rotor.problem.number_harmonics) * state_size
    print(f"    - State space dimension: {state_size}")
    print(f"    - HD matrix dimension: {hd_matrix_size}×{hd_matrix_size}")

    # Step 2: Create continuation analysis object
    print("\n[2/4] Initializing HD continuation analysis...")
    cont_analysis = ContinuationAnalysis(rotor)
    print(f"  [OK] Continuation analysis object created")

    # Step 3: Run continuation analysis over full RPM range
    print("\n[3/4] Running HD continuation analysis...")
    cont_analysis = cont_analysis.continuation()

    print(f"\n  [OK] HD continuation analysis completed")
    print(f"    - Number of solution points: {len(cont_analysis.modal_solution)}")

    # Check stability
    max_damping_overall = max([np.max(sol.damping) for sol in cont_analysis.modal_solution])
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

    all_damping = np.concatenate([sol.damping for sol in cont_analysis.modal_solution])
    all_frequency = np.concatenate([sol.frequency for sol in cont_analysis.modal_solution])

    # Damping limits - include all dampings
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

    # Get number of eigenvalues
    n_eigenvalues = len(cont_analysis.modal_solution[0].damping)

    # Plot all modes for damping
    modes_to_plot_damping = list(range(1, n_eigenvalues + 1))

    # For frequency, only plot modes with positive frequencies
    modes_with_positive_freq = set()
    for sol in cont_analysis.modal_solution:
        for i, freq in enumerate(sol.frequency):
            if freq > 0:
                modes_with_positive_freq.add(i + 1)  # Convert to 1-based

    modes_to_plot_frequency = sorted(list(modes_with_positive_freq))

    # Plot 1: Damping by mode order with continuation tracking
    if len(modes_to_plot_damping) > 0:
        fig1 = plotter.plot_damping_order(
            cont_analysis.modal_solution,
            modes=modes_to_plot_damping,
            xlimits=(rpm_min, rpm_max),
            ylimits=(damping_min, damping_max)
        )
        plt.figure(fig1.number)
        plt.title(f"Damping vs RPM (HD Continuation, N={rotor.problem.number_harmonics})")
        plt.tight_layout()

    # Plot 2: Frequency by mode order
    if len(modes_to_plot_frequency) > 0:
        fig2 = plotter.plot_frequency_order(
            cont_analysis.modal_solution,
            modes=modes_to_plot_frequency,
            xlimits=(rpm_min, rpm_max),
            ylimits=(freq_min, freq_max)
        )
        plt.figure(fig2.number)
        plt.title(f"Frequency vs RPM (HD Continuation, N={rotor.problem.number_harmonics})")
        plt.tight_layout()

    print("\n" + "="*60)
    print("HD Continuation Analysis Complete!")
    print("="*60)
    plt.show()


if __name__ == "__main__":
    main()
