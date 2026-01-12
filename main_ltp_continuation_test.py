"""Main file to test LTP stability analysis with continuation method.

This script performs LTP stability analysis using the continuation method
for tracking Floquet multiplier branches across the parameter range. This
avoids mode crossing issues and provides smooth tracking of characteristic
exponents.

LTP Continuation vs Standard LTP:
---------------------------------
- Standard LTP: Solves monodromy matrix eigenvalue problem independently at each RPM
  → Can have mode crossing/jumping issues
  → Each solution is independent

- LTP Continuation: Uses predictor-corrector to track individual eigenvalue branches
  → Smooth tracking of each Floquet multiplier branch
  → Handles near-degeneracies and avoided crossings
  → More computationally expensive but more robust

Algorithm:
----------
1. Predictor: Use sensitivity analysis to predict next Floquet multipliers
   - Compute dM/dΩ (monodromy matrix sensitivity)
   - Extrapolate eigenvalues linearly

2. Corrector: Newton-Raphson iteration to correct prediction
   - Iterate until convergence within tolerance
   - Find closest eigenvalue to prediction

Usage:
------
python main_ltp_continuation_test.py

Output:
-------
- Console: Analysis progress and key results
- Plots: Damping and frequency vs RPM with smooth mode tracking

Requirements:
-------------
- LTP solver with continuation enabled
- Sufficient computation time (monodromy + sensitivity at each point)
"""

import numpy as np
import matplotlib.pyplot as plt
from src.rotor_build import RotorBuild
from src.stability_analysis.continuation_analysis import ContinuationAnalysis
from src.utils.plotting import MyPlot


def main():
    """Run LTP continuation stability analysis and plot results.

    This function performs a complete LTP continuation analysis workflow:
    1. Build rotor system with LTP solver and continuation enabled
    2. Initialize continuation analysis object
    3. Compute Floquet multipliers over RPM range using continuation
    4. Display results and stability assessment
    5. Generate and display damping and frequency plots with mode tracking

    The continuation method tracks individual Floquet multiplier branches
    smoothly across the parameter space, avoiding mode crossing issues.

    Outputs:
    --------
    - Console output with progress information and key results
    - Damping plot: Real part of characteristic exponents vs RPM (with lines)
    - Frequency plot: Imaginary part of characteristic exponents vs RPM (with lines)

    Notes:
    ------
    - Computation time is higher than standard LTP due to sensitivity calculations
    - Provides smooth eigenvalue tracking essential for bifurcation analysis
    - Tolerance and max iterations can be adjusted in problem definition
    """

    print("="*60)
    print("LTP Continuation Stability Analysis")
    print("="*60)

    # Step 1: Build rotor system
    print("\n[1/4] Building rotor system...")
    rotor = RotorBuild.build_all()

    # Override to ensure LTP solver with continuation is used
    rotor.problem.required_solver = "LTP"
    rotor.problem.continuation = "YES"

    print(f"  [OK] Rotor system built successfully")
    print(f"    - Number of blades: {rotor.problem.number_blades}")
    print(f"    - Solver type: {rotor.problem.required_solver}")
    print(f"    - Continuation: {rotor.problem.continuation}")
    print(f"    - Damper connection: {rotor.problem.damper_connection}")
    print(f"    - RPM range: {rotor.problem.lower_rotor_RPM} - {rotor.problem.higher_rotor_RPM}")
    print(f"    - Number of points: {rotor.problem.number_points}")
    print(f"    - Continuation tolerance: {rotor.problem.continuation_tolerance}")
    print(f"    - Max iterations: {rotor.problem.continuation_max_iter}")

    # Step 2: Create continuation analysis object
    print("\n[2/4] Initializing LTP continuation analysis...")
    cont_analysis = ContinuationAnalysis(rotor)
    print(f"  [OK] Continuation analysis object created")

    # Step 3: Run continuation analysis over full RPM range
    print("\n[3/4] Running LTP continuation analysis...")
    print(f"  Computing monodromy matrices and sensitivities...")
    print(f"  This will take longer than standard LTP due to sensitivity calculations...")
    cont_analysis = cont_analysis.continuation()

    print(f"\n  [OK] LTP continuation analysis completed")
    print(f"    - Number of solution points: {len(cont_analysis.modal_solution)}")
    print(f"    - Number of characteristic exponents per point: {len(cont_analysis.modal_solution[0].damping)}")

    # Print sample results
    print(f"\n  Sample results at first operating point:")
    print(f"    - RPM: {cont_analysis.modal_solution[0].OMEGA_RPM:.2f}")
    print(f"    - Period T: {cont_analysis.modal_solution[0].T:.6f} s")
    print(f"    - First 6 damping values: {cont_analysis.modal_solution[0].damping[:6]}")
    print(f"    - First 6 frequency values: {cont_analysis.modal_solution[0].frequency[:6]}")

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

    # Get number of modes
    n_modes = len(cont_analysis.modal_solution[0].damping)
    print(f"\n  Total number of modes: {n_modes}")

    # Plot all modes for damping (using 1-based numbering)
    modes_to_plot_damping = list(range(1, n_modes + 1))

    # Identify modes with positive frequencies at any operating point
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
        plt.title("Damping vs RPM (LTP Continuation Method)")
        plt.tight_layout()

    # Plot 2: Frequency by mode order (only modes with positive frequencies)
    if len(modes_to_plot_frequency) > 0:
        fig2 = plotter.plot_frequency_order(
            cont_analysis.modal_solution,
            modes=modes_to_plot_frequency,
            xlimits=(rpm_min, rpm_max),
            ylimits=(freq_min, freq_max)
        )
        plt.figure(fig2.number)
        plt.title("Frequency vs RPM (LTP Continuation Method)")
        plt.tight_layout()

    print("\n" + "="*60)
    print("LTP Continuation Analysis Complete!")
    print("="*60)
    print("\nThe LTP continuation method combines Floquet theory with")
    print("predictor-corrector algorithms for robust eigenvalue tracking")
    print("across the parameter range, avoiding mode crossing issues.")

    # Show plots
    print("\nDisplaying plots...")
    plt.show()


if __name__ == "__main__":
    main()
