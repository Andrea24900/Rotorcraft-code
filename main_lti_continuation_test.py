"""Main file to test LTI stability analysis with continuation method.

This script performs a complete LTI stability analysis using the continuation
method for tracking eigenvalue branches, avoiding mode crossing issues.
"""

import numpy as np
import matplotlib.pyplot as plt
from src.rotor_build import RotorBuild
from src.stability_analysis.continuation_analysis import ContinuationAnalysis
from src.utils.plotting import MyPlot


def main():
    """Run LTI continuation stability analysis and plot results."""

    print("="*60)
    print("LTI Continuation Stability Analysis")
    print("="*60)

    # Step 1: Build rotor system
    print("\n[1/4] Building rotor system...")
    rotor = RotorBuild.build_all()

    # Override to ensure LTI solver with continuation is used
    rotor.problem.required_solver = "LTI"
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
    print("\n[2/4] Initializing continuation analysis...")
    cont_analysis = ContinuationAnalysis(rotor)
    print(f"  [OK] Continuation analysis object created")

    # Step 3: Run continuation analysis over full RPM range
    print("\n[3/4] Running continuation analysis...")
    print(f"  This may take several minutes due to iterative corrections...")
    cont_analysis = cont_analysis.continuation()

    print(f"\n  [OK] Continuation analysis completed")
    print(f"    - Number of solution points: {len(cont_analysis.modal_solution)}")
    print(f"    - Number of eigenvalues per point: {len(cont_analysis.modal_solution[0].damping)}")

    # Print sample results
    print(f"\n  Sample results at first operating point:")
    print(f"    - RPM: {cont_analysis.modal_solution[0].OMEGA_RPM:.2f}")
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

    # For 4-blade rotor, we expect 12 eigenvalues
    n_modes = len(cont_analysis.modal_solution[0].damping)
    print(f"\n  Total number of modes: {n_modes}")

    # Plot all modes for damping (using 1-based numbering)
    modes_to_plot_damping = list(range(1, min(13, n_modes + 1)))  # Modes 1-12

    # For frequency, only plot modes with positive frequencies
    modes_with_positive_freq = set()
    for sol in cont_analysis.modal_solution:
        for i, freq in enumerate(sol.frequency):
            if freq > 0 and i < 12:
                modes_with_positive_freq.add(i + 1)  # Convert to 1-based

    modes_to_plot_frequency = sorted(list(modes_with_positive_freq))

    print(f"  Modes with positive frequencies: {modes_to_plot_frequency}")
    print(f"  Plotting all dampings for modes: {modes_to_plot_damping}")

    # Plot 1: Damping by mode order (all 12 modes) with continuation tracking
    if len(modes_to_plot_damping) > 0:
        print(f"\n  - Creating damping plot for all {len(modes_to_plot_damping)} modes...")
        print(f"    Note: Continuation method tracks individual eigenvalue branches")
        fig1 = plotter.plot_damping_order(
            cont_analysis.modal_solution,
            modes=modes_to_plot_damping,
            xlimits=(rpm_min, rpm_max),
            ylimits=(damping_min, damping_max)
        )
        plt.figure(fig1.number)
        plt.title("Damping vs RPM (Continuation Method)")
        plt.tight_layout()

    # Plot 2: Frequency by mode order (only modes with positive frequencies)
    if len(modes_to_plot_frequency) > 0:
        print(f"  - Creating frequency plot for {len(modes_to_plot_frequency)} modes with positive frequencies...")
        fig2 = plotter.plot_frequency_order(
            cont_analysis.modal_solution,
            modes=modes_to_plot_frequency,
            xlimits=(rpm_min, rpm_max),
            ylimits=(freq_min, freq_max)
        )
        plt.figure(fig2.number)
        plt.title("Frequency vs RPM (Continuation Method)")
        plt.tight_layout()

    print("\n" + "="*60)
    print("Continuation Analysis Complete!")
    print("="*60)
    print("\nThe continuation method provides smooth eigenvalue tracking")
    print("across the parameter range, avoiding mode crossing issues.")

    # Show plots
    print("\nDisplaying plots...")
    plt.show()


if __name__ == "__main__":
    main()
