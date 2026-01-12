"""Main file to test LTI stability analysis with base rotor configuration.

This script performs a complete LTI stability analysis on the base rotor
without using continuation methods, and generates damping and frequency plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from src.rotor_build import RotorBuild
from src.stability_analysis.lti_stability import LTIStability
from src.utils.plotting import MyPlot


def main():
    """Run LTI stability analysis and plot results."""

    print("="*60)
    print("LTI Stability Analysis - Base Rotor Configuration")
    print("="*60)

    # Step 1: Build rotor system
    print("\n[1/4] Building rotor system...")
    rotor = RotorBuild.build_all()

    # Override to ensure LTI solver is used
    rotor.problem.required_solver = "LTI"
    rotor.problem.continuation = "NO"

    print(f"  [OK] Rotor system built successfully")
    print(f"    - Number of blades: {rotor.problem.number_blades}")
    print(f"    - Solver type: {rotor.problem.required_solver}")
    print(f"    - Damper connection: {rotor.problem.damper_connection}")
    print(f"    - RPM range: {rotor.problem.lower_rotor_RPM} - {rotor.problem.higher_rotor_RPM}")
    print(f"    - Number of points: {rotor.problem.number_points}")

    # Step 2: Create LTI stability analysis object
    print("\n[2/4] Initializing LTI stability analysis...")
    lti_analysis = LTIStability(rotor)
    print(f"  [OK] LTI analysis object created")

    # Step 3: Run eigenvalue analysis over full RPM range
    print("\n[3/4] Computing eigenvalues over RPM range...")
    print(f"  This may take a moment...")
    lti_analysis = lti_analysis.eigen_full_range()

    print(f"  [OK] Eigenvalue analysis completed")
    print(f"    - Number of solution points: {len(lti_analysis.modal_solution)}")
    print(f"    - Number of eigenvalues per point: {len(lti_analysis.modal_solution[0].damping)}")

    # Print sample results
    print(f"\n  Sample results at first operating point:")
    print(f"    - RPM: {lti_analysis.modal_solution[0].OMEGA_RPM:.2f}")
    print(f"    - First 6 damping values: {lti_analysis.modal_solution[0].damping[:6]}")
    print(f"    - First 6 frequency values: {lti_analysis.modal_solution[0].frequency[:6]}")

    # Check stability
    max_damping_overall = max([np.max(sol.damping) for sol in lti_analysis.modal_solution])
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

    all_damping = np.concatenate([sol.damping for sol in lti_analysis.modal_solution])
    all_frequency = np.concatenate([sol.frequency for sol in lti_analysis.modal_solution])

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

    # Determine number of modes from the system
    n_modes = len(lti_analysis.modal_solution[0].damping)
    print(f"\n  Total number of modes: {n_modes}")

    # Plot all modes (using 1-based numbering)
    modes_to_plot_damping = list(range(1, n_modes + 1))

    # Identify modes with positive frequencies at any operating point
    modes_with_positive_freq = set()
    for sol in lti_analysis.modal_solution:
        for i, freq in enumerate(sol.frequency):
            if freq > 0:
                modes_with_positive_freq.add(i + 1)  # Convert to 1-based

    modes_to_plot_frequency = sorted(list(modes_with_positive_freq))

    # Plot 1: Damping (generic - dots only, no lines)
    fig1 = plotter.plot_damping_generic(
        lti_analysis.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(damping_min, damping_max)
    )
    plt.figure(fig1.number)
    plt.title("Damping vs RPM")
    plt.tight_layout()

    # Plot 2: Frequency (generic - dots only, no lines)
    fig2 = plotter.plot_frequency_generic(
        lti_analysis.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(freq_min, freq_max)
    )
    plt.figure(fig2.number)
    plt.title("Frequency vs RPM")
    plt.tight_layout()

    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)
    plt.show()


if __name__ == "__main__":
    main()
