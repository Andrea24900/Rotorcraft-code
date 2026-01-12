"""Main file to test LTP stability analysis with base rotor configuration.

This script performs a complete Linear Time-Periodic (LTP) stability analysis
on the base rotor system without using continuation methods. It uses Floquet
theory to compute characteristic exponents and generates damping/frequency plots.

LTP Stability Analysis:
----------------------
LTP analysis is used for rotor systems with time-periodic coefficients, such as:
- Multi-blade systems with blade passing effects
- Rotors with asymmetric properties (anisotropy)
- Systems with periodic external excitation

The analysis uses Floquet theory, which characterizes stability through:
1. Monodromy Matrix (M): The state transition matrix over one period
2. Characteristic Multipliers (μ): Eigenvalues of M
3. Characteristic Exponents (λ): Related to multipliers by μ = e^(λT)
   - Real part (σ): Indicates stability (σ < 0 → stable, σ > 0 → unstable)
   - Imaginary part (ω): Indicates oscillation frequency

Comparison with LTI:
-------------------
- LTI (Linear Time-Invariant): For systems with constant coefficients
  → Uses standard eigenvalue analysis on A matrix
- LTP (Linear Time-Periodic): For systems with periodic coefficients
  → Uses monodromy matrix and Floquet theory
- LTP reduces to LTI when A(t) = constant

Usage:
-----
python main_ltp_test.py

Output:
------
- Console: Analysis progress and key results
- Plots: Damping and frequency vs RPM diagrams

Requirements:
------------
- Base rotor configuration with LTP solver settings
- Sufficient computation time (monodromy matrix integration is expensive)
"""

import numpy as np
import matplotlib.pyplot as plt
from src.rotor_build import RotorBuild
from src.stability_analysis.ltp_stability import LTPStability
from src.utils.plotting import MyPlot


def main():
    """Run LTP stability analysis and plot results.

    This function performs a complete LTP stability analysis workflow:
    1. Build rotor system with LTP solver configuration
    2. Initialize LTP stability analysis object
    3. Compute characteristic exponents over RPM range using Floquet theory
    4. Display results and stability assessment
    5. Generate and display damping and frequency plots

    The analysis computes the monodromy matrix at each RPM point by
    integrating the fundamental solution matrix over one period. This
    is computationally expensive but provides accurate stability
    characterization for time-periodic systems.

    Outputs:
    -------
    - Console output with progress information and key results
    - Damping plot: Real part of characteristic exponents vs RPM
    - Frequency plot: Imaginary part of characteristic exponents vs RPM

    Notes:
    -----
    - Computation time scales with number of points and system size
    - Stability is determined by the maximum real characteristic exponent (damping)
    """

    print("="*60)
    print("LTP Stability Analysis - Base Rotor Configuration")
    print("="*60)

    # Step 1: Build rotor system
    print("\n[1/4] Building rotor system...")
    rotor = RotorBuild.build_all()

    # Override to ensure LTP solver is used
    rotor.problem.required_solver = "LTP"
    rotor.problem.continuation = "NO"

    print(f"  [OK] Rotor system built successfully")
    print(f"    - Number of blades: {rotor.problem.number_blades}")
    print(f"    - Solver type: {rotor.problem.required_solver}")
    print(f"    - Damper connection: {rotor.problem.damper_connection}")
    print(f"    - RPM range: {rotor.problem.lower_rotor_RPM} - {rotor.problem.higher_rotor_RPM}")
    print(f"    - Number of points: {rotor.problem.number_points}")

    # Step 2: Create LTP stability analysis object
    print("\n[2/4] Initializing LTP stability analysis...")
    ltp_analysis = LTPStability(rotor)
    print(f"  [OK] LTP analysis object created")

    # Step 3: Run Floquet analysis over full RPM range
    print("\n[3/4] Computing characteristic exponents over RPM range...")
    print(f"  This may take several minutes (computing monodromy matrices)...")
    ltp_analysis = ltp_analysis.ltp_full_range()

    print(f"  [OK] Floquet analysis completed")
    print(f"    - Number of solution points: {len(ltp_analysis.modal_solution)}")
    print(f"    - Number of characteristic exponents per point: {len(ltp_analysis.modal_solution[0].damping)}")

    # Print sample results
    print(f"\n  Sample results at first operating point:")
    print(f"    - RPM: {ltp_analysis.modal_solution[0].OMEGA_RPM:.2f}")
    print(f"    - Period T: {ltp_analysis.modal_solution[0].T:.6f} s")
    print(f"    - First 6 real characteristic exponents (damping): {ltp_analysis.modal_solution[0].damping[:6]}")
    print(f"    - First 6 imaginary characteristic exponents (frequency): {ltp_analysis.modal_solution[0].frequency[:6]}")

    # Check stability
    max_damping_overall = max([np.max(sol.damping) for sol in ltp_analysis.modal_solution])
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

    all_damping = np.concatenate([sol.damping for sol in ltp_analysis.modal_solution])
    all_frequency = np.concatenate([sol.frequency for sol in ltp_analysis.modal_solution])

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
    n_modes = len(ltp_analysis.modal_solution[0].damping)
    print(f"\n  Total number of modes: {n_modes}")

    # Plot all modes (using 1-based numbering)
    modes_to_plot_damping = list(range(1, min(13, n_modes + 1)))  # Modes 1-12 for damping

    # For frequency, only plot modes with positive frequencies
    # Identify modes with positive frequencies at any operating point
    modes_with_positive_freq = set()
    for sol in ltp_analysis.modal_solution:
        for i, freq in enumerate(sol.frequency):
            if freq > 0 and i < 12:
                modes_with_positive_freq.add(i + 1)  # Convert to 1-based

    modes_to_plot_frequency = sorted(list(modes_with_positive_freq))

    print(f"  Modes with positive frequencies: {modes_to_plot_frequency}")
    print(f"  Plotting all dampings for modes: {modes_to_plot_damping}")

    # Plot 1: Damping (generic - dots only, no lines)
    print(f"\n  - Creating damping plot (all modes, dots only)...")
    fig1 = plotter.plot_damping_generic(
        ltp_analysis.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(damping_min, damping_max)
    )
    plt.figure(fig1.number)
    plt.title("Damping vs RPM (LTP - Floquet Theory)")
    plt.tight_layout()

    # Plot 2: Frequency (generic - dots only, no lines)
    print(f"  - Creating frequency plot (all modes, dots only)...")
    fig2 = plotter.plot_frequency_generic(
        ltp_analysis.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(freq_min, freq_max)
    )
    plt.figure(fig2.number)
    plt.title("Frequency vs RPM (LTP - Floquet Theory)")
    plt.tight_layout()

    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)

    # Show plots
    print("\nDisplaying plots...")
    # Save figures at 600 DPI

    fig1.savefig("results_figures/ltp_damping.png", dpi=600, bbox_inches='tight')

    fig2.savefig("results_figures/ltp_frequency.png", dpi=600, bbox_inches='tight')

    print("\n  Figures saved to results_figures/")


    plt.show()


if __name__ == "__main__":
    main()
