"""Example script for parametric stability analysis.

This script demonstrates how to use the ParametricAnalysis class to sweep
a parameter (damper coefficient) across a range of values while computing
the stability at each (OMEGA, parameter) point.

The example shows:
1. Setting up a parametric sweep for damper coefficient
2. Running the analysis with LTI solver (valid for symmetric systems)
3. Plotting stability maps, frequency maps, and damping slices
4. Demonstrating the periodicity check with an invalid configuration
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from src.config.problem_definition import ProblemDefinition, create_default_problem
from src.analysis.parameter_analysis import ParametricAnalysis
from src.utils.plotting import MyPlot

# Output directory for figures
FIGURES_DIR = "results_figures"


def main():
    """Run parametric analysis examples."""

    # Create output directory if it doesn't exist
    os.makedirs(FIGURES_DIR, exist_ok=True)

    print("=" * 70)
    print("Parametric Stability Analysis - Damper Coefficient Sweep")
    print("=" * 70)

    # =========================================================================
    # Example 1: Sweep damper coefficient with LTI solver (symmetric system)
    # =========================================================================
    print("\n[Example 1] Damper coefficient sweep with LTI solver")
    print("-" * 70)

    # Create base problem with all dampers active (symmetric -> LTI valid)
    problem = create_default_problem()
    problem.rotor_characteristics.damper_activation = "ALL"
    problem.rotor_characteristics.number_blades = 1
    problem.number_harmonics = 1

    print(f"Configuration:")
    print(f"  - Number of blades: {problem.rotor_characteristics.number_blades}")
    print(f"  - Damper activation: {problem.rotor_characteristics.damper_activation}")
    print(f"  - Damper connection: {problem.rotor_characteristics.damper_connection}")

    # Check periodicity
    is_periodic = ParametricAnalysis.check_system_periodicity(problem)
    print(f"  - System is periodic: {is_periodic}")
    print(f"  - Valid solvers: {'LTP, HD' if is_periodic else 'LTI, LTP, HD'}")

    # Create parametric analysis object
    analysis = ParametricAnalysis(problem)

    # Define sweep parameters
    omega_range_RPM = np.linspace(50, 400, 40)  # 40 points from 50 to 400 RPM
    damping_range = np.linspace(0,4067, 40)  # 25 points from 500 to 8000 Nms/rad

    # Run the sweep
    print(f"\nRunning parametric sweep...")
    analysis.sweep(
        solver_type="LTI",
        omega_values_RPM=omega_range_RPM,
        parameter_name="nominal_damping_Cd",
        parameter_values=damping_range,
        verbose=True
    )

    # =========================================================================
    # Generate plots using MyPlot class
    # =========================================================================
    print("\n[Generating Plots]")
    print("-" * 70)

    # Plot 1: Stability map (contour plot of damping)
    print("  Creating stability map...")
    fig1 = MyPlot.plot_parametric_stability_map(
        analysis,
        figsize=(12, 8),
        cmap="RdYlGn_r",  # Red = unstable, Green = stable
        show_stability_boundary=True,
        title="Ground Resonance Stability Map\n(Damper Coefficient Variation)"
    )
    fig1.tight_layout()

    # Plot 2: Frequency map
    print("  Creating frequency map...")
    fig2 = MyPlot.plot_parametric_frequency_map(
        analysis,
        figsize=(12, 8),
        cmap="viridis",
        title="Least Damped Mode Frequency\n(Damper Coefficient Variation)"
    )
    fig2.tight_layout()

    # Plot 3: Damping slices (line plots at selected parameter values)
    print("  Creating damping slices...")
    fig3 = MyPlot.plot_parametric_damping_slices(
        analysis,
        param_indices=[0, 6, 12, 18, 24],  # 5 evenly spaced indices
        figsize=(12, 6),
        title="Damping vs Rotor Speed at Selected Damper Coefficients"
    )
    fig3.tight_layout()

    # Plot 4: Stability boundary with filled regions
    print("  Creating stability boundary plot...")
    fig4 = MyPlot.plot_parametric_stability_boundary(
        analysis,
        figsize=(12, 6),
        fill_unstable=True,
        title="Stability Boundary\n(Green = Stable, Red = Unstable)"
    )
    fig4.tight_layout()

    # =========================================================================
    # Example 2: Demonstrate periodicity check error
    # =========================================================================
    print("\n" + "=" * 70)
    print("[Example 2] Demonstrating periodicity check with ODI configuration")
    print("-" * 70)

    # Create problem with One Damper Inoperative (periodic system)
    problem_odi = create_default_problem()
    problem_odi.rotor_characteristics.damper_activation = "ODI"

    print(f"Configuration:")
    print(f"  - Damper activation: {problem_odi.rotor_characteristics.damper_activation}")

    is_periodic_odi = ParametricAnalysis.check_system_periodicity(problem_odi)
    print(f"  - System is periodic: {is_periodic_odi}")

    # Try to validate LTI solver (should fail)
    print(f"\nAttempting to validate LTI solver for periodic system...")
    try:
        ParametricAnalysis.validate_solver_for_problem("LTI", problem_odi, raise_error=True)
    except ValueError as e:
        print(f"  [EXPECTED ERROR] {type(e).__name__}:")
        for line in str(e).split('\n'):
            print(f"    {line}")

    # Show that LTP/HD are valid
    print(f"\nValidating LTP solver for periodic system...")
    is_valid = ParametricAnalysis.validate_solver_for_problem("LTP", problem_odi, raise_error=False)
    print(f"  LTP solver valid: {is_valid}")

    is_valid = ParametricAnalysis.validate_solver_for_problem("HD", problem_odi, raise_error=False)
    print(f"  HD solver valid: {is_valid}")

    # =========================================================================
    # Example 3: Custom parameter sweep with HD vs LTP comparison
    # =========================================================================
    print("\n" + "=" * 70)
    print("[Example 3] Custom damper activation vector sweep (HD vs LTP comparison)")
    print("-" * 70)

    # Create problem with CUSTOM damper activation
    problem_custom = create_default_problem()
    problem_custom.rotor_characteristics.damper_activation = "CUSTOM"
    problem_custom.rotor_characteristics.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0])
    print(f"Configuration:")
    print(f"  - Damper activation: {problem_custom.rotor_characteristics.damper_activation}")
    print(f"  - Initial activation vector: {problem_custom.rotor_characteristics.damper_activation_vector}")

    # Define a custom parameter setter that modifies the first damper's effectiveness
    def set_first_damper_ratio(problem: ProblemDefinition, ratio: float):
        """Set the first damper to a fraction of nominal, others remain at 100%."""
        problem.rotor_characteristics.damper_activation = "CUSTOM"
        problem.rotor_characteristics.damper_activation_vector = np.array([ratio, 1.0, 1.0, 1.0])
        

    # Sweep parameters
    omega_sweep = np.linspace(100, 350, 40)
    damper_ratios = np.linspace(0.0, 1.0, 40)  # C1 varies from 0 to Cd

    # -------------------------------------------------------------------------
    # Run HD analysis
    # -------------------------------------------------------------------------
    print(f"\n[3a] Running HD solver sweep...")
    analysis_hd = ParametricAnalysis(problem_custom)
    analysis_hd.sweep(
        solver_type="HD",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_1_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_first_damper_ratio,
        verbose=True
    )

    # -------------------------------------------------------------------------
    # Run LTP analysis for comparison
    # -------------------------------------------------------------------------
    print(f"\n[3b] Running LTP solver sweep for comparison...")
    problem_custom_ltp = create_default_problem()
    problem_custom_ltp.rotor_characteristics.damper_activation = "CUSTOM"
    problem_custom_ltp.rotor_characteristics.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0])

    analysis_ltp = ParametricAnalysis(problem_custom_ltp)
    analysis_ltp.sweep(
        solver_type="LTP",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_1_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_first_damper_ratio,
        verbose=True
    )

    # -------------------------------------------------------------------------
    # Plot HD results with LTP boundary overlay
    # -------------------------------------------------------------------------
    print("\n  Creating stability map with HD contours and LTP boundary overlay...")
    fig5 = MyPlot.plot_parametric_stability_map(
        analysis_hd,
        figsize=(12, 8),
        reference_analysis=analysis_ltp,
        reference_label="LTP stability boundary",
        reference_color="blue",
        reference_linestyle="-",
        title="",
        ylabel=r"$C_1/C_d$ $\mathrm{[-]}$"
    )
    fig5.tight_layout()

    # =========================================================================
    # Save figures
    # =========================================================================
    print("\n[Saving Figures]")
    print("-" * 70)

    # Save stability map (damper coefficient)
    fig1_path = os.path.join(FIGURES_DIR, "parametric_stability_map_Cd.png")
    fig1.savefig(fig1_path, dpi=600, bbox_inches='tight')
    print(f"  Saved: {fig1_path}")

    # Save effectiveness ratio stability map (HD with LTP boundary)
    fig5_path = os.path.join(FIGURES_DIR, "parametric_stability_map_effectiveness_ratio.png")
    fig5.savefig(fig5_path, dpi=600, bbox_inches='tight')
    print(f"  Saved: {fig5_path}")

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)
    print("\nGenerated figures:")
    print("  Fig 1: Stability map (damper coefficient sweep)")
    print("  Fig 2: Frequency map (damper coefficient sweep)")
    print("  Fig 3: Damping slices at selected parameter values")
    print("  Fig 4: Stability boundary with filled regions")
    print("  Fig 5: Stability map (first damper effectiveness sweep)")

    print(f"\nSaved figures to '{FIGURES_DIR}/':")
    print(f"  - parametric_stability_map_Cd.png")
    print(f"  - parametric_stability_map_effectiveness_ratio.png")

    # Find critical damping at a specific RPM
    critical_rpm = 250.0
    critical_damping = analysis.find_stability_boundary(critical_rpm)
    if critical_damping is not None:
        print(f"\nCritical damper coefficient at {critical_rpm} RPM: {critical_damping:.1f} Nms/rad")
    else:
        print(f"\nNo stability boundary crossing found at {critical_rpm} RPM")

    plt.show()


if __name__ == "__main__":
    main()
