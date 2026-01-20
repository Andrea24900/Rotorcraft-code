# this file contains tests to run, generate plots and save them in a special folder "extra-figures"

# Instructions:

#run the test and generate the plots which compare the C1/Cd ratio for the two methods (as the last example in the main_parametric_analysis_test)
#generate in this file the call for all these tests
# separate clearly the tests

# common plot characteristics:
# no title
# same options as in the test example for the parametric analysis


# tests to run with parametric analysis
# H2B configuration with 4 blades (1 harmonic), 5 blades (2 harmonics) and 7 blades (3 harmonics)
# H2B configuration with 4 blades (1 harmonic) with custom damper activation: two dampers (in adjacent and opposite position in the rotor) fail with a C1 and C2 decreasing equally

import os
import numpy as np
import matplotlib.pyplot as plt
from src.config.problem_definition import ProblemDefinition, create_default_problem
from src.analysis.parameter_analysis import ParametricAnalysis
from src.utils.plotting import MyPlot

# Output directory for figures
FIGURES_DIR = "extra_figures"


# =============================================================================
# TEST 1: H2B configuration with 4 blades (1 harmonic)
# =============================================================================
def test_1_h2b_4blades_1harmonic():
    """H2B configuration with 4 blades, 1 harmonic - C1/Cd ratio comparison (HD vs LTP)."""

    print("\n" + "=" * 70)
    print("TEST 1: H2B - 4 blades, 1 harmonic")
    print("=" * 70)

    # Create problem with CUSTOM damper activation
    problem_hd = create_default_problem()
    problem_hd.damper_connection = "H2B"
    problem_hd.number_blades = 4
    problem_hd.number_harmonics = 1
    problem_hd.damper_activation = "CUSTOM"
    problem_hd.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0])

    print(f"Configuration:")
    print(f"  - Number of blades: {problem_hd.number_blades}")
    print(f"  - Number of harmonics: {problem_hd.number_harmonics}")
    print(f"  - Damper connection: {problem_hd.damper_connection}")

    # Define parameter setter for first damper ratio
    def set_first_damper_ratio(problem: ProblemDefinition, ratio: float):
        """Set the first damper to a fraction of nominal, others remain at 100%."""
        problem.damper_activation = "CUSTOM"
        problem.damper_activation_vector = np.array([ratio, 1.0, 1.0, 1.0])

    # Sweep parameters
    omega_sweep = np.linspace(150, 350, 30)
    damper_ratios = np.linspace(0.0, 1.0, 30)

    # Run HD analysis
    print(f"\n  Running HD solver sweep...")
    analysis_hd = ParametricAnalysis(problem_hd)
    analysis_hd.sweep(
        solver_type="HD",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_1_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_first_damper_ratio,
        verbose=True
    )

    # Run LTP analysis for comparison
    print(f"\n  Running LTP solver sweep...")
    problem_ltp = create_default_problem()
    problem_ltp.damper_connection = "H2B"
    problem_ltp.number_blades = 4
    problem_ltp.number_harmonics = 1
    problem_ltp.damper_activation = "CUSTOM"
    problem_ltp.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0])

    analysis_ltp = ParametricAnalysis(problem_ltp)
    analysis_ltp.sweep(
        solver_type="LTP",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_1_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_first_damper_ratio,
        verbose=True
    )

    # Plot stability map with HD contours and LTP boundary overlay
    print("\n  Creating stability map...")
    fig = MyPlot.plot_parametric_stability_map(
        analysis_hd,
        figsize=(8, 6),
        reference_analysis=analysis_ltp,
        reference_label="LTP stability boundary",
        reference_color="blue",
        reference_linestyle="-",
        title="",
        ylabel=r"$C_1/C_d$ $\mathrm{[-]}$"
    )
    fig.tight_layout()

    # Save figure (PNG and PDF)
    fig_path_png = os.path.join(FIGURES_DIR, "test1_h2b_4blades_1harmonic.png")
    fig_path_pdf = os.path.join(FIGURES_DIR, "test1_h2b_4blades_1harmonic.pdf")
    fig.savefig(fig_path_png, dpi=600, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    # print(f"  Saved: {fig_path_png}")
    # print(f"  Saved: {fig_path_pdf}")

    return fig


# =============================================================================
# TEST 2: H2B configuration with 5 blades (2 harmonics)
# =============================================================================
def test_2_h2b_5blades_2harmonics():
    """H2B configuration with 5 blades, 2 harmonics - C1/Cd ratio comparison (HD vs LTP)."""

    print("\n" + "=" * 70)
    print("TEST 2: H2B - 5 blades, 2 harmonics")
    print("=" * 70)

    # Create problem with CUSTOM damper activation
    problem_hd = create_default_problem()
    problem_hd.damper_connection = "H2B"
    problem_hd.number_blades = 5
    problem_hd.number_harmonics = 2
    problem_hd.damper_activation = "CUSTOM"
    problem_hd.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0, 1.0])

    print(f"Configuration:")
    print(f"  - Number of blades: {problem_hd.number_blades}")
    print(f"  - Number of harmonics: {problem_hd.number_harmonics}")
    print(f"  - Damper connection: {problem_hd.damper_connection}")

    # Define parameter setter for first damper ratio
    def set_first_damper_ratio(problem: ProblemDefinition, ratio: float):
        """Set the first damper to a fraction of nominal, others remain at 100%."""
        problem.damper_activation = "CUSTOM"
        problem.damper_activation_vector = np.array([ratio, 1.0, 1.0, 1.0, 1.0])

    # Sweep parameters
    omega_sweep = np.linspace(150, 350, 30)
    damper_ratios = np.linspace(0.0, 1.0, 30)

    # Run HD analysis
    print(f"\n  Running HD solver sweep...")
    analysis_hd = ParametricAnalysis(problem_hd)
    analysis_hd.sweep(
        solver_type="HD",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_1_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_first_damper_ratio,
        verbose=True
    )

    # Run LTP analysis for comparison
    print(f"\n  Running LTP solver sweep...")
    problem_ltp = create_default_problem()
    problem_ltp.damper_connection = "H2B"
    problem_ltp.number_blades = 5
    problem_ltp.number_harmonics = 2
    problem_ltp.damper_activation = "CUSTOM"
    problem_ltp.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0, 1.0])

    analysis_ltp = ParametricAnalysis(problem_ltp)
    analysis_ltp.sweep(
        solver_type="LTP",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_1_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_first_damper_ratio,
        verbose=True
    )

    # Plot stability map with HD contours and LTP boundary overlay
    print("\n  Creating stability map...")
    fig = MyPlot.plot_parametric_stability_map(
        analysis_hd,
        figsize=(8, 6),
        reference_analysis=analysis_ltp,
        reference_label="LTP stability boundary",
        reference_color="blue",
        reference_linestyle="-",
        title="",
        ylabel=r"$C_1/C_d$ $\mathrm{[-]}$"
    )
    fig.tight_layout()

    # Save figure (PNG and PDF)
    fig_path_png = os.path.join(FIGURES_DIR, "test2_h2b_5blades_2harmonics.png")
    fig_path_pdf = os.path.join(FIGURES_DIR, "test2_h2b_5blades_2harmonics.pdf")
    fig.savefig(fig_path_png, dpi=600, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    # print(f"  Saved: {fig_path_png}")
    # print(f"  Saved: {fig_path_pdf}")

    return fig


# =============================================================================
# TEST 3: H2B configuration with 7 blades (3 harmonics)
# =============================================================================
def test_3_h2b_7blades_3harmonics():
    """H2B configuration with 7 blades, 3 harmonics - C1/Cd ratio comparison (HD vs LTP)."""

    print("\n" + "=" * 70)
    print("TEST 3: H2B - 7 blades, 3 harmonics")
    print("=" * 70)

    # Create problem with CUSTOM damper activation
    problem_hd = create_default_problem()
    problem_hd.damper_connection = "H2B"
    problem_hd.number_blades = 7
    problem_hd.number_harmonics = 3
    problem_hd.damper_activation = "CUSTOM"
    problem_hd.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    print(f"Configuration:")
    print(f"  - Number of blades: {problem_hd.number_blades}")
    print(f"  - Number of harmonics: {problem_hd.number_harmonics}")
    print(f"  - Damper connection: {problem_hd.damper_connection}")

    # Define parameter setter for first damper ratio
    def set_first_damper_ratio(problem: ProblemDefinition, ratio: float):
        """Set the first damper to a fraction of nominal, others remain at 100%."""
        problem.damper_activation = "CUSTOM"
        problem.damper_activation_vector = np.array([ratio, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    # Sweep parameters
    omega_sweep = np.linspace(150, 350, 30)
    damper_ratios = np.linspace(0.0, 1.0, 30)

    # Run HD analysis
    print(f"\n  Running HD solver sweep...")
    analysis_hd = ParametricAnalysis(problem_hd)
    analysis_hd.sweep(
        solver_type="HD",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_1_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_first_damper_ratio,
        verbose=True
    )

    # Run LTP analysis for comparison
    print(f"\n  Running LTP solver sweep...")
    problem_ltp = create_default_problem()
    problem_ltp.damper_connection = "H2B"
    problem_ltp.number_blades = 7
    problem_ltp.number_harmonics = 3
    problem_ltp.damper_activation = "CUSTOM"
    problem_ltp.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    analysis_ltp = ParametricAnalysis(problem_ltp)
    analysis_ltp.sweep(
        solver_type="LTP",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_1_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_first_damper_ratio,
        verbose=True
    )

    # Plot stability map with HD contours and LTP boundary overlay
    print("\n  Creating stability map...")
    fig = MyPlot.plot_parametric_stability_map(
        analysis_hd,
        figsize=(8, 6),
        reference_analysis=analysis_ltp,
        reference_label="LTP stability boundary",
        reference_color="blue",
        reference_linestyle="-",
        title="",
        ylabel=r"$C_1/C_d$ $\mathrm{[-]}$"
    )
    fig.tight_layout()

    # Save figure (PNG and PDF)
    fig_path_png = os.path.join(FIGURES_DIR, "test3_h2b_7blades_3harmonics.png")
    fig_path_pdf = os.path.join(FIGURES_DIR, "test3_h2b_7blades_3harmonics.pdf")
    fig.savefig(fig_path_png, dpi=600, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    # print(f"  Saved: {fig_path_png}")
    # print(f"  Saved: {fig_path_pdf}")

    return fig


# =============================================================================
# TEST 4: H2B 4 blades - Two ADJACENT dampers failing (C1 and C2 decrease equally)
# =============================================================================
def test_4_h2b_4blades_adjacent_dampers_fail():
    """H2B configuration with 4 blades - Two adjacent dampers (1 and 2) failing equally."""

    print("\n" + "=" * 70)
    print("TEST 4: H2B - 4 blades, Two ADJACENT dampers failing (dampers 1 & 2)")
    print("=" * 70)

    # Create problem with CUSTOM damper activation
    problem_hd = create_default_problem()
    problem_hd.damper_connection = "H2B"
    problem_hd.number_blades = 4
    problem_hd.number_harmonics = 1
    problem_hd.damper_activation = "CUSTOM"
    problem_hd.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0])

    print(f"Configuration:")
    print(f"  - Number of blades: {problem_hd.number_blades}")
    print(f"  - Number of harmonics: {problem_hd.number_harmonics}")
    print(f"  - Damper connection: {problem_hd.damper_connection}")
    print(f"  - Failing dampers: 1 and 2 (adjacent)")

    # Define parameter setter for two adjacent dampers failing equally
    def set_adjacent_dampers_ratio(problem: ProblemDefinition, ratio: float):
        """Set dampers 1 and 2 (adjacent) to a fraction of nominal, others remain at 100%."""
        problem.damper_activation = "CUSTOM"
        problem.damper_activation_vector = np.array([ratio, ratio, 1.0, 1.0])

    # Sweep parameters
    omega_sweep = np.linspace(150, 350, 30)
    damper_ratios = np.linspace(0.0, 1.0, 30)

    # Run HD analysis
    print(f"\n  Running HD solver sweep...")
    analysis_hd = ParametricAnalysis(problem_hd)
    analysis_hd.sweep(
        solver_type="HD",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_12_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_adjacent_dampers_ratio,
        verbose=True
    )

    # Run LTP analysis for comparison
    print(f"\n  Running LTP solver sweep...")
    problem_ltp = create_default_problem()
    problem_ltp.damper_connection = "H2B"
    problem_ltp.number_blades = 4
    problem_ltp.number_harmonics = 1
    problem_ltp.damper_activation = "CUSTOM"
    problem_ltp.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0])

    analysis_ltp = ParametricAnalysis(problem_ltp)
    analysis_ltp.sweep(
        solver_type="LTP",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_12_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_adjacent_dampers_ratio,
        verbose=True
    )

    # Plot stability map with HD contours and LTP boundary overlay
    print("\n  Creating stability map...")
    fig = MyPlot.plot_parametric_stability_map(
        analysis_hd,
        figsize=(8, 6),
        reference_analysis=analysis_ltp,
        reference_label="LTP stability boundary",
        reference_color="blue",
        reference_linestyle="-",
        title="",
        ylabel=r"$C_{1,2}/C_d$ $\mathrm{[-]}$"
    )
    fig.tight_layout()

    # Save figure (PNG and PDF)
    fig_path_png = os.path.join(FIGURES_DIR, "test4_h2b_4blades_adjacent_dampers_fail.png")
    fig_path_pdf = os.path.join(FIGURES_DIR, "test4_h2b_4blades_adjacent_dampers_fail.pdf")
    fig.savefig(fig_path_png, dpi=600, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    # print(f"  Saved: {fig_path_png}")
    # print(f"  Saved: {fig_path_pdf}")

    return fig


# =============================================================================
# TEST 5: H2B 4 blades - Two OPPOSITE dampers failing (C1 and C3 decrease equally)
# =============================================================================
def test_5_h2b_4blades_opposite_dampers_fail():
    """H2B configuration with 4 blades - Two opposite dampers (1 and 3) failing equally."""

    print("\n" + "=" * 70)
    print("TEST 5: H2B - 4 blades, Two OPPOSITE dampers failing (dampers 1 & 3)")
    print("=" * 70)

    # Create problem with CUSTOM damper activation
    problem_hd = create_default_problem()
    problem_hd.damper_connection = "H2B"
    problem_hd.number_blades = 4
    problem_hd.number_harmonics = 1
    problem_hd.damper_activation = "CUSTOM"
    problem_hd.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0])

    print(f"Configuration:")
    print(f"  - Number of blades: {problem_hd.number_blades}")
    print(f"  - Number of harmonics: {problem_hd.number_harmonics}")
    print(f"  - Damper connection: {problem_hd.damper_connection}")
    print(f"  - Failing dampers: 1 and 3 (opposite)")

    # Define parameter setter for two opposite dampers failing equally
    def set_opposite_dampers_ratio(problem: ProblemDefinition, ratio: float):
        """Set dampers 1 and 3 (opposite) to a fraction of nominal, others remain at 100%."""
        problem.damper_activation = "CUSTOM"
        problem.damper_activation_vector = np.array([ratio, 1.0, ratio, 1.0])

    # Sweep parameters
    omega_sweep = np.linspace(150, 350, 30)
    damper_ratios = np.linspace(0.0, 1.0, 30)

    # Run HD analysis
    print(f"\n  Running HD solver sweep...")
    analysis_hd = ParametricAnalysis(problem_hd)
    analysis_hd.sweep(
        solver_type="HD",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_13_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_opposite_dampers_ratio,
        verbose=True
    )

    # Run LTP analysis for comparison
    print(f"\n  Running LTP solver sweep...")
    problem_ltp = create_default_problem()
    problem_ltp.damper_connection = "H2B"
    problem_ltp.number_blades = 4
    problem_ltp.number_harmonics = 1
    problem_ltp.damper_activation = "CUSTOM"
    problem_ltp.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0])

    analysis_ltp = ParametricAnalysis(problem_ltp)
    analysis_ltp.sweep(
        solver_type="LTP",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_13_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_opposite_dampers_ratio,
        verbose=True
    )

    # Plot stability map with HD contours and LTP boundary overlay
    print("\n  Creating stability map...")
    fig = MyPlot.plot_parametric_stability_map(
        analysis_hd,
        figsize=(8, 6),
        reference_analysis=analysis_ltp,
        reference_label="LTP stability boundary",
        reference_color="blue",
        reference_linestyle="-",
        title="",
        ylabel=r"$C_{1,3}/C_d$ $\mathrm{[-]}$"
    )
    fig.tight_layout()

    # Save figure (PNG and PDF)
    fig_path_png = os.path.join(FIGURES_DIR, "test5_h2b_4blades_opposite_dampers_fail.png")
    fig_path_pdf = os.path.join(FIGURES_DIR, "test5_h2b_4blades_opposite_dampers_fail.pdf")
    fig.savefig(fig_path_png, dpi=600, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    # print(f"  Saved: {fig_path_png}")
    # print(f"  Saved: {fig_path_pdf}")

    return fig


def test_1_b2b_4blades_1harmonic():
    """B2B configuration with 4 blades, 1 harmonic - C1/Cd ratio comparison (HD vs LTP)."""

    print("\n" + "=" * 70)
    print("TEST 6: B2B - 4 blades, 1 harmonic")
    print("=" * 70)

    # Create problem with CUSTOM damper activation
    problem_hd = create_default_problem()
    problem_hd.damper_connection = "B2B"
    problem_hd.number_blades = 4
    problem_hd.number_harmonics = 1
    problem_hd.damper_activation = "CUSTOM"
    problem_hd.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0])

    print(f"Configuration:")
    print(f"  - Number of blades: {problem_hd.number_blades}")
    print(f"  - Number of harmonics: {problem_hd.number_harmonics}")
    print(f"  - Damper connection: {problem_hd.damper_connection}")

    # Define parameter setter for first damper ratio
    def set_first_damper_ratio(problem: ProblemDefinition, ratio: float):
        """Set the first damper to a fraction of nominal, others remain at 100%."""
        problem.damper_activation = "CUSTOM"
        problem.damper_activation_vector = np.array([ratio, 1.0, 1.0, 1.0])

    # Sweep parameters
    omega_sweep = np.linspace(150, 350, 30)
    damper_ratios = np.linspace(0.0, 1.0, 30)

    # Run HD analysis
    print(f"\n  Running HD solver sweep...")
    analysis_hd = ParametricAnalysis(problem_hd)
    analysis_hd.sweep(
        solver_type="HD",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_1_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_first_damper_ratio,
        verbose=True
    )

    # Run LTP analysis for comparison
    print(f"\n  Running LTP solver sweep...")
    problem_ltp = create_default_problem()
    problem_ltp.damper_connection = "B2B"
    problem_ltp.number_blades = 4
    problem_ltp.number_harmonics = 1
    problem_ltp.damper_activation = "CUSTOM"
    problem_ltp.damper_activation_vector = np.array([1.0, 1.0, 1.0, 1.0])

    analysis_ltp = ParametricAnalysis(problem_ltp)
    analysis_ltp.sweep(
        solver_type="LTP",
        omega_values_RPM=omega_sweep,
        parameter_name="damper_1_ratio",
        parameter_values=damper_ratios,
        parameter_setter=set_first_damper_ratio,
        verbose=True
    )

    # Plot stability map with HD contours and LTP boundary overlay
    print("\n  Creating stability map...")
    fig = MyPlot.plot_parametric_stability_map(
        analysis_hd,
        figsize=(8, 6),
        reference_analysis=analysis_ltp,
        reference_label="LTP stability boundary",
        reference_color="blue",
        reference_linestyle="-",
        title="",
        ylabel=r"$C_1/C_d$ $\mathrm{[-]}$"
    )
    fig.tight_layout()

    # Save figure (PNG and PDF)
    fig_path_png = os.path.join(FIGURES_DIR, "test1_b2b_4blades_1harmonic.png")
    fig_path_pdf = os.path.join(FIGURES_DIR, "test1_b2b_4blades_1harmonic.pdf")
    fig.savefig(fig_path_png, dpi=600, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    # print(f"  Saved: {fig_path_png}")
    # print(f"  Saved: {fig_path_pdf}")

    return fig

# =============================================================================
# MAIN - Run all tests
# =============================================================================
def main():
    """Run all tests and generate plots."""

    # Create output directory if it doesn't exist
    os.makedirs(FIGURES_DIR, exist_ok=True)

    print("=" * 70)
    print("RUNNING ALL PARAMETRIC ANALYSIS TESTS")
    print(f"Output directory: {FIGURES_DIR}")
    print("=" * 70)

    # Run all tests
    #fig1 = test_1_h2b_4blades_1harmonic()
    #fig2 = test_2_h2b_5blades_2harmonics()
    #fig3 = test_3_h2b_7blades_3harmonics()
    fig6 = test_1_b2b_4blades_1harmonic()
    #broken
    #fig4 = test_4_h2b_4blades_adjacent_dampers_fail()
    #fig5 = test_5_h2b_4blades_opposite_dampers_fail()

    # Summary
    print("\n" + "=" * 70)
    print("ALL TESTS COMPLETE!")
    print("=" * 70)
    print(f"\nGenerated figures in '{FIGURES_DIR}/' (PNG and PDF):")
    print("  - test1_h2b_4blades_1harmonic")
    print("  - test2_h2b_5blades_2harmonics")
    print("  - test3_h2b_7blades_3harmonics")
    print("  - test4_h2b_4blades_adjacent_dampers_fail")
    print("  - test5_h2b_4blades_opposite_dampers_fail")

    plt.show()


if __name__ == "__main__":
    main()
