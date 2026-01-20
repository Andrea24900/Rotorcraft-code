"""LTI stability analysis test."""

import numpy as np
import matplotlib.pyplot as plt
from src.config.problem_definition import create_default_problem
from src.rotor_build import RotorBuild
from src.stability_analysis.lti_stability import LTIStability
from src.utils.plotting import MyPlot


def main():
    """Run LTI stability analysis and plot results."""

    print("=" * 60)
    print("LTI Stability Analysis")
    print("=" * 60)

    # Build rotor system - set parameters BEFORE building
    problem = create_default_problem()
    problem.number_blades = 4
    problem.damper_activation = "ALL"
    problem.required_solver = "LTI"
    problem.continuation = "NO"

    rotor = RotorBuild(problem)

    print(f"  Blades: {problem.number_blades}, Solver: {problem.required_solver}")
    print(f"  RPM range: {problem.lower_rotor_RPM}-{problem.higher_rotor_RPM}, Points: {problem.number_points}")

    # Run analysis
    print("\nRunning eigenvalue analysis...")
    lti_analysis = LTIStability(rotor).eigen_full_range()

    # Check stability
    max_damping = max([np.max(sol.damping) for sol in lti_analysis.modal_solution])
    status = "STABLE" if max_damping < 0 else "UNSTABLE"
    print(f"  Result: {status} (max damping: {max_damping:.4f})")

    # Generate plots
    print("\nGenerating plots...")
    plotter = MyPlot()

    rpm_min, rpm_max = problem.lower_rotor_RPM, problem.higher_rotor_RPM
    all_damping = np.concatenate([sol.damping for sol in lti_analysis.modal_solution])
    all_frequency = np.concatenate([sol.frequency for sol in lti_analysis.modal_solution])

    damping_min = np.min(all_damping) - 0.5
    damping_max = max(1.0, np.max(all_damping) + 0.5)

    positive_freq = all_frequency[all_frequency > 0]
    freq_max = np.max(positive_freq) + 5 if len(positive_freq) > 0 else 50

    fig1 = plotter.plot_damping_generic(
        lti_analysis.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(damping_min, damping_max)
    )
    plt.figure(fig1.number)
    plt.title("LTI Analysis: Damping vs RPM")
    plt.tight_layout()

    fig2 = plotter.plot_frequency_generic(
        lti_analysis.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(0, freq_max)
    )
    plt.figure(fig2.number)
    plt.title("LTI Analysis: Frequency vs RPM")
    plt.tight_layout()

    fig1.savefig("results_figures/lti_damping.png", dpi=600, bbox_inches='tight')
    fig2.savefig("results_figures/lti_frequency.png", dpi=600, bbox_inches='tight')
    print("  Figures saved to results_figures/")

    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    plt.show()


if __name__ == "__main__":
    main()
