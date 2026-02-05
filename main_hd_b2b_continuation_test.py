"""HD stability analysis with continuation method for B2B (blade-to-blade) damper configuration."""

import numpy as np
import matplotlib.pyplot as plt
from src.config.problem_definition import ProblemDefinition, RotorCharacteristics
from src.rotor_build import RotorBuild
from src.stability_analysis.continuation_analysis import ContinuationAnalysis
from src.utils.plotting import MyPlot


def main():
    """Run HD continuation stability analysis for B2B configuration and plot results."""

    print("=" * 60)
    print("HD Continuation Stability Analysis - B2B Configuration")
    print("=" * 60)

    # Build rotor system - must use ProblemDefinition constructor for proper B2B initialization
    problem = ProblemDefinition(
        rotor_characteristics=RotorCharacteristics(
            number_blades=4, damper_connection="B2B", damper_activation="ODI"),
        required_solver="HD",
        number_harmonics=1,
        hd_use_complex=True,
        continuation="YES"
    )

    rotor = RotorBuild(problem)

    state_size = 2 * problem.rotor_characteristics.number_blades + 4
    hd_size = (1 + 2 * problem.number_harmonics) * state_size

    print(f"  Damper connection: {problem.rotor_characteristics.damper_connection}")
    print(f"  Blades: {problem.rotor_characteristics.number_blades}, Harmonics: {problem.number_harmonics}")
    print(f"  HD matrix: {hd_size}x{hd_size}")
    print(f"  RPM range: {problem.lower_rotor_RPM}-{problem.higher_rotor_RPM}, Points: {problem.number_points}")
    print(f"  Continuation tol: {problem.continuation_tolerance}, max iter: {problem.continuation_max_iter}")

    # Run analysis
    print("\nRunning HD continuation analysis for B2B...")
    cont_analysis = ContinuationAnalysis(rotor).continuation()

    # Check stability
    max_damping = max([np.max(sol.damping) for sol in cont_analysis.modal_solution])
    status = "STABLE" if max_damping < 0 else "UNSTABLE"
    print(f"  Result: {status} (max damping: {max_damping:.4f})")

    # Generate plots
    print("\nGenerating plots...")
    plotter = MyPlot()

    rpm_min, rpm_max = problem.lower_rotor_RPM, problem.higher_rotor_RPM
    all_damping = np.concatenate([sol.damping for sol in cont_analysis.modal_solution])
    all_frequency = np.concatenate([sol.frequency for sol in cont_analysis.modal_solution])

    damping_min = np.min(all_damping) - 0.5
    damping_max = max(1.0, np.max(all_damping) + 0.5)

    positive_freq = all_frequency[all_frequency > 0]
    freq_max = np.max(positive_freq) + 5 if len(positive_freq) > 0 else 50

    n_modes = len(cont_analysis.modal_solution[0].damping)
    modes_damping = list(range(1, n_modes + 1))

    modes_freq = set()
    for sol in cont_analysis.modal_solution:
        for i, freq in enumerate(sol.frequency):
            if freq > 0:
                modes_freq.add(i + 1)
    modes_freq = sorted(list(modes_freq))

    fig1 = plotter.plot_damping_order(
        cont_analysis.modal_solution,
        modes=modes_damping,
        xlimits=(rpm_min, rpm_max),
        ylimits=(damping_min, damping_max)
    )
    plt.figure(fig1.number)
    plt.title(f"HD Continuation B2B: Damping vs RPM (N={problem.number_harmonics})")
    plt.tight_layout()

    fig2 = plotter.plot_frequency_order(
        cont_analysis.modal_solution,
        modes=modes_freq,
        xlimits=(rpm_min, rpm_max),
        ylimits=(0, freq_max)
    )
    plt.figure(fig2.number)
    plt.title(f"HD Continuation B2B: Frequency vs RPM (N={problem.number_harmonics})")
    plt.tight_layout()

    fig1.savefig("results_figures/hd_b2b_cont_damping.png", dpi=600, bbox_inches='tight')
    fig2.savefig("results_figures/hd_b2b_cont_frequency.png", dpi=600, bbox_inches='tight')
    print("  Figures saved to results_figures/")

    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    plt.show()


if __name__ == "__main__":
    main()
