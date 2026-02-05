"""HD (Harmonic Decomposition) stability analysis test."""

import numpy as np
import matplotlib.pyplot as plt
from src.config.problem_definition import create_default_problem
from src.rotor_build import RotorBuild
from src.stability_analysis.hd_stability import HDStability
from src.utils.plotting import MyPlot


def main():
    """Run HD stability analysis and plot results."""

    print("=" * 60)
    print("HD (Harmonic Decomposition) Stability Analysis")
    print("=" * 60)

    # Build rotor system - set parameters BEFORE building
    problem = create_default_problem()
    problem.rotor_characteristics.number_blades = 4
    problem.rotor_characteristics.damper_activation = "ODI"
    problem.required_solver = "HD"
    problem.number_harmonics = 1

    rotor = RotorBuild(problem)

    state_size = 2 * problem.rotor_characteristics.number_blades + 4
    hd_size = (1 + 2 * problem.number_harmonics) * state_size

    print(f"  Blades: {problem.rotor_characteristics.number_blades}, Harmonics: {problem.number_harmonics}")
    print(f"  HD matrix: {hd_size}x{hd_size}")
    print(f"  RPM range: {problem.lower_rotor_RPM}-{problem.higher_rotor_RPM}, Points: {problem.number_points}")

    # Run analysis
    print("\nRunning HD analysis...")
    hd_analysis = HDStability(rotor).hd_full_range()

    # Check stability
    max_damping = max([np.max(sol.damping) for sol in hd_analysis.modal_solution])
    status = "STABLE" if max_damping < 0 else "UNSTABLE"
    print(f"  Result: {status} (max damping: {max_damping:.4f})")

    # Generate plots
    print("\nGenerating plots...")
    plotter = MyPlot()

    rpm_min, rpm_max = problem.lower_rotor_RPM, problem.higher_rotor_RPM
    all_damping = np.concatenate([sol.damping for sol in hd_analysis.modal_solution])
    all_frequency = np.concatenate([sol.frequency for sol in hd_analysis.modal_solution])

    damping_min = np.min(all_damping) - 0.5
    damping_max = max(1.0, np.max(all_damping) + 0.5)

    positive_freq = all_frequency[all_frequency > 0]
    freq_max = np.max(positive_freq) + 5 if len(positive_freq) > 0 else 50

    fig1 = plotter.plot_damping_generic(
        hd_analysis.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(damping_min, damping_max)
    )
    plt.figure(fig1.number)
    plt.title(f"HD Analysis: Damping vs RPM (N={problem.number_harmonics})")
    plt.tight_layout()

    fig2 = plotter.plot_frequency_generic(
        hd_analysis.modal_solution,
        xlimits=(rpm_min, rpm_max),
        ylimits=(0, freq_max)
    )
    plt.figure(fig2.number)
    plt.title(f"HD Analysis: Frequency vs RPM (N={problem.number_harmonics})")
    plt.tight_layout()

    fig1.savefig("results_figures/hd_damping.png", dpi=600, bbox_inches='tight')
    fig2.savefig("results_figures/hd_frequency.png", dpi=600, bbox_inches='tight')
    print("  Figures saved to results_figures/")

    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    plt.show()


if __name__ == "__main__":
    main()
