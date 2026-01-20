"""Tests to generate stability plots for various blade configurations.

Generates plots in extra_figures/ folder (PNG and PDF):
1. 4 blades: LTI, LTP, HD, LTI-cont, LTP-cont, HD-cont
2. 5 blades: LTP (ALL, ODI), HD with 2 harmonics (ALL, ODI)
3. 7 blades: LTP (ALL, ODI), HD with 3 harmonics (ALL, ODI)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from src.config.problem_definition import create_default_problem
from src.rotor_build import RotorBuild
from src.stability_analysis.lti_stability import LTIStability
from src.stability_analysis.ltp_stability import LTPStability
from src.stability_analysis.hd_stability import HDStability
from src.stability_analysis.continuation_analysis import ContinuationAnalysis
from src.utils.plotting import MyPlot

FIGURES_DIR = "extra_figures"


def save_figures(fig1, fig2, base_name):
    """Save damping and frequency figures in PNG and PDF."""
    fig1.savefig(os.path.join(FIGURES_DIR, f"{base_name}_damping.png"), dpi=600, bbox_inches='tight')
    fig1.savefig(os.path.join(FIGURES_DIR, f"{base_name}_damping.pdf"), bbox_inches='tight')
    fig2.savefig(os.path.join(FIGURES_DIR, f"{base_name}_frequency.png"), dpi=600, bbox_inches='tight')
    fig2.savefig(os.path.join(FIGURES_DIR, f"{base_name}_frequency.pdf"), bbox_inches='tight')
    print(f"  Saved: {base_name}_damping.png/pdf, {base_name}_frequency.png/pdf")


def plot_generic(modal_solution, problem, title_prefix):
    """Generate generic damping and frequency plots."""
    plotter = MyPlot()
    rpm_min, rpm_max = problem.lower_rotor_RPM, problem.higher_rotor_RPM

    all_damping = np.concatenate([sol.damping for sol in modal_solution])
    all_frequency = np.concatenate([sol.frequency for sol in modal_solution])

    damping_min = np.min(all_damping) - 0.5
    damping_max = max(1.0, np.max(all_damping) + 0.5)

    positive_freq = all_frequency[all_frequency > 0]
    freq_max = np.max(positive_freq) + 5 if len(positive_freq) > 0 else 50

    fig1 = plotter.plot_damping_generic(modal_solution, xlimits=(rpm_min, rpm_max), ylimits=(damping_min, damping_max))
    plt.figure(fig1.number)
    plt.title(f"{title_prefix}: Damping vs RPM")
    plt.tight_layout()

    fig2 = plotter.plot_frequency_generic(modal_solution, xlimits=(rpm_min, rpm_max), ylimits=(0, freq_max))
    plt.figure(fig2.number)
    plt.title(f"{title_prefix}: Frequency vs RPM")
    plt.tight_layout()

    return fig1, fig2


def plot_continuation(modal_solution, problem, title_prefix):
    """Generate continuation damping and frequency plots with mode tracking."""
    plotter = MyPlot()
    rpm_min, rpm_max = problem.lower_rotor_RPM, problem.higher_rotor_RPM

    all_damping = np.concatenate([sol.damping for sol in modal_solution])
    all_frequency = np.concatenate([sol.frequency for sol in modal_solution])

    damping_min = np.min(all_damping) - 0.5
    damping_max = max(1.0, np.max(all_damping) + 0.5)

    positive_freq = all_frequency[all_frequency > 0]
    freq_max = np.max(positive_freq) + 5 if len(positive_freq) > 0 else 50

    n_modes = len(modal_solution[0].damping)
    modes_damping = list(range(1, n_modes + 1))

    modes_freq = set()
    for sol in modal_solution:
        for i, freq in enumerate(sol.frequency):
            if freq > 0:
                modes_freq.add(i + 1)
    modes_freq = sorted(list(modes_freq))

    fig1 = plotter.plot_damping_order(modal_solution, modes=modes_damping, xlimits=(rpm_min, rpm_max), ylimits=(damping_min, damping_max))
    plt.figure(fig1.number)
    plt.title(f"{title_prefix}: Damping vs RPM")
    plt.tight_layout()

    fig2 = plotter.plot_frequency_order(modal_solution, modes=modes_freq, xlimits=(rpm_min, rpm_max), ylimits=(0, freq_max))
    plt.figure(fig2.number)
    plt.title(f"{title_prefix}: Frequency vs RPM")
    plt.tight_layout()

    return fig1, fig2


# =============================================================================
# 4 BLADES TESTS
# =============================================================================

def test_4blades_lti():
    """4 blades - LTI analysis."""
    print("\n[4 blades] LTI analysis...")
    problem = create_default_problem()
    problem.number_blades = 4
    problem.damper_activation = "ALL"
    problem.required_solver = "LTI"

    rotor = RotorBuild(problem)
    analysis = LTIStability(rotor).eigen_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "4 blades LTI")
    save_figures(fig1, fig2, "4blades_lti")
    return fig1, fig2


def test_4blades_ltp():
    """4 blades - LTP analysis (ODI)."""
    print("\n[4 blades] LTP analysis (ODI)...")
    problem = create_default_problem()
    problem.number_blades = 4
    problem.damper_activation = "ODI"
    problem.required_solver = "LTP"

    rotor = RotorBuild(problem)
    analysis = LTPStability(rotor).ltp_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "4 blades LTP (ODI)")
    save_figures(fig1, fig2, "4blades_ltp_odi")
    return fig1, fig2


def test_4blades_hd():
    """4 blades - HD analysis (ODI, 1 harmonic)."""
    print("\n[4 blades] HD analysis (ODI, N=1)...")
    problem = create_default_problem()
    problem.number_blades = 4
    problem.damper_activation = "ODI"
    problem.required_solver = "HD"
    problem.number_harmonics = 1

    rotor = RotorBuild(problem)
    analysis = HDStability(rotor).hd_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "4 blades HD (ODI, N=1)")
    save_figures(fig1, fig2, "4blades_hd_odi")
    return fig1, fig2


def test_4blades_lti_cont():
    """4 blades - LTI continuation analysis."""
    print("\n[4 blades] LTI continuation analysis...")
    problem = create_default_problem()
    problem.number_blades = 4
    problem.damper_activation = "ALL"
    problem.required_solver = "LTI"
    problem.continuation = "YES"

    rotor = RotorBuild(problem)
    analysis = ContinuationAnalysis(rotor).continuation()

    fig1, fig2 = plot_continuation(analysis.modal_solution, problem, "4 blades LTI Cont")
    save_figures(fig1, fig2, "4blades_lti_cont")
    return fig1, fig2


def test_4blades_ltp_cont():
    """4 blades - LTP continuation analysis (ODI)."""
    print("\n[4 blades] LTP continuation analysis (ODI)...")
    problem = create_default_problem()
    problem.number_blades = 4
    problem.damper_activation = "ODI"
    problem.required_solver = "LTP"
    problem.continuation = "YES"

    rotor = RotorBuild(problem)
    analysis = ContinuationAnalysis(rotor).continuation()

    fig1, fig2 = plot_continuation(analysis.modal_solution, problem, "4 blades LTP Cont (ODI)")
    save_figures(fig1, fig2, "4blades_ltp_cont_odi")
    return fig1, fig2


def test_4blades_hd_cont():
    """4 blades - HD continuation analysis (ODI, 1 harmonic)."""
    print("\n[4 blades] HD continuation analysis (ODI, N=1)...")
    problem = create_default_problem()
    problem.number_blades = 4
    problem.damper_activation = "ODI"
    problem.required_solver = "HD"
    problem.hd_use_complex = True
    problem.continuation = "YES"
    problem.number_harmonics = 1

    rotor = RotorBuild(problem)
    analysis = ContinuationAnalysis(rotor).continuation()

    fig1, fig2 = plot_continuation(analysis.modal_solution, problem, "4 blades HD Cont (ODI, N=1)")
    save_figures(fig1, fig2, "4blades_hd_cont_odi")
    return fig1, fig2


# =============================================================================
# 5 BLADES TESTS
# =============================================================================

def test_5blades_ltp_all():
    """5 blades - LTP analysis (ALL)."""
    print("\n[5 blades] LTP analysis (ALL)...")
    problem = create_default_problem()
    problem.number_blades = 5
    problem.damper_activation = "ALL"
    problem.required_solver = "LTP"

    rotor = RotorBuild(problem)
    analysis = LTPStability(rotor).ltp_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "5 blades LTP (ALL)")
    save_figures(fig1, fig2, "5blades_ltp_all")
    return fig1, fig2


def test_5blades_ltp_odi():
    """5 blades - LTP analysis (ODI)."""
    print("\n[5 blades] LTP analysis (ODI)...")
    problem = create_default_problem()
    problem.number_blades = 5
    problem.damper_activation = "ODI"
    problem.required_solver = "LTP"

    rotor = RotorBuild(problem)
    analysis = LTPStability(rotor).ltp_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "5 blades LTP (ODI)")
    save_figures(fig1, fig2, "5blades_ltp_odi")
    return fig1, fig2


def test_5blades_hd_all():
    """5 blades - HD analysis (ALL, 2 harmonics)."""
    print("\n[5 blades] HD analysis (ALL, N=2)...")
    problem = create_default_problem()
    problem.number_blades = 5
    problem.damper_activation = "ALL"
    problem.required_solver = "HD"
    problem.number_harmonics = 2

    rotor = RotorBuild(problem)
    analysis = HDStability(rotor).hd_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "5 blades HD (ALL, N=2)")
    save_figures(fig1, fig2, "5blades_hd_all")
    return fig1, fig2


def test_5blades_hd_odi():
    """5 blades - HD analysis (ODI, 2 harmonics)."""
    print("\n[5 blades] HD analysis (ODI, N=2)...")
    problem = create_default_problem()
    problem.number_blades = 5
    problem.damper_activation = "ODI"
    problem.required_solver = "HD"
    problem.number_harmonics = 2

    rotor = RotorBuild(problem)
    analysis = HDStability(rotor).hd_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "5 blades HD (ODI, N=2)")
    save_figures(fig1, fig2, "5blades_hd_odi")
    return fig1, fig2


# =============================================================================
# 7 BLADES TESTS
# =============================================================================

def test_7blades_ltp_all():
    """7 blades - LTP analysis (ALL)."""
    print("\n[7 blades] LTP analysis (ALL)...")
    problem = create_default_problem()
    problem.number_blades = 7
    problem.damper_activation = "ALL"
    problem.required_solver = "LTP"

    rotor = RotorBuild(problem)
    analysis = LTPStability(rotor).ltp_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "7 blades LTP (ALL)")
    save_figures(fig1, fig2, "7blades_ltp_all")
    return fig1, fig2


def test_7blades_ltp_odi():
    """7 blades - LTP analysis (ODI)."""
    print("\n[7 blades] LTP analysis (ODI)...")
    problem = create_default_problem()
    problem.number_blades = 7
    problem.damper_activation = "ODI"
    problem.required_solver = "LTP"

    rotor = RotorBuild(problem)
    analysis = LTPStability(rotor).ltp_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "7 blades LTP (ODI)")
    save_figures(fig1, fig2, "7blades_ltp_odi")
    return fig1, fig2


def test_7blades_hd_all():
    """7 blades - HD analysis (ALL, 3 harmonics)."""
    print("\n[7 blades] HD analysis (ALL, N=3)...")
    problem = create_default_problem()
    problem.number_blades = 7
    problem.damper_activation = "ALL"
    problem.required_solver = "HD"
    problem.number_harmonics = 3

    rotor = RotorBuild(problem)
    analysis = HDStability(rotor).hd_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "7 blades HD (ALL, N=3)")
    save_figures(fig1, fig2, "7blades_hd_all")
    return fig1, fig2


def test_7blades_hd_odi():
    """7 blades - HD analysis (ODI, 3 harmonics)."""
    print("\n[7 blades] HD analysis (ODI, N=3)...")
    problem = create_default_problem()
    problem.number_blades = 7
    problem.damper_activation = "ODI"
    problem.required_solver = "HD"
    problem.number_harmonics = 3

    rotor = RotorBuild(problem)
    analysis = HDStability(rotor).hd_full_range()

    fig1, fig2 = plot_generic(analysis.modal_solution, problem, "7 blades HD (ODI, N=3)")
    save_figures(fig1, fig2, "7blades_hd_odi")
    return fig1, fig2


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Run all tests."""
    os.makedirs(FIGURES_DIR, exist_ok=True)

    print("=" * 70)
    print("STABILITY ANALYSIS TESTS")
    print(f"Output directory: {FIGURES_DIR}")
    print("=" * 70)

    # 4 blades tests (6 plots)
    print("\n" + "=" * 70)
    print("4 BLADES TESTS")
    print("=" * 70)
    test_4blades_lti()
    test_4blades_ltp()
    test_4blades_hd()
    test_4blades_lti_cont()
    test_4blades_ltp_cont()
    test_4blades_hd_cont()

    # 5 blades tests
    print("\n" + "=" * 70)
    print("5 BLADES TESTS")
    print("=" * 70)
    test_5blades_ltp_all()
    test_5blades_ltp_odi()
    test_5blades_hd_all()
    test_5blades_hd_odi()

    # 7 blades tests
    print("\n" + "=" * 70)
    print("7 BLADES TESTS")
    print("=" * 70)
    test_7blades_ltp_all()
    test_7blades_ltp_odi()
    test_7blades_hd_all()
    test_7blades_hd_odi()

    print("\n" + "=" * 70)
    print("ALL TESTS COMPLETE!")
    print("=" * 70)
    print(f"\nGenerated figures in '{FIGURES_DIR}/':")
    print("  4 blades: lti, ltp_odi, hd_odi, lti_cont, ltp_cont_odi, hd_cont_odi")
    print("  5 blades: ltp_all, ltp_odi, hd_all (N=2), hd_odi (N=2)")
    print("  7 blades: ltp_all, ltp_odi, hd_all (N=3), hd_odi (N=3)")

    plt.show()


if __name__ == "__main__":
    main()
