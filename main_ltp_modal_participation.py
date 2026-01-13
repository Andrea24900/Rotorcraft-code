"""LTP Modal Participation Analysis with Interactive Plotting.

This script performs LTP (Linear Time-Periodic) continuation analysis and computes
modal participation factors, allowing you to interactively select which state and
mode to plot.

Modal Participation Analysis:
-----------------------------
Analyzes how each harmonic component contributes to each mode at each state.
For LTP methods, harmonics are extracted using Floquet theory with ODE integration
and FFT analysis.

Output:
-------
- Console: System info, stability results, available states/modes
- User Input: Interactive selection of state and mode to plot
- Plots: Modal participation factors (individual and summed harmonics)
"""

import numpy as np
import matplotlib.pyplot as plt
from src.rotor_build import RotorBuild
from src.stability_analysis.continuation_analysis import ContinuationAnalysis
from src.analysis.modal_participation import ModalParticipationAnalysis
from src.utils.plotting import MyPlot


def print_header(title):
    """Print a formatted header."""
    print("\n" + "="*70)
    print(f"  {title}")
    print("="*70)


def print_step(step_num, total_steps, description):
    """Print a formatted step indicator."""
    print(f"\n[{step_num}/{total_steps}] {description}")


def get_user_choice(prompt, min_val, max_val):
    """Get integer input from user within specified range.

    Args:
        prompt: Prompt message to display
        min_val: Minimum valid value (inclusive)
        max_val: Maximum valid value (inclusive)

    Returns:
        Selected integer value
    """
    while True:
        try:
            choice = int(input(prompt))
            if min_val <= choice <= max_val:
                return choice
            else:
                print(f"  [ERROR] Please enter a value between {min_val} and {max_val}")
        except ValueError:
            print("  [ERROR] Please enter a valid integer")
        except KeyboardInterrupt:
            print("\n\n  [CANCELLED] Exiting...")
            exit(0)


def main():
    """Run LTP modal participation analysis with interactive plotting."""

    print_header("LTP Modal Participation Analysis")

    # ========================================================================
    # Step 1: Build rotor system
    # ========================================================================
    print_step(1, 6, "Building rotor system...")

    rotor = RotorBuild.build_all()

    # Configure for LTP with continuation
    rotor.problem.required_solver = "LTP"
    rotor.problem.continuation = "YES"
    rotor.problem.number_harmonics = 2  # Using 2 harmonics for LTP analysis

    print(f"  [OK] Rotor system configured")
    print(f"    - Number of blades: {rotor.problem.number_blades}")
    print(f"    - Solver type: {rotor.problem.required_solver}")
    print(f"    - Number of harmonics (N): {rotor.problem.number_harmonics}")
    print(f"    - RPM range: {rotor.problem.lower_rotor_RPM} - {rotor.problem.higher_rotor_RPM}")
    print(f"    - Number of points: {rotor.problem.number_points}")

    # Calculate dimensions (LTP has original state space size)
    state_size = 2 * rotor.problem.number_blades + 4
    print(f"    - State space dimension: {state_size}")

    # ========================================================================
    # Step 2: Run LTP continuation analysis
    # ========================================================================
    print_step(2, 6, "Running LTP continuation analysis...")

    cont_analysis = ContinuationAnalysis(rotor)
    cont_analysis = cont_analysis.continuation()

    print(f"  [OK] LTP continuation completed")
    print(f"    - Solution points: {len(cont_analysis.modal_solution)}")

    # Check stability
    max_damping = max([np.max(sol.damping) for sol in cont_analysis.modal_solution])
    if max_damping < 0:
        print(f"    - System is STABLE (max damping: {max_damping:.4f})")
    else:
        print(f"    - System is UNSTABLE (max damping: {max_damping:.4f})")

    # ========================================================================
    # Step 3: Compute modal participation factors
    # ========================================================================
    print_step(3, 6, "Computing modal participation factors...")
    print("  [INFO] LTP method uses Floquet theory with ODE integration and FFT")
    print("  [INFO] This may take longer than HD method...")
    print("\n  [WARNING] LTP modal participation is dependent on the integer")
    print("  [WARNING] multiple of Omega added to the complex logarithm.")
    print("  [WARNING] Caution is advised when interpreting results.\n")

    modal_part_analysis = ModalParticipationAnalysis(cont_analysis)

    print(f"  [OK] Modal participation computed")
    print(f"    - Participation results: {len(modal_part_analysis.modal_participation)}")

    # Get dimensions
    first_result = modal_part_analysis.modal_participation[0]
    n_states, n_modes, n_harmonics = first_result.phi_jkn.shape

    print(f"    - Number of states (j): {n_states}")
    print(f"    - Number of modes (k): {n_modes}")
    print(f"    - Number of harmonics (n): {n_harmonics}")

    # Calculate harmonic range
    n_harm_side = n_harmonics // 2
    print(f"    - Harmonic range: n = {-n_harm_side} to {+n_harm_side}")

    # ========================================================================
    # Step 4: Display overview of all LTP modes
    # ========================================================================
    print_step(4, 6, "Displaying overview of all LTP modes...")

    plotter = MyPlot()
    rpm_min = rotor.problem.lower_rotor_RPM
    rpm_max = rotor.problem.higher_rotor_RPM

    # Calculate damping limits for overview plot
    all_damping = np.concatenate([sol.damping for sol in cont_analysis.modal_solution])
    damping_min = np.min(all_damping) - 0.5
    damping_max = max(1.0, np.max(all_damping) + 0.5)

    # Get total number of eigenvalues
    n_eigenvalues = len(cont_analysis.modal_solution[0].damping)

    # Create overview damping plot for all modes
    print("  - Plotting damping of all modes...")
    fig_overview = plotter.plot_damping_order(
        cont_analysis.modal_solution,
        modes=list(range(1, n_eigenvalues + 1)),  # All modes (1-based)
        xlimits=(rpm_min, rpm_max),
        ylimits=(damping_min, damping_max)
    )
    plt.figure(fig_overview.number)
    plt.suptitle(f'LTP Continuation: Damping of All Modes vs RPM (N={rotor.problem.number_harmonics} harmonics)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    # Display the overview plot non-blocking so user can review it
    plt.show(block=False)
    plt.pause(0.1)  # Brief pause to ensure plot is rendered

    print(f"  [OK] Overview plot displayed")
    print(f"    - Review the damping plot to help select mode of interest")
    print(f"    - Plot shows all {n_eigenvalues} eigenvalues")

    # ========================================================================
    # Step 5: Interactive selection of state and mode
    # ========================================================================
    print_step(5, 6, "Select state and mode to plot...")

    print("\n  Available states (1-indexed):")
    print(f"    1 to {n_states}")
    print("\n  State mapping:")
    print(f"    1-4: MBC velocity states (Z_dot)")
    print(f"    5-6: Hub velocity states (x_H,y_H)")
    print(f"    7-10: MBC states (Z)")
    print(f"    11-12: Hub states (x_H,y_H)")

    # Calculate number of mode pairs (complex conjugate pairs)
    n_mode_pairs = n_modes // 2

    print("\n  Available mode pairs (1-indexed):")
    print(f"    1 to {n_mode_pairs}")
    print(f"    (Each mode pair is a complex conjugate eigenvalue pair)")
    print(f"    (Mode pair k corresponds to eigenvalues {2*1-1},{2*1} ... {2*n_mode_pairs-1},{2*n_mode_pairs})")

    # Get user input (1-based, convert to 0-based internally)
    print("\n" + "-"*70)
    state_num = get_user_choice(
        f"\n  Select state number [1-{n_states}]: ",
        1, n_states
    )

    mode_pair_num = get_user_choice(
        f"  Select mode pair number [1-{n_mode_pairs}]: ",
        1, n_mode_pairs
    )

    # Convert to 0-based indexing for internal use
    state_idx = state_num - 1
    mode_idx = mode_pair_num - 1  # 0-based mode pair index

    print(f"\n  [OK] Selected: State {state_num}, Mode Pair {mode_pair_num}")

    # ========================================================================
    # Step 6: Generate plots
    # ========================================================================
    print_step(6, 6, "Generating modal participation plots...")

    # Plot 1: Individual harmonics
    print("  - Plotting individual harmonics...")
    fig1 = plotter.plot_mod_part(
        modal_participation=modal_part_analysis.modal_participation,
        state_index=state_idx,
        mode_index=mode_idx,
        xlimits=(rpm_min, rpm_max),
        ylimits=(0, 1),
        show_legend=True
    )

    plt.figure(fig1.number)
    plt.suptitle(f'Modal Participation: State {state_num}, Mode Pair {mode_pair_num} (Individual Harmonics)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    # Plot 2: Summed harmonics (±n combined)
    print("  - Plotting summed harmonics...")
    fig2 = plotter.plot_mod_part_sum(
        modal_participation=modal_part_analysis.modal_participation,
        state_index=state_idx,
        mode_index=mode_idx,
        xlimits=(rpm_min, rpm_max),
        ylimits=(0, 1),
        show_legend=True
    )

    plt.figure(fig2.number)
    plt.suptitle(f'Modal Participation: State {state_num}, Mode Pair {mode_pair_num} (Summed Harmonics)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    # Plot 3: LTP Mode Damping
    print("  - Plotting LTP mode damping...")
    fig3 = plotter.plot_damping_order(
        cont_analysis.modal_solution,
        modes=[mode_pair_num * 2 - 1, mode_pair_num * 2],  # Both conjugate pair eigenvalues
        xlimits=(rpm_min, rpm_max),
        ylimits=None  # Auto-scale
    )
    plt.figure(fig3.number)
    plt.suptitle(f'LTP Mode Damping: Mode Pair {mode_pair_num}',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    print(f"  [OK] Plots generated")

    # ========================================================================
    # Summary
    # ========================================================================
    print_header("Analysis Complete!")

    print("\n  Summary:")
    print(f"    - State analyzed: {state_num}")
    print(f"    - Mode pair analyzed: {mode_pair_num}")
    print(f"    - Eigenvalue indices: {mode_pair_num * 2 - 1}, {mode_pair_num * 2}")
    print(f"    - Harmonic range: n = {-n_harm_side} to {+n_harm_side}")
    print(f"    - RPM range: {rpm_min} - {rpm_max}")

    print("\n  Generated Plots:")
    print("    1. Individual Harmonics - Shows each n separately")
    print("    2. Summed Harmonics - Combines ±n contributions")
    print("    3. Mode Damping - Shows stability of the selected mode pair")

    print("\n  Interpretation Guide:")
    print("    - Each harmonic (n) shows its contribution to the selected mode")
    print("    - Participation factors sum to 1.0 across all harmonics")
    print("    - Line styles:")
    print("      • Solid lines: Negative harmonics (n < 0) or ±n sums")
    print("      • Dash-dot line: Zero harmonic (n = 0)")
    print("      • Dashed lines: Positive harmonics (n > 0)")
    print("    - Damping plot shows eigenvalue real parts (negative = stable)")

    print("\n  LTP Method Notes:")
    print("    - Uses Floquet theory with ODE integration over one period")
    print("    - Applies FFT to extract harmonic content from time-varying eigenvectors")
    print("    - More accurate for strongly time-periodic systems")
    print("    - Computationally more expensive than HD method")

    print("\n  Displaying plots...")
    print("  (Close plot windows to exit)")
    print("\n" + "="*70 + "\n")

    plt.show()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n" + "="*70)
        print("  [CANCELLED] Analysis interrupted by user")
        print("="*70 + "\n")
    except Exception as e:
        print("\n\n" + "="*70)
        print(f"  [ERROR] An error occurred: {str(e)}")
        print("="*70 + "\n")
        raise
