"""Main script for rotor stability analysis.

Converted from MATLAB CALL.m script.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from rotor_build import RotorBuild
from stability_analysis.base import StabilityAnalysis
from analysis.modal_participation import ModalParticipationAnalysis
from utils.plotting import MyPlot


def main():
    """Main analysis routine."""
    print("=" * 60)
    print("Rotor Dynamics Stability Analysis")
    print("=" * 60)
    
    # Build rotor system
    print("\n1. Building rotor system...")
    rotor = RotorBuild.build_all()
    print(f"   - Number of blades: {rotor.problem.number_blades}")
    print(f"   - Solver type: {rotor.problem.required_solver}")
    print(f"   - RPM range: {rotor.problem.lower_rotor_RPM} - {rotor.problem.higher_rotor_RPM}")
    
    # Run stability analysis
    print("\n2. Running stability analysis...")
    modes = StabilityAnalysis.run_stability(rotor)
    print(f"   - Number of operating points: {len(modes.modal_solution)}")
    if len(modes.modal_solution) > 0:
        print(f"   - Number of modes: {len(modes.modal_solution[0].damping)}")
    
    # Modal participation analysis
    if rotor.problem.required_solver in ["LTP", "HD"]:
        print("\n3. Computing modal participation...")
        modal_participation = ModalParticipationAnalysis(modes)
        print(f"   - Results computed for {len(modal_participation.modal_participation)} points")
    
    # Plotting
    print("\n4. Generating plots...")
    plotter = MyPlot()
    
    # Plot damping vs RPM for selected modes
    modes_to_plot = list(range(1, min(13, len(modes.modal_solution[0].damping) + 1)))
    print(f"   - Plotting modes: {modes_to_plot}")
    
    # fig = plotter.plot_damping_order(modes.modal_solution, modes=modes_to_plot)
    # Note: Plotting functionality needs to be fully implemented
    
    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
