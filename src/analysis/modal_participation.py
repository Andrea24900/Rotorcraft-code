"""Modal participation analysis.

Converted from MATLAB modal_participation_analysis.m classdef.
"""

import numpy as np
from typing import List
from dataclasses import dataclass, field


@dataclass
class ModalParticipation:
    """Container for modal participation at a single operating point."""
    
    OMEGA: float = 0.0
    OMEGA_RPM: float = 0.0
    phi_jkn: np.ndarray = field(default_factory=lambda: np.array([]))


class ModalParticipationAnalysis:
    """Modal participation analysis class.
    
    This class computes modal participation factors showing how each
    harmonic contributes to each mode.
    
    Attributes:
        modal_participation: List of modal participation results
        modal_solution: Modal solutions from stability analysis
        rotor_build: Rotor build object
    """
    
    def __init__(self, modes):
        """Initialize modal participation analysis.
        
        Args:
            modes: StabilityAnalysis object with computed modes
        """
        self.modal_solution = modes.modal_solution
        self.rotor_build = modes.rotor_build
        self.modal_participation: List[ModalParticipation] = []
        
        # Run analysis
        self._modal_part_all()
    
    def _modal_part_all(self):
        """Compute modal participation for all operating points."""
        for i in range(len(self.modal_solution)):
            if self.rotor_build.problem.required_solver == "LTP":
                phi_jkn_mat = self._mod_part_single_LTP(i)
            elif self.rotor_build.problem.required_solver == "HD":
                phi_jkn_mat = self._mod_part_single_HD(i)
            else:
                # LTI doesn't have modal participation analysis
                continue
            
            # Store results
            mod_part = ModalParticipation()
            mod_part.OMEGA = self.modal_solution[i].OMEGA
            mod_part.OMEGA_RPM = self.modal_solution[i].OMEGA_RPM
            mod_part.phi_jkn = phi_jkn_mat
            self.modal_participation.append(mod_part)
    
    def _mod_part_single_LTP(self, index: int) -> np.ndarray:
        """Compute modal participation for single LTP point.
        
        Args:
            index: Index in modal_solution
            
        Returns:
            Modal participation matrix phi_jkn
        """
        # TODO: Implement LTP modal participation
        # This requires:
        # 1. Computing time-varying eigenvectors V(t)
        # 2. FFT to get harmonic content
        # 3. Normalization
        
        # Placeholder
        return np.array([])
    
    def _mod_part_single_HD(self, index: int) -> np.ndarray:
        """Compute modal participation for single HD point.
        
        Args:
            index: Index in modal_solution
            
        Returns:
            Modal participation matrix phi_jkn
        """
        # TODO: Implement HD modal participation
        # This requires:
        # 1. Extracting harmonic components from expanded eigenvectors
        # 2. Computing complex coefficients
        # 3. Normalization
        
        # Placeholder
        return np.array([])


if __name__ == "__main__":
    # Test modal participation analysis
    from ..rotor_build import RotorBuild
    from ..stability_analysis.base import StabilityAnalysis
    
    print("Building rotor system...")
    rotor = RotorBuild.build_all()
    
    print("Running stability analysis...")
    modes = StabilityAnalysis.run_stability(rotor)
    
    print("Computing modal participation...")
    modal_part = ModalParticipationAnalysis(modes)
    
    print(f"Number of participation results: {len(modal_part.modal_participation)}")
