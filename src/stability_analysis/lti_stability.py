"""LTI stability analysis class.

Converted from MATLAB LTI_stability.m classdef.
"""

import numpy as np
from .base import StabilityAnalysis


class LTIStability(StabilityAnalysis):
    """Linear Time-Invariant stability analysis.
    
    This class provides methods for continuation analysis of classic LTI
    problems or analysis of a single OMEGA value.
    """
    
    def __init__(self, rotor_build):
        """Initialize LTI stability analysis.
        
        Args:
            rotor_build: RotorBuild object
        """
        super().__init__(rotor_build)
    
    def eigen_single_point(self, OMEGA: float) -> dict:
        """Compute eigenvalues and eigenvectors at single operating point.
        
        Args:
            OMEGA: Rotor speed (rad/s)
            
        Returns:
            Dictionary containing eigenvectors and eigenvalues
        """
        # Get state matrix at t=1 (arbitrary for LTI) and given OMEGA
        A = self.rotor_build.state_matrix_A_handles(1.0, OMEGA)
        
        # Compute eigenvalues and eigenvectors
        eigenvalues, eigenvectors = np.linalg.eig(A)
        
        eigensolution = {
            'eigenvectors': eigenvectors,
            'eigenvalues': eigenvalues
        }
        
        return eigensolution
    
    def eigen_full_range(self) -> 'LTIStability':
        """Compute eigenvalues over full RPM range.
        
        Returns:
            Self with populated modal_solution
        """
        # Assign OMEGA range
        self.assign_range_OMEGA()
        
        # Compute eigenvalues at each point
        for i in range(self.rotor_build.problem.number_points):
            eigensolution = self.eigen_single_point(self.modal_solution[i].OMEGA)
            
            # Extract damping (real part) and frequency (imaginary part)
            self.modal_solution[i].damping = np.real(eigensolution['eigenvalues'])
            self.modal_solution[i].frequency = np.imag(eigensolution['eigenvalues'])
            self.modal_solution[i].eigensolution = eigensolution
        
        return self


if __name__ == "__main__":
    # Test LTI stability analysis
    from ..rotor_build import RotorBuild
    
    print("Building rotor system...")
    rotor = RotorBuild.build_all()
    
    # Override to use LTI solver
    rotor.problem.required_solver = "LTI"
    
    print("Running LTI stability analysis...")
    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()
    
    print(f"Number of solutions: {len(lti.modal_solution)}")
    print(f"First point OMEGA: {lti.modal_solution[0].OMEGA_RPM:.2f} RPM")
    print(f"Number of eigenvalues: {len(lti.modal_solution[0].eigenvalues)}")
