"""HD (Harmonic Decomposition) stability analysis using Lopez and Prasad method.

Converted from MATLAB HD_stability.m classdef.
"""

import numpy as np
from .base import StabilityAnalysis
from .ltp_stability import LTPStability


class HDStability(LTPStability):
    """Harmonic Domain stability analysis using Harmonic Decomposition method.
    
    This class extends LTP stability to use Harmonic Decomposition method
    for analyzing periodic systems in the frequency domain.
    """
    
    def __init__(self, rotor_build):
        """Initialize HD stability analysis.
        
        Args:
            rotor_build: RotorBuild object
        """
        super().__init__(rotor_build)
    
    def HD_single_point(self, OMEGA: float, T: float) -> dict:
        """Compute eigenvalues using Harmonic Decomposition method.

        Args:
            OMEGA: Rotor speed (rad/s)
            T: Period (s)

        Returns:
            Dictionary containing eigenvectors and eigenvalues
        """
        # Create time-dependent state matrix handle
        A = lambda t: self.rotor_build.state_matrix_A_handles(t, OMEGA)

        # Time vector for sampling using problem definition
        time = np.linspace(0, T, self.rotor_build.problem.time_samples)

        # Compute Harmonic Decomposition matrix
        A_HD = StabilityAnalysis.HD_computer(
            A,
            time,
            self.rotor_build.problem.number_harmonics,
            OMEGA
        )

        # Compute eigenvalues and eigenvectors of expanded system
        eigenvalues, eigenvectors = np.linalg.eig(A_HD)

        eigensolution = {
            'eigevectors': eigenvectors,
            'eigenvalues': eigenvalues
        }

        return eigensolution
    
    def HD_full_range(self) -> 'HDStability':
        """Compute eigenvalues over full RPM range using HD method.
        
        Returns:
            Self with populated modal_solution
        """
        # Assign OMEGA range and periods
        self.assign_range_OMEGA()
        self.assign_period_T()
        
        # Compute eigenvalues at each point
        for i in range(self.rotor_build.problem.number_points):
            eigensolution = self.HD_single_point(
                self.modal_solution[i].OMEGA,
                self.modal_solution[i].T
            )
            
            # Store results
            self.modal_solution[i].damping = np.real(eigensolution['eigenvalues'])
            self.modal_solution[i].frequency = np.imag(eigensolution['eigenvalues'])
            self.modal_solution[i].eigensolution = eigensolution
        
        return self


if __name__ == "__main__":
    # Test HD stability analysis
    from ..rotor_build import RotorBuild
    
    print("Building rotor system...")
    rotor = RotorBuild.build_all()
    
    # Override to use HD solver
    rotor.problem.required_solver = "HD"
    rotor.problem.number_harmonics = 1
    
    print("Running HD stability analysis...")
    hd = HDStability(rotor)
    hd = hd.HD_full_range()
    
    print(f"Number of solutions: {len(hd.modal_solution)}")
    print(f"First point OMEGA: {hd.modal_solution[0].OMEGA_RPM:.2f} RPM")
    if len(hd.modal_solution) > 0:
        print(f"Number of eigenvalues: {len(hd.modal_solution[0].eigenvalues)}")
