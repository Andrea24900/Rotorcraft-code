"""LTP stability analysis class using Floquet theory.

Converted from MATLAB LTP_stability.m classdef.
"""

import numpy as np
from .base import StabilityAnalysis


class LTPStability(StabilityAnalysis):
    """Linear Time-Periodic stability analysis using Floquet theory.
    
    This class uses monodromy matrix computation to determine
    characteristic exponents of time-periodic systems.
    """
    
    def __init__(self, rotor_build):
        """Initialize LTP stability analysis.
        
        Args:
            rotor_build: RotorBuild object
        """
        super().__init__(rotor_build)
    
    def LTP_single_point(self, OMEGA: float, T: float) -> dict:
        """Compute characteristic multipliers at single operating point.
        
        Args:
            OMEGA: Rotor speed (rad/s)
            T: Period (s)
            
        Returns:
            Dictionary containing characteristic multipliers and exponents
        """
        # Create time-dependent state matrix handle
        A = lambda t: self.rotor_build.state_matrix_A_handles(t, OMEGA)
        
        # Compute monodromy matrix
        monodromy_matrix_M = self.monodromy_computer(A, T)
        
        # Compute characteristic multipliers (eigenvalues of M)
        char_multipliers = np.linalg.eigvals(monodromy_matrix_M)
        
        # Compute characteristic exponents
        # Real part: σ = (1/2T) * ln(|μ|²)
        real_char_exp = (1 / (2 * T)) * np.log(
            np.real(char_multipliers)**2 + np.imag(char_multipliers)**2
        )
        
        # Imaginary part: ω = (1/T) * atan2(Im(μ), Re(μ))
        imag_char_exp = (1 / T) * np.arctan2(
            np.imag(char_multipliers),
            np.real(char_multipliers)
        )
        
        char_solutions = {
            'char_multipliers': char_multipliers,
            'real_char_exp': real_char_exp,
            'imag_char_exp': imag_char_exp
        }
        
        return char_solutions
    
    def LTP_full_range(self) -> 'LTPStability':
        """Compute characteristic exponents over full RPM range.
        
        Returns:
            Self with populated modal_solution
        """
        # Assign OMEGA range and periods
        self.assign_range_OMEGA()
        self.assign_period_T()
        
        # Compute characteristic solutions at each point
        for i in range(self.rotor_build.problem.number_points):
            char_solution = self.LTP_single_point(
                self.modal_solution[i].OMEGA,
                self.modal_solution[i].T
            )
            
            # Store results
            self.modal_solution[i].char_solution = char_solution
            self.modal_solution[i].damping = char_solution['real_char_exp']
            self.modal_solution[i].frequency = char_solution['imag_char_exp']
        
        return self


if __name__ == "__main__":
    # Test LTP stability analysis
    from ..rotor_build import RotorBuild
    
    print("Building rotor system...")
    rotor = RotorBuild.build_all()
    
    # Override to use LTP solver
    rotor.problem.required_solver = "LTP"
    
    print("Running LTP stability analysis...")
    ltp = LTPStability(rotor)
    ltp = ltp.LTP_full_range()
    
    print(f"Number of solutions: {len(ltp.modal_solution)}")
    print(f"First point OMEGA: {ltp.modal_solution[0].OMEGA_RPM:.2f} RPM")
    if len(ltp.modal_solution) > 0:
        print(f"Number of characteristic exponents: {len(ltp.modal_solution[0].damping)}")
