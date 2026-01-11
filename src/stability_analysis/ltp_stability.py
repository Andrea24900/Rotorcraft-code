"""LTP stability analysis class using Floquet theory.

This module implements Linear Time-Periodic (LTP) stability analysis for rotor
systems with time-varying coefficients. It uses Floquet theory to determine
system stability through monodromy matrix computation.

Theoretical Background:
----------------------
For a linear time-periodic system: dx/dt = A(t)x, where A(t+T) = A(t)

Floquet theory states that solutions can be expressed as:
    x(t) = Φ(t)x₀, where Φ(t) is the fundamental solution matrix

The monodromy matrix M is defined as: M = Φ(T)
Its eigenvalues (characteristic multipliers μ) determine stability:
    - |μ| < 1: Stable
    - |μ| > 1: Unstable
    - |μ| = 1: Marginally stable

Characteristic exponents λ are related to multipliers by: μ = e^(λT)
    - Real part (σ): σ = (1/2T) * ln(|μ|²) → damping
    - Imaginary part (ω): ω = (1/T) * atan2(Im(μ), Re(μ)) → frequency

Converted from MATLAB LTP_stability.m classdef.
"""

import numpy as np
from .base import StabilityAnalysis


class LTPStability(StabilityAnalysis):
    """Linear Time-Periodic stability analysis using Floquet theory.

    This class analyzes the stability of rotor systems with time-periodic
    coefficients (e.g., due to blade passing effects, asymmetric rotor
    properties, or periodic excitation). It computes the monodromy matrix
    by integrating the fundamental solution matrix over one period and
    extracts characteristic exponents to determine system stability.

    Key Features:
    ------------
    - Monodromy matrix computation via numerical integration
    - Characteristic multiplier and exponent extraction
    - Support for multiple operating points (RPM sweep)
    - Compatible with continuation methods

    Attributes:
    ----------
    rotor_build : RotorBuild
        Rotor system configuration and state matrices
    modal_solution : list[ModalSolution]
        List of solutions at each operating point, containing:
        - OMEGA: Rotor speed (rad/s)
        - OMEGA_RPM: Rotor speed (RPM)
        - T: Period (s)
        - char_solution: Dict with multipliers and exponents
        - damping: Real part of characteristic exponents (σ)
        - frequency: Imaginary part of characteristic exponents (ω)

    Methods:
    -------
    ltp_single_point(OMEGA, T)
        Compute characteristic exponents at a single operating point
    ltp_full_range()
        Compute characteristic exponents over full RPM range

    Example:
    -------
    >>> rotor = RotorBuild.build_all()
    >>> rotor.problem.required_solver = "LTP"
    >>> ltp = LTPStability(rotor)
    >>> ltp = ltp.LTP_full_range()
    >>> # Check stability
    >>> max_damping = max([np.max(sol.damping) for sol in ltp.modal_solution])
    >>> is_stable = max_damping < 0
    """
    
    def __init__(self, rotor_build):
        """Initialize LTP stability analysis.
        
        Args:
            rotor_build: RotorBuild object
        """
        super().__init__(rotor_build)
    
    def ltp_single_point(self, OMEGA: float, T: float) -> dict:
        """Compute characteristic multipliers and exponents at single operating point.

        This method computes the monodromy matrix M = Φ(T) by integrating each
        column of the fundamental solution matrix over one period. The eigenvalues
        of M (characteristic multipliers) are then converted to characteristic
        exponents which represent the system's damping and frequency.

        Algorithm:
        ---------
        1. Create time-dependent state matrix handle A(t, OMEGA)
        2. Compute monodromy matrix M by integrating: dΦ/dt = A(t)Φ, Φ(0) = I
        3. Extract characteristic multipliers: μ = eig(M)
        4. Convert to characteristic exponents:
           - Real part (damping): σ = (1/2T) * ln(|μ|²)
           - Imaginary part (frequency): ω = (1/T) * atan2(Im(μ), Re(μ))

        Args:
        ----
        OMEGA : float
            Rotor speed (rad/s)
        T : float
            Period of oscillation (s), typically T = 2π/OMEGA

        Returns:
        -------
        dict
            Dictionary containing:
            - 'char_multipliers': Complex array of characteristic multipliers (μ)
            - 'real_char_exp': Real part of characteristic exponents (σ) - damping
            - 'imag_char_exp': Imaginary part of characteristic exponents (ω) - frequency

        Notes:
        -----
        - Negative real exponents (σ < 0) indicate stable modes
        - Positive real exponents (σ > 0) indicate unstable modes
        - Imaginary exponents represent oscillation frequencies
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
    
    def ltp_full_range(self) -> 'LTPStability':
        """Compute characteristic exponents over full RPM range.

        This method performs a sweep of operating points across the specified
        RPM range, computing characteristic exponents at each point using
        Floquet theory. The results are stored in the modal_solution attribute.

        Process:
        -------
        1. Generate array of OMEGA values from problem definition
        2. Compute corresponding periods T = 2π/OMEGA
        3. For each operating point:
           - Compute monodromy matrix
           - Extract characteristic multipliers and exponents
           - Store damping and frequency values

        Returns:
        -------
        LTPStability
            Self with populated modal_solution list containing results at
            each operating point

        Attributes Modified:
        -------------------
        modal_solution : list[ModalSolution]
            Populated with characteristic exponents at each RPM point.
            Each entry contains:
            - OMEGA, OMEGA_RPM: Operating speed
            - T: Period
            - char_solution: Full solution dictionary
            - damping: Real characteristic exponents (stability indicators)
            - frequency: Imaginary characteristic exponents (oscillation frequencies)

        Notes:
        -----
        - This method can be computationally expensive for many points
        - Each point requires integration of n columns of the fundamental matrix
        - For n DOF system, expect n complex conjugate pairs of exponents
        - Use continuation methods for more efficient tracking of specific modes
        """
        # Assign OMEGA range and periods
        self.assign_range_OMEGA()
        self.assign_period_T()
        
        # Compute characteristic solutions at each point
        for i in range(self.rotor_build.problem.number_points):
            char_solution = self.ltp_single_point(
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
    ltp = ltp.ltp_full_range()
    
    print(f"Number of solutions: {len(ltp.modal_solution)}")
    print(f"First point OMEGA: {ltp.modal_solution[0].OMEGA_RPM:.2f} RPM")
    if len(ltp.modal_solution) > 0:
        print(f"Number of characteristic exponents: {len(ltp.modal_solution[0].damping)}")
