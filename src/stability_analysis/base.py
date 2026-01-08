"""Base stability analysis class.

Converted from MATLAB stability_analysis.m classdef.
"""

import numpy as np
from typing import Callable, List, Dict, Any
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field


@dataclass
class ModalSolution:
    """Container for modal solution at a single operating point."""
    
    OMEGA: float = 0.0  # Rotor speed (rad/s)
    OMEGA_RPM: float = 0.0  # Rotor speed (RPM)
    T: float = 0.0  # Period (s)
    damping: np.ndarray = field(default_factory=lambda: np.array([]))
    frequency: np.ndarray = field(default_factory=lambda: np.array([]))
    eigensolution: Dict[str, Any] = field(default_factory=dict)
    char_solution: Dict[str, Any] = field(default_factory=dict)


class StabilityAnalysis:
    """Base class for stability analysis.
    
    This class contains main functions for various stability analyses
    including LTI, LTP, and HD methods.
    
    Attributes:
        rotor_build: Rotor build object with state matrices
        modal_solution: List of modal solutions at different operating points
    """
    
    def __init__(self, rotor_build):
        """Initialize stability analysis.
        
        Args:
            rotor_build: RotorBuild object
        """
        self.rotor_build = rotor_build
        self.modal_solution: List[ModalSolution] = []
    
    @staticmethod
    def run_stability(rotor_build):
        """Run full stability analysis according to problem definition.
        
        This is the main entry point that selects the appropriate analysis
        method based on the problem configuration.
        
        Args:
            rotor_build: RotorBuild object
            
        Returns:
            StabilityAnalysis subclass instance with computed results
        """
        from .lti_stability import LTIStability
        from .ltp_stability import LTPStability
        from .hd_stability import HDStability
        
        solver_type = rotor_build.problem.required_solver
        use_continuation = rotor_build.problem.continuation == "YES"
        
        if solver_type == "LTI":
            if use_continuation:
                # TODO: Implement continuation analysis
                raise NotImplementedError("Continuation analysis not yet implemented")
            else:
                obj = LTIStability(rotor_build)
                obj = obj.eigen_full_range()
                
        elif solver_type == "LTP":
            if use_continuation:
                raise NotImplementedError("Continuation analysis not yet implemented")
            else:
                obj = LTPStability(rotor_build)
                obj = obj.LTP_full_range()
                
        elif solver_type == "HD":
            if use_continuation:
                raise NotImplementedError("Continuation analysis not yet implemented")
            else:
                obj = HDStability(rotor_build)
                obj = obj.HD_full_range()
        else:
            raise ValueError(f"Unknown solver type: {solver_type}")
        
        return obj
    
    @staticmethod
    def monodromy_computer(
        A_handle_time: Callable[[float], np.ndarray],
        period_T: float
    ) -> np.ndarray:
        """Compute monodromy matrix for LTP system.
        
        The monodromy matrix M is computed by integrating the fundamental
        solution matrix Φ(t) over one period: M = Φ(T).
        
        Args:
            A_handle_time: Function handle A(t) returning state matrix
            period_T: Period of oscillation
            
        Returns:
            Monodromy matrix M
        """
        n = A_handle_time(0).shape[0]
        phi_0_matrix = np.eye(n)
        monodromy_matrix_M = np.zeros((n, n))
        
        def odefun_monodromy(t, phi_vec, A_func):
            """ODE function for fundamental solution matrix."""
            phi = phi_vec.reshape((n, n))
            dphi_dt = A_func(t) @ phi
            return dphi_dt.flatten()
        
        # Integrate each column of the fundamental solution matrix
        for j in range(n):
            phi_0 = phi_0_matrix[:, j]
            
            sol = solve_ivp(
                lambda t, phi: odefun_monodromy(t, phi, A_handle_time),
                [0, period_T],
                phi_0.flatten() if len(phi_0.shape) == 1 else phi_0,
                method='DOP853',
                rtol=1e-6,
                atol=1e-8
            )
            
            # Extract final value
            phi_T = sol.y[:, -1].reshape((n, 1))
            monodromy_matrix_M[:, j] = phi_T.flatten()
        
        return monodromy_matrix_M
    
    @staticmethod
    def HD_computer(
        A_handle_time: Callable[[float], np.ndarray],
        time: np.ndarray,
        number_harmonics: int,
        OMEGA: float
    ) -> np.ndarray:
        """Compute Hill-determinant matrix for HD analysis.
        
        This method constructs the expanded Hill matrix by computing
        Fourier coefficients of the time-periodic state matrix A(t).
        
        Args:
            A_handle_time: Function handle A(t)
            time: Time vector for sampling
            number_harmonics: Number of harmonics to include
            OMEGA: Fundamental frequency
            
        Returns:
            Hill-determinant matrix A_HD
        """
        state_number = A_handle_time(0).shape[0]
        size_A_HD = (1 + 2 * number_harmonics) * state_number
        A_HD = np.zeros((size_A_HD, size_A_HD))
        
        # Compute time realizations of A(t)
        time_realisation_A = np.zeros((state_number, state_number, len(time)))
        for i, t in enumerate(time):
            time_realisation_A[:, :, i] = A_handle_time(t)
        
        # Compute FFT coefficients for each element
        A_coeff = np.zeros((state_number, state_number, 2 * number_harmonics + 1))
        
        for j in range(state_number):
            for k in range(state_number):
                # FFT of A[j,k](t)
                full_transform = np.fft.fft(time_realisation_A[j, k, :]) / len(time)
                
                # Extract coefficients: A_0, A_1c, A_1s, A_2c, A_2s, ...
                A_coeff[j, k, 0] = full_transform[0].real  # A_0
                
                index_coeff = 1
                for n in range(1, number_harmonics + 1):
                    A_coeff[j, k, index_coeff] = 2 * full_transform[n].real  # A_nc
                    A_coeff[j, k, index_coeff + 1] = -2 * full_transform[n].imag  # A_ns
                    index_coeff += 2
        
        # TODO: Implement H-matrix assignment functions
        # This is a placeholder for the complete HD matrix construction
        # The full implementation requires the H assignment helper functions
        
        return A_HD
    
    def assign_range_OMEGA(self) -> 'StabilityAnalysis':
        """Assign range of rotor speeds to analyze.
        
        Returns:
            Self (for method chaining)
        """
        rotor_OMEGA_vector = np.linspace(
            self.rotor_build.problem.lower_rotor_RPM / 60 * 2 * np.pi,
            self.rotor_build.problem.higher_rotor_RPM / 60 * 2 * np.pi,
            self.rotor_build.problem.number_points
        )
        
        rotor_OMEGA_vector_RPM = np.linspace(
            self.rotor_build.problem.lower_rotor_RPM,
            self.rotor_build.problem.higher_rotor_RPM,
            self.rotor_build.problem.number_points
        )
        
        # Reverse if solving right to left
        if self.rotor_build.problem.solution_direction == "R2L":
            rotor_OMEGA_vector = rotor_OMEGA_vector[::-1]
            rotor_OMEGA_vector_RPM = rotor_OMEGA_vector_RPM[::-1]
        
        # Initialize modal solution list
        self.modal_solution = []
        for i in range(self.rotor_build.problem.number_points):
            modal_sol = ModalSolution()
            modal_sol.OMEGA = rotor_OMEGA_vector[i]
            modal_sol.OMEGA_RPM = rotor_OMEGA_vector_RPM[i]
            self.modal_solution.append(modal_sol)
        
        return self
    
    def assign_period_T(self) -> 'StabilityAnalysis':
        """Assign period T = 2π/Ω for each operating point.
        
        Returns:
            Self (for method chaining)
        """
        for modal_sol in self.modal_solution:
            modal_sol.T = 2 * np.pi / modal_sol.OMEGA
        
        return self


if __name__ == "__main__":
    # Test stability analysis base class
    from ..rotor_build import RotorBuild
    
    print("Building rotor system...")
    rotor = RotorBuild.build_all()
    
    print("Initializing stability analysis...")
    stability = StabilityAnalysis(rotor)
    
    print("Assigning OMEGA range...")
    stability = stability.assign_range_OMEGA()
    print(f"Number of points: {len(stability.modal_solution)}")
    print(f"First OMEGA: {stability.modal_solution[0].OMEGA:.2f} rad/s")
    print(f"Last OMEGA: {stability.modal_solution[-1].OMEGA:.2f} rad/s")
