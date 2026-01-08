"""Rotor build class for constructing state-space matrices.

Converted from MATLAB rotor_build.m classdef.
"""

import numpy as np
from typing import Callable
from .config.problem_definition import ProblemDefinition, create_default_problem
from .utils.matrix_repositories import build_mass_matrix, build_damping_matrix


class RotorBuild:
    """Main rotor class containing methods to build state-space matrices.
    
    This class contains all methods to build the state-space matrix
    for given rotor characteristics.
    
    Attributes:
        problem: Problem definition with rotor configuration
        state_matrix_A_handles: Function handle for state matrix A(t, OMEGA)
    """
    
    def __init__(self, problem: ProblemDefinition):
        """Initialize rotor build with problem definition.
        
        Args:
            problem: Problem definition containing rotor characteristics
        """
        self.problem = problem
        self.state_matrix_A_handles = None
    
    @staticmethod
    def build_all() -> 'RotorBuild':
        """Build rotor from default problem definition.
        
        This is the main constructor method that:
        1. Loads problem definition
        2. Creates rotor object
        3. Adds state matrix A to rotor
        
        Returns:
            Fully constructed RotorBuild object
        """
        # Load problem definition
        problem = create_default_problem()
        
        # Create rotor object
        obj = RotorBuild(problem)
        
        # Add state matrix A
        obj = obj.add_A_to_rotor()
        
        return obj
    
    def build_mass_matrix(self) -> np.ndarray:
        """Build mass matrix from rotor characteristics.
        
        Returns:
            Mass matrix M
        """
        return build_mass_matrix(self.problem)
    
    def build_damping_stiffness_matrices(
        self
    ) -> tuple[Callable[[float, float], np.ndarray], Callable[[float, float], np.ndarray]]:
        """Build damping and stiffness matrices as function handles.
        
        Returns:
            Tuple of (damping_matrix_C, stiffness_matrix_K) function handles
        """
        return build_damping_matrix(self.problem)
    
    def build_state_matrix(self) -> Callable[[float, float], np.ndarray]:
        """Build state-space matrix A as a function handle.
        
        The state matrix has the form:
        A = [[-M^-1 * C,  -M^-1 * K],
             [    I,           0    ]]
        
        Returns:
            Function handle A(t, rotor_OMEGA) that returns state matrix
        """
        # Get mass matrix (constant)
        M = build_mass_matrix(self.problem)
        M_inv = np.linalg.inv(M)
        
        # Get damping and stiffness function handles
        C_func, K_func = build_damping_matrix(self.problem)
        
        # Get size
        n = M.shape[0]
        
        def state_matrix_A(t: float, rotor_OMEGA: float) -> np.ndarray:
            """State matrix as function of time and rotor speed.
            
            Args:
                t: Time
                rotor_OMEGA: Rotor angular velocity (rad/s)
                
            Returns:
                State matrix A (2n x 2n)
            """
            C = C_func(t, rotor_OMEGA)
            K = K_func(t, rotor_OMEGA)
            
            # Upper row: [-M^-1 * C, -M^-1 * K]
            upper_left = -M_inv @ C
            upper_right = -M_inv @ K
            
            # Lower row: [I, 0]
            lower_left = np.eye(n)
            lower_right = np.zeros((n, n))
            
            # Combine into full state matrix
            A = np.block([
                [upper_left, upper_right],
                [lower_left, lower_right]
            ])
            
            return A
        
        return state_matrix_A
    
    def add_A_to_rotor(self) -> 'RotorBuild':
        """Add state matrix A to rotor object.
        
        Returns:
            Self (for method chaining)
        """
        state_matrix_A = self.build_state_matrix()
        self.state_matrix_A_handles = state_matrix_A
        return self


if __name__ == "__main__":
    # Test rotor build
    print("Building rotor system...")
    rotor = RotorBuild.build_all()
    
    print(f"Problem configuration:")
    print(f"  Number of blades: {rotor.problem.number_blades}")
    print(f"  Solver: {rotor.problem.required_solver}")
    print(f"  Damper connection: {rotor.problem.damper_connection}")
    
    # Test state matrix evaluation
    if rotor.state_matrix_A_handles is not None:
        A = rotor.state_matrix_A_handles(0.0, 25.0)
        print(f"\nState matrix A shape: {A.shape}")
        print(f"State matrix A condition number: {np.linalg.cond(A):.2e}")
