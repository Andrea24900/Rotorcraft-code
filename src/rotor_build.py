"""Rotor build class for constructing state-space matrices.

Converted from MATLAB rotor_build.m classdef.
Builds state-space representation for rotorcraft dynamics analysis.
"""

import numpy as np
from typing import Callable
from .config.problem_definition import ProblemDefinition, create_default_problem
from .utils.matrix_generation_GR import build_mass_matrix, build_damping_matrix


class RotorBuild:
    """Main rotor class for building state-space matrices.

    This class constructs the state-space representation of a rotorcraft system:
    - Mass matrix M(t, Ω) (time and speed dependent, via MBC transformation)
    - Damping matrix C(t, Ω) (time and speed dependent)
    - Stiffness matrix K(t, Ω) (time and speed dependent)
    - State matrix A(t, Ω) in standard form

    Attributes:
        problem: Problem definition with rotor configuration
        mass_matrix_func: Function M(t, Ω) returning mass matrix
        mass_matrix: Cached mass matrix M evaluated at reference point
        mass_matrix_inv: Inverse of cached mass matrix M
        state_matrix_function: Function A(t, Ω) returning state matrix

    Example:
        >>> rotor = RotorBuild.build_all()
        >>> A = rotor.state_matrix_function(t=0.0, rotor_OMEGA=25.0)
        >>> eigenvalues = np.linalg.eigvals(A)
    """

    def __init__(self, problem: ProblemDefinition, auto_build: bool = True):
        """Initialize rotor build with problem definition.

        Args:
            problem: Problem definition containing rotor characteristics
            auto_build: If True, automatically build state matrix (default: True)
        """
        self.problem = problem
        self.mass_matrix_func: Callable[[float, float], np.ndarray] | None = None
        self.mass_matrix: np.ndarray | None = None
        self.mass_matrix_inv: np.ndarray | None = None
        self.state_matrix_function: Callable[[float, float], np.ndarray] | None = None
        
        if auto_build:
            self._initialize_matrices()
    
    def _initialize_matrices(self) -> None:
        """Initialize mass matrix and state matrix function.

        This private method:
        1. Gets the mass matrix function M(t, Ω)
        2. Evaluates and caches M at reference point (t=0, Ω=1)
        3. Computes and caches the inverse M^-1
        4. Creates the state matrix function A(t, Ω)
        """
        # Get mass matrix function handle
        (self.mass_matrix_func,self.mass_matrix_rot) = build_mass_matrix(self.problem)
        # Evaluate at reference point for cached inverse (M is constant in MBC for LTI)
        self.mass_matrix = self.mass_matrix_func(0.0,1.0)
        self.mass_matrix_inv = np.linalg.inv(self.mass_matrix)
        self.state_matrix_function = self._build_state_matrix_function()

    @staticmethod
    def build_all() -> 'RotorBuild':
        """Build rotor from default problem definition.

        This factory method creates a fully initialized RotorBuild
        with default configuration.

        Returns:
            Fully constructed RotorBuild object with default settings
        """
        problem = create_default_problem()
        return RotorBuild(problem, auto_build=True)
    
    def get_mass_matrix(self) -> np.ndarray:
        """Get the mass matrix.

        Returns:
            Mass matrix M (n x n)

        Raises:
            RuntimeError: If matrices haven't been initialized
        """
        if self.mass_matrix is None:
            raise RuntimeError("Mass matrix not initialized. Call _initialize_matrices() first.")
        return self.mass_matrix

    def get_damping_stiffness_functions(
        self
    ) -> tuple[Callable[[float, float], np.ndarray], Callable[[float, float], np.ndarray]]:
        """Get damping and stiffness matrix functions.

        Returns:
            Tuple of (C_func, K_func) where:
                - C_func(t, Ω): Returns damping matrix C
                - K_func(t, Ω): Returns stiffness matrix K
        """
        return build_damping_matrix(self.problem)

    def _build_state_matrix_function(self) -> Callable[[float, float], np.ndarray]:
        """Build state-space matrix function A(t, Ω).

        The state matrix has the standard first-order form:
            A = [[-M⁻¹C,  -M⁻¹K],
                 [  I,      0   ]]

        where:
            - M: Mass matrix (n x n, constant)
            - C: Damping matrix (n x n, function of t and Ω)
            - K: Stiffness matrix (n x n, function of t and Ω)
            - I: Identity matrix (n x n)
            - 0: Zero matrix (n x n)

        Returns:
            Function A(t, Ω) that returns the state matrix (2n x 2n)
        """
        # Use cached mass matrix inverse
        M_inv = self.mass_matrix_inv
        if M_inv is None:
            raise RuntimeError("Mass matrix inverse not initialized")

        # Get damping and stiffness function handles
        C_func, K_func = build_damping_matrix(self.problem, mass_matrix_rot=self.mass_matrix_rot)
        # Get system size
        n = M_inv.shape[0]

        def state_matrix_A(t: float, rotor_OMEGA: float) -> np.ndarray:
            """Compute state matrix at given time and rotor speed.

            Args:
                t: Time (seconds)
                rotor_OMEGA: Rotor angular velocity (rad/s)

            Returns:
                State matrix A (2n x 2n)
            """
            # Evaluate time-varying matrices
            C = C_func(t, rotor_OMEGA)
            K = K_func(t, rotor_OMEGA)

            # Build state matrix blocks
            # Upper row: [-M⁻¹C, -M⁻¹K]
            upper_left = -M_inv @ C
            upper_right = -M_inv @ K

            # Lower row: [I, 0]
            lower_left = np.eye(n)
            lower_right = np.zeros((n, n))

            # Assemble full state matrix
            A = np.block([
                [upper_left, upper_right],
                [lower_left, lower_right]
            ])

            return A

        return state_matrix_A


if __name__ == "__main__":
    # Test rotor build
    print("Building rotor system...")
    rotor = RotorBuild.build_all()

    print(f"\nProblem configuration:")
    print(f"  Number of blades: {rotor.problem.number_blades}")
    print(f"  Solver: {rotor.problem.required_solver}")
    print(f"  Damper connection: {rotor.problem.damper_connection}")

    # Test mass matrix
    M = rotor.get_mass_matrix()
    print(f"\nMass matrix shape: {M.shape}")
    print(f"Mass matrix condition number: {np.linalg.cond(M):.2e}")

    # Test state matrix evaluation
    if rotor.state_matrix_function is not None:
        A = rotor.state_matrix_function(t=0.0, rotor_OMEGA=25.0)
        print(f"\nState matrix A shape: {A.shape}")
        print(f"State matrix A condition number: {np.linalg.cond(A):.2e}")
        print(f"State matrix structure verified: {A.shape[0] == 2 * M.shape[0]}")
