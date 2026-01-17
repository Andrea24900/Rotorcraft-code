"""Rotor build class for constructing state-space matrices.

Converted from MATLAB rotor_build.m classdef.
Builds state-space representation for rotorcraft dynamics analysis.
"""

import numpy as np
from typing import Callable
from .config.problem_definition import ProblemDefinition, create_default_problem
from .utils.matrix_generation_GR import build_MCK_matrices


class RotorBuild:
    """Main rotor class for building state-space matrices.

    Constructs the state-space representation A(t, Ω) = [[-M⁻¹C, -M⁻¹K], [I, 0]]
    where M, C, K are computed together sharing MBC intermediate matrices.

    Attributes:
        problem: Problem definition with rotor configuration
        MCK_func: Function (M, C, K)(t, Ω) returning all system matrices
        mass_matrix: Cached mass matrix M evaluated at reference point
        state_matrix_function: Function A(t, Ω) returning state matrix
    """

    def __init__(self, problem: ProblemDefinition, auto_build: bool = True):
        """Initialize rotor build with problem definition.

        Args:
            problem: Problem definition containing rotor characteristics
            auto_build: If True, automatically build state matrix (default: True)
        """
        self.problem = problem
        self.MCK_func: Callable[[float, float], tuple] | None = None
        self.state_matrix_function: Callable[[float, float], np.ndarray] | None = None

        if auto_build:
            self._initialize_matrices()

    def _initialize_matrices(self) -> None:
        """Initialize MCK function and state matrix function.

        This private method:
        1. Gets the MCK matrix function (M, C, K)(t, Ω)
        2. Creates the state matrix function A(t, Ω)
        """
        self.MCK_func = build_MCK_matrices(self.problem)
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

    def _build_state_matrix_function(self) -> Callable[[float, float], np.ndarray]:
        """Build state-space matrix function A(t, Ω).

        The state matrix has the standard first-order form:
            A = [[-M⁻¹C,  -M⁻¹K],
                 [  I,      0   ]]

        Returns:
            Function A(t, Ω) that returns the state matrix (2n x 2n)
        """
        MCK_func = self.MCK_func
        if MCK_func is None:
            raise RuntimeError("MCK function not initialized")

        def state_matrix_A(t: float, rotor_OMEGA: float) -> np.ndarray:
            """Compute state matrix at given time and rotor speed."""
            # Get M, C, K in one call (shared MBC computation)
            M, C, K = MCK_func(t, rotor_OMEGA)
            M_inv = np.linalg.inv(M)
            n = M.shape[0]

            # Assemble state matrix A = [[-M⁻¹C, -M⁻¹K], [I, 0]]
            A = np.block([
                [-M_inv @ C, -M_inv @ K],
                [np.eye(n), np.zeros((n, n))]
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

    # Test state matrix evaluation
    if rotor.state_matrix_function is not None:
        A = rotor.state_matrix_function(t=0.0, rotor_OMEGA=25.0)
        print(f"\nState matrix A shape: {A.shape}")
        print(f"State matrix A condition number: {np.linalg.cond(A):.2e}")
