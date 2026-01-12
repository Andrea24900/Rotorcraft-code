"""Continuation analysis for tracking eigenvalues across parameter space.

Converted from MATLAB continuation_analysis.m classdef.
This module implements predictor-corrector continuation methods for
bifurcation tracking in rotor dynamics stability analysis.
"""

import numpy as np
from typing import Callable, Tuple
from .base import StabilityAnalysis


class ContinuationAnalysis(StabilityAnalysis):
    """Continuation analysis subclass for eigenvalue tracking.

    Uses predictor-corrector methods to track eigenvalue branches
    as parameters (typically rotor speed) vary. This avoids mode
    crossing issues that can occur with standard eigenvalue analysis.
    """

    def __init__(self, rotor_build):
        """Initialize continuation analysis.

        Args:
            rotor_build: RotorBuild object
        """
        super().__init__(rotor_build)

    @staticmethod
    def _build_augmented_matrix(
        eigenvector: np.ndarray,
        eigenvalue: complex,
        A: np.ndarray,
        problem_size: int
    ) -> np.ndarray:
        """Build augmented system matrix for continuation.

        Constructs the augmented matrix used in both predictor (sensitivity)
        and corrector steps:
        [v,  λI - A]
        [0,  2v^H  ]

        Args:
            eigenvector: Current eigenvector (n×1)
            eigenvalue: Current eigenvalue (scalar)
            A: System matrix (n×n)
            problem_size: Size of the system

        Returns:
            Augmented matrix ((n+1)×(n+1))
        """
        return np.vstack([
            np.hstack([eigenvector.reshape(-1, 1),
                      eigenvalue * np.eye(problem_size) - A]),
            np.hstack([np.array([[0]]),
                      2 * eigenvector.conj().reshape(1, -1)])
        ])

    @staticmethod
    def _compute_residual(
        eigenvector: np.ndarray,
        eigenvalue: complex,
        A: np.ndarray
    ) -> np.ndarray:
        """Compute residual for corrector iterations.

        The residual measures how well the current iterate satisfies:
        1. Eigenvalue equation: Av = λv
        2. Normalization constraint: v^H·v = 1

        Args:
            eigenvector: Current eigenvector iterate
            eigenvalue: Current eigenvalue iterate
            A: System matrix

        Returns:
            Residual vector ((n+1)×1)
        """
        return np.vstack([
            (-eigenvalue * eigenvector + A @ eigenvector).reshape(-1, 1),
            np.array([[1 - np.vdot(eigenvector, eigenvector)]])
        ])

    @staticmethod
    def _solve_augmented_system(
        A_matrix: np.ndarray,
        b_vector: np.ndarray
    ) -> np.ndarray:
        """Solve augmented linear system with conditioning check.

        Uses direct solve for well-conditioned systems, pseudo-inverse
        for ill-conditioned systems.

        Args:
            A_matrix: Augmented system matrix
            b_vector: Right-hand side vector

        Returns:
            Solution vector
        """
        if np.linalg.cond(A_matrix) > 1/np.finfo(float).eps:
            return np.linalg.pinv(A_matrix) @ b_vector
        else:
            return np.linalg.solve(A_matrix, b_vector)

    def continuation(self) -> 'ContinuationAnalysis':
        """Run continuation analysis across parameter range.

        Returns:
            Self with computed modal solutions
        """
        # Initialize OMEGA range
        self = self.assign_range_OMEGA()

        # Assign periods for LTP solver
        if self.rotor_build.problem.required_solver == "LTP":
            self = self.assign_period_T()

        solver_type = self.rotor_build.problem.required_solver

        # STEP 1: Compute initial eigenvalue problem at first operating point
        print(f"Starting continuation analysis with {solver_type} solver...")

        if solver_type == "LTI":
            A, dM_dOMEGA, problem_size = self._setup_LTI_initial()

        elif solver_type == "LTP":
            A, dM_dOMEGA, problem_size = self._setup_LTP_initial()

        elif solver_type == "HD":
            A, dM_dOMEGA, problem_size = self._setup_HD_initial()
        else:
            raise ValueError(f"Unknown solver type: {solver_type}")

        # Solve initial eigenvalue problem
        eigenvalues, eigenvectors = np.linalg.eig(A)

        # Sort eigenvalues
        indices = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[indices]
        eigenvectors = eigenvectors[:, indices]

        # Store initial solution
        self.modal_solution[0].eigensolution = {
            'eigenvalues': eigenvalues.copy(),
            'eigenvectors': eigenvectors.copy()
        }

        # Extract damping and frequency
        self._extract_damping_frequency(0)

        # STEP 1.5: Predict next eigenvalues using sensitivity
        delta_OMEGA = abs(self.modal_solution[1].OMEGA - self.modal_solution[0].OMEGA)

        next_eigenvalue = np.zeros(problem_size, dtype=complex)
        next_eigenvector = np.zeros((problem_size, problem_size), dtype=complex)

        for j in range(problem_size):
            current_eigenvalue = self.modal_solution[0].eigensolution['eigenvalues'][j]
            current_eigenvector = self.modal_solution[0].eigensolution['eigenvectors'][:, j]

            # Build sensitivity system
            A_sens = self._build_augmented_matrix(
                current_eigenvector, current_eigenvalue, A, problem_size
            )

            b_sens = np.vstack([
                (dM_dOMEGA @ current_eigenvector).reshape(-1, 1),
                np.array([[0]])
            ])

            # Solve sensitivity system
            x_sens = self._solve_augmented_system(A_sens, b_sens)

            # Predict next values
            next_eigenvalue[j] = x_sens[0, 0] * delta_OMEGA + current_eigenvalue
            next_eigenvector[:, j] = x_sens[1:, 0] * delta_OMEGA + current_eigenvector

        # STEP 2: Iterate over remaining operating points
        for i in range(1, self.rotor_build.problem.number_points):
            print(f"  Point {i+1}/{self.rotor_build.problem.number_points}", end='\r')

            # Compute matrix and sensitivity at current operating point
            if solver_type == "LTI":
                A, dM_dOMEGA = self._compute_LTI_matrices(i)

            elif solver_type == "LTP":
                A, dM_dOMEGA = self._compute_LTP_matrices(i)

            elif solver_type == "HD":
                A, dM_dOMEGA, problem_size = self._compute_HD_matrices(i)

            # Initialize storage for this operating point
            self.modal_solution[i].eigensolution = {
                'eigenvalues': np.zeros(problem_size, dtype=complex),
                'eigenvectors': np.zeros((problem_size, problem_size), dtype=complex)
            }

            # Iterate to correct each predicted eigenvalue
            for j in range(problem_size):
                iter_eigenvalue = next_eigenvalue[j]
                iter_eigenvector = next_eigenvector[:, j]

                # Compute initial residual
                residual_b = self._compute_residual(iter_eigenvector, iter_eigenvalue, A)

                iter_norm = np.linalg.norm(np.hstack([iter_eigenvalue, iter_eigenvector]))
                iter_number = 0

                # Corrector iterations
                while (np.linalg.norm(residual_b) / iter_norm >
                       self.rotor_build.problem.continuation_tolerance and
                       iter_number < self.rotor_build.problem.continuation_max_iter):

                    # Build corrector system
                    A_con = self._build_augmented_matrix(
                        iter_eigenvector, iter_eigenvalue, A, problem_size
                    )

                    # Solve correction system
                    x_con = self._solve_augmented_system(A_con, residual_b)

                    # Update iterative values
                    iter_eigenvalue = x_con[0, 0] + iter_eigenvalue
                    iter_eigenvector = x_con[1:, 0] + iter_eigenvector
                    iter_number += 1

                    # Recompute residual
                    residual_b = self._compute_residual(iter_eigenvector, iter_eigenvalue, A)

                # Store converged values
                self.modal_solution[i].eigensolution['eigenvalues'][j] = iter_eigenvalue
                self.modal_solution[i].eigensolution['eigenvectors'][:, j] = iter_eigenvector

                # Predict next eigenvalue for next operating point (if not last point)
                if i < self.rotor_build.problem.number_points - 1:
                    current_eigenvalue = iter_eigenvalue
                    current_eigenvector = iter_eigenvector

                    # Build sensitivity system
                    A_sens = self._build_augmented_matrix(
                        current_eigenvector, current_eigenvalue, A, problem_size
                    )

                    b_sens = np.vstack([
                        (dM_dOMEGA @ current_eigenvector).reshape(-1, 1),
                        np.array([[0]])
                    ])

                    # Solve sensitivity system
                    x_sens = self._solve_augmented_system(A_sens, b_sens)

                    # Predict next values
                    next_eigenvalue[j] = x_sens[0, 0] * delta_OMEGA + current_eigenvalue
                    next_eigenvector[:, j] = x_sens[1:, 0] * delta_OMEGA + current_eigenvector

            # Extract damping and frequency for this operating point
            self._extract_damping_frequency(i)

        print(f"\n  Continuation completed!")
        return self

    def _setup_LTI_initial(self) -> Tuple[np.ndarray, np.ndarray, int]:
        """Setup initial LTI problem."""
        A_OMEGA = lambda omega: self.rotor_build.state_matrix_function(0, omega)
        A = A_OMEGA(self.modal_solution[0].OMEGA)
        problem_size = A.shape[0]

        dM_dOMEGA = self.numerical_sensitivity(
            'SS', A_OMEGA,
            self.rotor_build.problem.step_h,
            self.modal_solution[0].OMEGA
        )

        return A, dM_dOMEGA, problem_size

    def _setup_LTP_initial(self) -> Tuple[np.ndarray, np.ndarray, int]:
        """Setup initial LTP problem."""
        T = 2 * np.pi / self.modal_solution[0].OMEGA
        A_time_OMEGA = lambda t, omega: self.rotor_build.state_matrix_function(t, omega)
        A_time = lambda t: A_time_OMEGA(t, self.modal_solution[0].OMEGA)

        A = self.monodromy_computer(A_time, T)
        problem_size = A.shape[0]

        dM_dOMEGA = self.numerical_sensitivity(
            'MON', A_time_OMEGA,
            self.rotor_build.problem.step_h,
            self.modal_solution[0].OMEGA
        )

        return A, dM_dOMEGA, problem_size

    def _setup_HD_initial(self) -> Tuple[np.ndarray, np.ndarray, int]:
        """Setup initial HD problem."""
        T = 2 * np.pi / self.modal_solution[0].OMEGA
        time = np.linspace(0, T, self.rotor_build.problem.time_samples)

        A_time_OMEGA = lambda t, omega: self.rotor_build.state_matrix_function(t, omega)
        A_time = lambda t: A_time_OMEGA(t, self.modal_solution[0].OMEGA)

        A = self.hd_computer(A_time, time,
                            self.rotor_build.problem.number_harmonics,
                            self.modal_solution[0].OMEGA)
        problem_size = A.shape[0]

        dM_dOMEGA = self.numerical_sensitivity(
            'HD', A_time_OMEGA,
            self.rotor_build.problem.step_h,
            self.modal_solution[0].OMEGA,
            wrt_omega=True,
            HD_parameters=[self.rotor_build.problem.time_samples,
                          self.rotor_build.problem.number_harmonics]
        )

        return A, dM_dOMEGA, problem_size

    def _compute_LTI_matrices(self, i: int) -> Tuple[np.ndarray, np.ndarray]:
        """Compute LTI matrices at operating point i."""
        A_OMEGA = lambda omega: self.rotor_build.state_matrix_function(0, omega)
        A = A_OMEGA(self.modal_solution[i].OMEGA)

        dM_dOMEGA = self.numerical_sensitivity(
            'SS', A_OMEGA,
            self.rotor_build.problem.step_h,
            self.modal_solution[i].OMEGA
        )

        return A, dM_dOMEGA

    def _compute_LTP_matrices(self, i: int) -> Tuple[np.ndarray, np.ndarray]:
        """Compute LTP matrices at operating point i."""
        T = 2 * np.pi / self.modal_solution[i].OMEGA
        A_time_OMEGA = lambda t, omega: self.rotor_build.state_matrix_function(t, omega)
        A_time = lambda t: A_time_OMEGA(t, self.modal_solution[i].OMEGA)

        A = self.monodromy_computer(A_time, T)

        dM_dOMEGA = self.numerical_sensitivity(
            'MON', A_time_OMEGA,
            self.rotor_build.problem.step_h,
            self.modal_solution[i].OMEGA
        )

        return A, dM_dOMEGA

    def _compute_HD_matrices(self, i: int) -> Tuple[np.ndarray, np.ndarray, int]:
        """Compute HD matrices at operating point i."""
        T = 2 * np.pi / self.modal_solution[i].OMEGA
        time = np.linspace(0, T, self.rotor_build.problem.time_samples)

        A_time_OMEGA = lambda t, omega: self.rotor_build.state_matrix_function(t, omega)
        A_time = lambda t: A_time_OMEGA(t, self.modal_solution[i].OMEGA)

        A = self.hd_computer(A_time, time,
                            self.rotor_build.problem.number_harmonics,
                            self.modal_solution[i].OMEGA)
        problem_size = A.shape[0]

        dM_dOMEGA = self.numerical_sensitivity(
            'HD', A_time_OMEGA,
            self.rotor_build.problem.step_h,
            self.modal_solution[i].OMEGA,
            wrt_omega=True,
            HD_parameters=[self.rotor_build.problem.time_samples,
                          self.rotor_build.problem.number_harmonics]
        )

        return A, dM_dOMEGA, problem_size

    def _extract_damping_frequency(self, i: int):
        """Extract damping and frequency from eigenvalues at index i."""
        eigenvalues = self.modal_solution[i].eigensolution['eigenvalues']

        if self.rotor_build.problem.required_solver == "LTP":
            T = self.modal_solution[i].T
            # Floquet multipliers to exponents
            self.modal_solution[i].damping = (1 / (2 * T)) * np.log(
                np.real(eigenvalues)**2 + np.imag(eigenvalues)**2
            )
            self.modal_solution[i].frequency = (1 / T) * np.arctan2(
                np.imag(eigenvalues), np.real(eigenvalues)
            )
        else:
            # LTI or HD: eigenvalues are already in continuous time
            self.modal_solution[i].damping = np.real(eigenvalues)
            self.modal_solution[i].frequency = np.imag(eigenvalues)

    @staticmethod
    def _centered_difference(f_plus_h: np.ndarray, f_min_h: np.ndarray, h: float) -> np.ndarray:
        """Compute centered finite difference: (f(x+h) - f(x-h)) / (2h).

        Args:
            f_plus_h: Function evaluated at x + h
            f_min_h: Function evaluated at x - h
            h: Step size

        Returns:
            Centered difference approximation of derivative
        """
        return (f_plus_h - f_min_h) / (2 * h)

    @staticmethod
    def numerical_sensitivity(
        matrix_type: str,
        A_handle_all: Callable,
        step_size_h: float,
        OMEGA: float,
        wrt_omega: bool = True,
        par: float = None,
        HD_parameters: list = None
    ) -> np.ndarray:
        """Compute numerical sensitivity matrix for continuation.

        Uses centered finite difference approximation for derivatives.

        Args:
            matrix_type: 'SS' (state space), 'MON' (monodromy), or 'HD' (harmonic decomposition)
            A_handle_all: Function handle for state matrix
            step_size_h: Step size for finite differences
            OMEGA: Angular frequency
            wrt_omega: If True, compute derivative w.r.t. OMEGA; if False, w.r.t. par
            par: Parameter value (required if wrt_omega=False)
            HD_parameters: [time_samples, number_harmonics] for HD method

        Returns:
            Sensitivity matrix dM/d(parameter)
        """
        # Validate parameter requirement for parametric continuation
        if not wrt_omega and par is None:
            raise ValueError("par must be provided when wrt_omega=False")

        if matrix_type == 'SS':
            # State-space sensitivity
            if wrt_omega:
                # Variation with respect to OMEGA
                A_plus_h = A_handle_all(OMEGA + step_size_h)
                A_min_h = A_handle_all(OMEGA - step_size_h)
            else:
                # Variation with respect to parameter
                A_plus_h = A_handle_all(par + step_size_h, OMEGA)
                A_min_h = A_handle_all(par - step_size_h, OMEGA)

            sensitivity = ContinuationAnalysis._centered_difference(A_plus_h, A_min_h, step_size_h)

        elif matrix_type == 'HD':
            # Harmonic decomposition sensitivity
            number_time_instants = HD_parameters[0]
            number_harmonics = HD_parameters[1]

            if wrt_omega:
                # Variation with respect to OMEGA
                time_min_h = np.linspace(0, 2 * np.pi / (OMEGA - step_size_h), number_time_instants)
                A_min_h = lambda t: A_handle_all(t, OMEGA - step_size_h)
                A_HD_min_h = StabilityAnalysis.hd_computer(
                    A_min_h, time_min_h, number_harmonics, OMEGA - step_size_h
                )

                time_plus_h = np.linspace(0, 2 * np.pi / (OMEGA + step_size_h), number_time_instants)
                A_plus_h = lambda t: A_handle_all(t, OMEGA + step_size_h)
                A_HD_plus_h = StabilityAnalysis.hd_computer(
                    A_plus_h, time_plus_h, number_harmonics, OMEGA + step_size_h
                )
            else:
                # Variation with respect to parameter
                T = 2 * np.pi / OMEGA
                time = np.linspace(0, T, number_time_instants)

                A_plus_h = lambda t: A_handle_all(t, par + step_size_h, OMEGA)
                A_HD_plus_h = StabilityAnalysis.hd_computer(
                    A_plus_h, time, number_harmonics, OMEGA
                )

                A_min_h = lambda t: A_handle_all(t, par - step_size_h, OMEGA)
                A_HD_min_h = StabilityAnalysis.hd_computer(
                    A_min_h, time, number_harmonics, OMEGA
                )

            sensitivity = ContinuationAnalysis._centered_difference(A_HD_plus_h, A_HD_min_h, step_size_h)

        elif matrix_type == 'MON':
            # Monodromy matrix sensitivity
            if wrt_omega:
                # Variation with respect to OMEGA
                T_plus_h = 2 * np.pi / (OMEGA + step_size_h)
                A_plus_h = lambda t: A_handle_all(t, OMEGA + step_size_h)
                M_plus_h = StabilityAnalysis.monodromy_computer(A_plus_h, T_plus_h)

                T_min_h = 2 * np.pi / (OMEGA - step_size_h)
                A_min_h = lambda t: A_handle_all(t, OMEGA - step_size_h)
                M_min_h = StabilityAnalysis.monodromy_computer(A_min_h, T_min_h)
            else:
                # Variation with respect to parameter
                T = 2 * np.pi / OMEGA

                A_plus_h = lambda t: A_handle_all(t, par + step_size_h, OMEGA)
                M_plus_h = StabilityAnalysis.monodromy_computer(A_plus_h, T)

                A_min_h = lambda t: A_handle_all(t, par - step_size_h, OMEGA)
                M_min_h = StabilityAnalysis.monodromy_computer(A_min_h, T)

            sensitivity = ContinuationAnalysis._centered_difference(M_plus_h, M_min_h, step_size_h)

        else:
            raise ValueError(f"Unknown matrix type: {matrix_type}")

        return sensitivity
