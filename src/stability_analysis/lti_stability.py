"""LTI stability analysis class."""

import numpy as np
from typing import Dict
from .base import StabilityAnalysis


class LTIStability(StabilityAnalysis):
    """Linear Time-Invariant stability analysis.

    Eigenvalues λ = σ + jω of the state matrix A determine stability:
    - σ < 0: stable, σ > 0: unstable
    - ω: oscillation frequency (rad/s)

    Returns from eigen_single_point:
        {'eigenvectors': ndarray, 'eigenvalues': ndarray}
    """

    def __init__(self, rotor_build) -> None:
        """Initialize LTI stability analysis."""
        super().__init__(rotor_build)
    
    def eigen_single_point(self, OMEGA: float) -> Dict[str, np.ndarray]:
        """Compute eigenvalues and eigenvectors at single operating point.

        Args:
            OMEGA: Rotor speed (rad/s). Must be positive.

        Returns:
            Dict with 'eigenvectors' (n,n) and 'eigenvalues' (n,).
            Eigenvalues are complex: real=damping, imag=frequency.

        Raises:
            TypeError: If OMEGA is not numeric.
            ValueError: If OMEGA <= 0.
            RuntimeError: If computation fails.
        """
        # Input validation
        if not isinstance(OMEGA, (int, float, np.number)):
            raise TypeError(f"OMEGA must be numeric, got {type(OMEGA).__name__}")

        if OMEGA <= 0:
            raise ValueError(f"OMEGA must be positive, got {OMEGA}")

        # Compute state matrix
        try:
            A = self.rotor_build.state_matrix_function(1, OMEGA)
        except Exception as e:
            raise RuntimeError(
                f"Failed to compute state matrix at OMEGA={OMEGA}: {e}"
            ) from e

        # Validate state matrix
        if not isinstance(A, np.ndarray):
            raise ValueError(
                f"State matrix must be numpy array, got {type(A).__name__}"
            )

        if A.ndim != 2:
            raise ValueError(
                f"State matrix must be 2D array, got {A.ndim}D array"
            )

        if A.shape[0] != A.shape[1]:
            raise ValueError(
                f"State matrix must be square, got shape {A.shape}"
            )

        # Compute eigenvalues and eigenvectors
        try:
            eigenvalues, eigenvectors = np.linalg.eig(A)
        except np.linalg.LinAlgError as e:
            raise RuntimeError(
                f"Eigenvalue computation failed at OMEGA={OMEGA}: {e}"
            ) from e

        eigensolution = {
            'eigenvectors': eigenvectors,
            'eigenvalues': eigenvalues
        }

        return eigensolution
    
    def eigen_full_range(self) -> 'LTIStability':
        """Compute eigenvalues over full RPM range.

        Populates modal_solution with damping, frequency, and eigensolution
        at each operating point.

        Returns:
            Self with populated modal_solution list.

        Raises:
            RuntimeError: If eigenvalue computation fails at any point.
        """
        # Initialize OMEGA range and modal_solution list
        self = self.assign_range_OMEGA()

        # Compute eigenvalues at each operating point
        for i in range(self.rotor_build.problem.number_points):
            try:
                eigensolution = self.eigen_single_point(self.modal_solution[i].OMEGA)

                self.modal_solution[i].damping = np.real(eigensolution['eigenvalues'])
                self.modal_solution[i].frequency = np.imag(eigensolution['eigenvalues'])
                self.modal_solution[i].eigensolution = eigensolution

            except (ValueError, RuntimeError) as e:
                raise RuntimeError(
                    f"Eigenvalue computation failed at point {i+1}/{self.rotor_build.problem.number_points}, "
                    f"OMEGA={self.modal_solution[i].OMEGA:.2f} rad/s "
                    f"({self.modal_solution[i].OMEGA_RPM:.2f} RPM): {e}"
                ) from e

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
    if len(lti.modal_solution) > 0:
        print(f"Number of eigenvalues: {len(lti.modal_solution[0].damping)}")
        print(f"Sample eigenvalues at first point:")
        print(f"  Damping: {lti.modal_solution[0].damping[:3]}")
        print(f"  Frequency: {lti.modal_solution[0].frequency[:3]}")
