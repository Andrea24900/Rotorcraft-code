"""LTI stability analysis class.

Converted from MATLAB LTI_stability.m classdef.
"""

import numpy as np
from typing import Dict
from .base import StabilityAnalysis


class LTIStability(StabilityAnalysis):
    """Linear Time-Invariant stability analysis.

    This class provides methods for stability analysis of Linear Time-Invariant
    (LTI) rotor systems. For LTI systems, the state matrix A does not depend on
    time, and eigenvalue analysis directly reveals system stability.

    The eigenvalues of the state matrix A are continuous-time poles of the system:
    - Real part (damping): σ < 0 indicates stability, σ > 0 indicates instability
    - Imaginary part (frequency): ω represents oscillation frequency (rad/s)
    - Complex conjugate pairs represent oscillatory modes

    Attributes:
    ----------
    rotor_build : RotorBuild
        Rotor system configuration containing state matrix functions
    modal_solution : list[ModalSolution]
        List of modal solutions at each operating point, each containing:
        - OMEGA: Rotor speed (rad/s)
        - OMEGA_RPM: Rotor speed (RPM)
        - damping: Real parts of eigenvalues (stability indicators)
        - frequency: Imaginary parts of eigenvalues (oscillation frequencies)
        - eigensolution: Dict with eigenvectors and eigenvalues

    Notes:
    -----
    - LTI analysis is computationally efficient compared to LTP/HD methods
    - Suitable for rotors with time-invariant coefficients
    - For time-periodic systems, use LTPStability or HDStability instead
    - Natural frequency: ωₙ = √(σ² + ω²)
    - Damping ratio: ζ = -σ / √(σ² + ω²)

    Example:
    -------
    >>> rotor = RotorBuild.build_all()
    >>> rotor.problem.required_solver = "LTI"
    >>> lti = LTIStability(rotor)
    >>> lti = lti.eigen_full_range()
    >>> # Check stability at each point
    >>> for sol in lti.modal_solution:
    >>>     max_damping = np.max(sol.damping)
    >>>     print(f"RPM={sol.OMEGA_RPM:.0f}, Stable={max_damping < 0}")
    """

    def __init__(self, rotor_build) -> None:
        """Initialize LTI stability analysis.

        Args:
            rotor_build: RotorBuild object containing rotor configuration,
                        state matrix functions, and problem definition
        """
        super().__init__(rotor_build)
    
    def eigen_single_point(self, OMEGA: float) -> Dict[str, np.ndarray]:
        """Compute eigenvalues and eigenvectors at single operating point.

        This method performs eigenvalue analysis of the LTI state matrix A
        at a specified rotor speed. The eigenvalues λ = σ + jω represent
        the continuous-time poles of the linearized rotor system.

        Args:
            OMEGA: Rotor speed (rad/s). Must be positive.

        Returns:
            Dictionary containing:
            - 'eigenvectors': np.ndarray of shape (n, n) containing eigenvectors
              as columns. Each column corresponds to a mode shape.
            - 'eigenvalues': np.ndarray of shape (n,) containing eigenvalues.
              Real part is damping, imaginary part is frequency.

        Raises:
            TypeError: If OMEGA is not numeric
            ValueError: If OMEGA <= 0
            RuntimeError: If state matrix computation or eigenvalue
                         decomposition fails

        Notes:
        -----
        - Eigenvalues represent continuous-time poles: λ = σ + jω
        - Stability condition: σ < 0 for all eigenvalues
        - Complex conjugate pairs represent oscillatory modes
        - Uses numpy.linalg.eig for general (non-symmetric) matrices
        - Time argument in state_matrix_function is set to 1 (arbitrary for LTI)

        Example:
        -------
        >>> lti = LTIStability(rotor_build)
        >>> eigensol = lti.eigen_single_point(25.0)  # 25 rad/s
        >>> eigenvalues = eigensol['eigenvalues']
        >>> damping = np.real(eigenvalues)
        >>> frequency = np.imag(eigenvalues)
        >>> is_stable = np.all(damping < 0)
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

        This method performs eigenvalue analysis at each operating point
        across the specified RPM range, storing results in modal_solution.
        The method assumes the system is Linear Time-Invariant (LTI), so
        eigenvalues represent continuous-time poles of the system.

        Process:
        -------
        1. Generate array of OMEGA values from problem definition using
           assign_range_OMEGA() method (inherited from base class)
        2. For each operating point:
           - Compute state matrix A(OMEGA) via rotor_build
           - Extract eigenvalues (poles) and eigenvectors (modes)
           - Store damping (real parts) and frequency (imaginary parts)

        Returns:
        -------
        LTIStability
            Self with populated modal_solution list containing results at
            each operating point.

            Each entry in modal_solution contains:
            - OMEGA: Rotor speed (rad/s)
            - OMEGA_RPM: Rotor speed (RPM)
            - damping: Real part of eigenvalues (σ) - negative = stable
            - frequency: Imaginary part of eigenvalues (ω) - oscillation frequency
            - eigensolution: Dict with full eigenvectors and eigenvalues

        Attributes Modified:
        -------------------
        modal_solution : list[ModalSolution]
            Populated with eigenvalue solutions at each RPM point.
            Length equals self.rotor_build.problem.number_points.

        Raises:
        ------
        RuntimeError
            If eigenvalue computation fails at any operating point

        Notes:
        -----
        - For stable rotor, all damping values should be negative
        - Complex conjugate eigenvalue pairs represent oscillatory modes
        - Damping ratio: ζ = -σ / √(σ² + ω²)
        - Natural frequency: ωₙ = √(σ² + ω²)
        - Method uses method chaining pattern for fluent interface
        - Significantly faster than LTP (monodromy) or HD methods

        Example:
        -------
        >>> rotor = RotorBuild.build_all()
        >>> rotor.problem.required_solver = "LTI"
        >>> rotor.problem.lower_rotor_RPM = 0
        >>> rotor.problem.higher_rotor_RPM = 1000
        >>> rotor.problem.number_points = 100
        >>> lti = LTIStability(rotor)
        >>> lti = lti.eigen_full_range()
        >>> # Check stability at all points
        >>> max_damping = max([np.max(sol.damping) for sol in lti.modal_solution])
        >>> is_stable = max_damping < 0
        >>> print(f"System stable: {is_stable}")
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
