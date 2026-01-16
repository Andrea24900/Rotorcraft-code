"""HD (Harmonic Decomposition) stability analysis using Lopez and Prasad method."""

import numpy as np
from typing import Dict
from .base import StabilityAnalysis


class HDStability(StabilityAnalysis):
    """Harmonic Decomposition (HD) stability analysis.

    Alternative to LTP for time-periodic systems. Expands A(t) using Fourier series
    and solves an augmented algebraic eigenvalue problem (no ODE integration).

    HD matrix size: (1 + 2N) × n where N = number_harmonics, n = state dimension.
    Eigenvalues λ = σ + jω directly indicate stability (σ < 0 stable).

    Returns from hd_single_point:
        {'eigenvectors': ndarray, 'eigenvalues': ndarray}

    Requires problem attributes: number_harmonics, time_samples, hd_use_complex.
    """

    def __init__(self, rotor_build) -> None:
        """Initialize HD stability analysis."""
        super().__init__(rotor_build)
    
    def hd_single_point(self, OMEGA: float) -> Dict[str, np.ndarray]:
        """Compute eigenvalues using Harmonic Decomposition method.

        Args:
            OMEGA: Rotor speed (rad/s). Must be positive.

        Returns:
            Dict with 'eigenvectors' (n_hd, n_hd) and 'eigenvalues' (n_hd,).
            Eigenvalues are complex: real=damping, imag=frequency.

        Raises:
            TypeError: If OMEGA is not numeric.
            ValueError: If OMEGA <= 0 or config invalid.
            AttributeError: If required problem attributes missing.
            RuntimeError: If computation fails.
        """
        # Input validation
        if not isinstance(OMEGA, (int, float, np.number)):
            raise TypeError(f"OMEGA must be numeric, got {type(OMEGA).__name__}")

        if OMEGA <= 0:
            raise ValueError(f"OMEGA must be positive, got {OMEGA}")

        # Compute period from OMEGA
        T = 2 * np.pi / OMEGA

        # Validate problem configuration
        if not hasattr(self.rotor_build.problem, 'number_harmonics'):
            raise AttributeError(
                "rotor_build.problem must have 'number_harmonics' attribute"
            )

        if not hasattr(self.rotor_build.problem, 'time_samples'):
            raise AttributeError(
                "rotor_build.problem must have 'time_samples' attribute"
            )

        if not hasattr(self.rotor_build.problem, 'hd_use_complex'):
            raise AttributeError(
                "rotor_build.problem must have 'hd_use_complex' attribute"
            )

        number_harmonics = self.rotor_build.problem.number_harmonics
        time_samples = self.rotor_build.problem.time_samples

        # Validate configuration values
        if not isinstance(number_harmonics, (int, np.integer)) or number_harmonics < 1:
            raise ValueError(
                f"number_harmonics must be integer >= 1, got {number_harmonics}"
            )

        if not isinstance(time_samples, (int, np.integer)) or time_samples < 10:
            raise ValueError(
                f"time_samples must be integer >= 10, got {time_samples}"
            )

        # Create time-dependent state matrix handle
        try:
            A = lambda t: self.rotor_build.state_matrix_function(t, OMEGA)
        except Exception as e:
            raise RuntimeError(
                f"Failed to create state matrix handle at OMEGA={OMEGA}: {e}"
            ) from e

        # Time vector for sampling using problem definition
        try:
            time = np.linspace(0, T, time_samples)
        except Exception as e:
            raise RuntimeError(
                f"Failed to create time vector: {e}"
            ) from e

        # Compute Harmonic Decomposition matrix
        try:
            A_HD = StabilityAnalysis.hd_computer(
                A,
                time,
                number_harmonics,
                OMEGA,
                use_complex=self.rotor_build.problem.hd_use_complex
            )
        except Exception as e:
            raise RuntimeError(
                f"HD matrix computation failed at OMEGA={OMEGA:.2f} rad/s, "
                f"T={T:.6f}s with {number_harmonics} harmonics and "
                f"{time_samples} samples: {e}"
            ) from e

        # Validate HD matrix
        if not isinstance(A_HD, np.ndarray):
            raise RuntimeError(
                f"HD matrix must be numpy array, got {type(A_HD).__name__}"
            )

        if A_HD.ndim != 2:
            raise RuntimeError(
                f"HD matrix must be 2D array, got {A_HD.ndim}D array"
            )

        if A_HD.shape[0] != A_HD.shape[1]:
            raise RuntimeError(
                f"HD matrix must be square, got shape {A_HD.shape}"
            )

        # Compute eigenvalues and eigenvectors of expanded system
        try:
            eigenvalues, eigenvectors = np.linalg.eig(A_HD)
        except np.linalg.LinAlgError as e:
            raise RuntimeError(
                f"Eigenvalue computation failed for HD matrix at "
                f"OMEGA={OMEGA:.2f} rad/s: {e}"
            ) from e

        eigensolution = {
            'eigenvectors': eigenvectors,
            'eigenvalues': eigenvalues
        }

        return eigensolution
    
    def hd_full_range(self) -> 'HDStability':
        """Compute eigenvalues over full RPM range using HD method.

        Populates modal_solution with damping, frequency, and eigensolution
        at each operating point.

        Returns:
            Self with populated modal_solution list.

        Raises:
            RuntimeError: If computation fails at any point.
            AttributeError: If problem config missing required attributes.
        """
        # Assign OMEGA range
        self.assign_range_OMEGA()

        # Compute eigenvalues at each point
        for i in range(self.rotor_build.problem.number_points):
            try:
                eigensolution = self.hd_single_point(self.modal_solution[i].OMEGA)

                # Store results
                self.modal_solution[i].damping = np.real(eigensolution['eigenvalues'])
                self.modal_solution[i].frequency = np.imag(eigensolution['eigenvalues'])
                self.modal_solution[i].eigensolution = eigensolution

            except (ValueError, RuntimeError, TypeError, AttributeError) as e:
                raise RuntimeError(
                    f"HD analysis failed at point {i+1}/{self.rotor_build.problem.number_points}, "
                    f"OMEGA={self.modal_solution[i].OMEGA:.2f} rad/s "
                    f"({self.modal_solution[i].OMEGA_RPM:.2f} RPM): {e}"
                ) from e

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
    hd = hd.hd_full_range()
    
    print(f"Number of solutions: {len(hd.modal_solution)}")
    print(f"First point OMEGA: {hd.modal_solution[0].OMEGA_RPM:.2f} RPM")
    if len(hd.modal_solution) > 0:
        print(f"Number of eigenvalues: {len(hd.modal_solution[0].damping)}")
        print(f"Sample eigenvalues at first point:")
        print(f"  Damping: {hd.modal_solution[0].damping[:3]}")
        print(f"  Frequency: {hd.modal_solution[0].frequency[:3]}")
