"""LTP stability analysis class using Floquet theory."""

import numpy as np
from typing import Dict
from .base import StabilityAnalysis


class LTPStability(StabilityAnalysis):
    """Linear Time-Periodic stability analysis using Floquet theory.

    For systems with dx/dt = A(t)x where A(t+T) = A(t).
    Computes monodromy matrix M = Φ(T) and extracts characteristic multipliers μ.

    Stability: |μ| < 1 stable, |μ| > 1 unstable.
    Exponents: σ = (1/2T)ln(|μ|²), ω = (1/T)atan2(Im(μ), Re(μ))

    Returns from ltp_single_point:
        {'char_multipliers': ndarray, 'real_char_exp': ndarray, 'imag_char_exp': ndarray}
    """

    def __init__(self, rotor_build) -> None:
        """Initialize LTP stability analysis."""
        super().__init__(rotor_build)
    
    def ltp_single_point(self, OMEGA: float) -> Dict[str, np.ndarray]:
        """Compute characteristic multipliers and exponents at single operating point.

        Args:
            OMEGA: Rotor speed (rad/s). Must be positive.

        Returns:
            Dict with 'char_multipliers' (n,), 'real_char_exp' (n,), 'imag_char_exp' (n,).
            real_char_exp is damping (σ < 0 stable), imag_char_exp is frequency.

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

        # Compute period from OMEGA
        T = 2 * np.pi / OMEGA

        # Create time-dependent state matrix handle A(OMEGA, t)
        try:
            A = lambda Om, t: self.rotor_build.state_matrix_function(t, Om)
        except Exception as e:
            raise RuntimeError(
                f"Failed to create state matrix handle at OMEGA={OMEGA}: {e}"
            ) from e

        # Compute monodromy matrix
        try:
            monodromy_matrix_M = self.monodromy_computer(
                A, OMEGA,
                rtol=self.rotor_build.problem.ode_rtol,
                atol=self.rotor_build.problem.ode_atol
            )
        except Exception as e:
            raise RuntimeError(
                f"Monodromy matrix computation failed at OMEGA={OMEGA:.2f} rad/s, "
                f"T={T:.6f}s: {e}"
            ) from e

        # Validate monodromy matrix
        if not isinstance(monodromy_matrix_M, np.ndarray):
            raise RuntimeError(
                f"Monodromy matrix must be numpy array, got {type(monodromy_matrix_M).__name__}"
            )

        if monodromy_matrix_M.ndim != 2 or monodromy_matrix_M.shape[0] != monodromy_matrix_M.shape[1]:
            raise RuntimeError(
                f"Monodromy matrix must be square 2D array, got shape {monodromy_matrix_M.shape}"
            )

        # Compute characteristic multipliers (eigenvalues of M)
        try:
            char_multipliers = np.linalg.eigvals(monodromy_matrix_M)
        except np.linalg.LinAlgError as e:
            raise RuntimeError(
                f"Failed to compute eigenvalues of monodromy matrix at "
                f"OMEGA={OMEGA:.2f} rad/s: {e}"
            ) from e

        # Compute characteristic exponents
        # Real part: σ = (1/2T) * ln(|μ|²)
        magnitude_squared = np.real(char_multipliers)**2 + np.imag(char_multipliers)**2

        # Avoid log of zero or negative values
        magnitude_squared = np.maximum(magnitude_squared, 1e-300)

        real_char_exp = (1 / (2 * T)) * np.log(magnitude_squared)

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

        Populates modal_solution with damping (real_char_exp), frequency (imag_char_exp),
        and char_solution at each operating point.

        Returns:
            Self with populated modal_solution list.

        Raises:
            RuntimeError: If computation fails at any point.
        """
        # Assign OMEGA range
        self.assign_range_OMEGA()

        # Compute characteristic solutions at each point
        for i in range(self.rotor_build.problem.number_points):
            try:
                char_solution = self.ltp_single_point(self.modal_solution[i].OMEGA)

                # Store results
                self.modal_solution[i].char_solution = char_solution
                self.modal_solution[i].damping = char_solution['real_char_exp']
                self.modal_solution[i].frequency = char_solution['imag_char_exp']

            except (ValueError, RuntimeError, TypeError) as e:
                raise RuntimeError(
                    f"LTP analysis failed at point {i+1}/{self.rotor_build.problem.number_points}, "
                    f"OMEGA={self.modal_solution[i].OMEGA:.2f} rad/s "
                    f"({self.modal_solution[i].OMEGA_RPM:.2f} RPM): {e}"
                ) from e

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
        print(f"Sample characteristic exponents at first point:")
        print(f"  Damping: {ltp.modal_solution[0].damping[:3]}")
        print(f"  Frequency: {ltp.modal_solution[0].frequency[:3]}")
