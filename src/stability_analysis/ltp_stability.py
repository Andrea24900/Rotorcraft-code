"""LTP stability analysis class using Floquet theory.

This module implements Linear Time-Periodic (LTP) stability analysis for rotor
systems with time-varying coefficients. It uses Floquet theory to determine
system stability through monodromy matrix computation.

Theoretical Background:
----------------------
For a linear time-periodic system: dx/dt = A(t)x, where A(t+T) = A(t)

Floquet theory states that solutions can be expressed as:
    x(t) = Φ(t)x₀, where Φ(t) is the transition matrix

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
from typing import Dict
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
    >>> ltp = ltp.ltp_full_range()
    >>> # Check stability
    >>> max_damping = max([np.max(sol.damping) for sol in ltp.modal_solution])
    >>> is_stable = max_damping < 0
    """

    def __init__(self, rotor_build) -> None:
        """Initialize LTP stability analysis.

        Args:
            rotor_build: RotorBuild object containing rotor configuration,
                        state matrix functions, and problem definition
        """
        super().__init__(rotor_build)
    
    def ltp_single_point(self, OMEGA: float, T: float) -> Dict[str, np.ndarray]:
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
            Rotor speed (rad/s). Must be positive.
        T : float
            Period of oscillation (s), typically T = 2π/OMEGA. Must be positive.

        Returns:
        -------
        Dict[str, np.ndarray]
            Dictionary containing:
            - 'char_multipliers': Complex array of shape (n,) with characteristic
              multipliers (μ). These are eigenvalues of monodromy matrix.
            - 'real_char_exp': Real array of shape (n,) with real part of
              characteristic exponents (σ) - damping indicators
            - 'imag_char_exp': Real array of shape (n,) with imaginary part of
              characteristic exponents (ω) - oscillation frequencies

        Raises:
        ------
        TypeError
            If OMEGA or T are not numeric
        ValueError
            If OMEGA <= 0 or T <= 0
        RuntimeError
            If monodromy matrix computation or eigenvalue extraction fails

        Notes:
        -----
        - Negative real exponents (σ < 0) indicate stable modes
        - Positive real exponents (σ > 0) indicate unstable modes
        - Imaginary exponents represent oscillation frequencies
        - For rotors, typically T = 2π/OMEGA (one revolution period)
        - Monodromy computation can be expensive for large systems
        - Uses inherited monodromy_computer from StabilityAnalysis base class

        Example:
        -------
        >>> ltp = LTPStability(rotor_build)
        >>> OMEGA = 25.0  # rad/s
        >>> T = 2 * np.pi / OMEGA
        >>> char_sol = ltp.ltp_single_point(OMEGA, T)
        >>> multipliers = char_sol['char_multipliers']
        >>> damping = char_sol['real_char_exp']
        >>> is_stable = np.all(damping < 0)
        """
        # Input validation
        if not isinstance(OMEGA, (int, float, np.number)):
            raise TypeError(f"OMEGA must be numeric, got {type(OMEGA).__name__}")

        if OMEGA <= 0:
            raise ValueError(f"OMEGA must be positive, got {OMEGA}")

        if not isinstance(T, (int, float, np.number)):
            raise TypeError(f"T must be numeric, got {type(T).__name__}")

        if T <= 0:
            raise ValueError(f"T (period) must be positive, got {T}")

        # Validate consistency between OMEGA and T
        expected_T = 2 * np.pi / OMEGA
        if not np.isclose(T, expected_T, rtol=0.1):
            import warnings
            warnings.warn(
                f"Period T={T:.6f}s does not match expected value "
                f"2π/OMEGA={expected_T:.6f}s (OMEGA={OMEGA:.4f} rad/s). "
                f"Ensure T is correct for your application.",
                UserWarning,
                stacklevel=2
            )

        # Create time-dependent state matrix handle
        try:
            A = lambda t: self.rotor_build.state_matrix_function(t, OMEGA)
        except Exception as e:
            raise RuntimeError(
                f"Failed to create state matrix handle at OMEGA={OMEGA}: {e}"
            ) from e

        # Compute monodromy matrix
        try:
            monodromy_matrix_M = self.monodromy_computer(A, T)
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

        This method performs a sweep of operating points across the specified
        RPM range, computing characteristic exponents at each point using
        Floquet theory. The results are stored in the modal_solution attribute.

        Process:
        -------
        1. Generate array of OMEGA values from problem definition using
           assign_range_OMEGA() method (inherited from base class)
        2. Compute corresponding periods T = 2π/OMEGA for each point
        3. For each operating point:
           - Compute monodromy matrix via numerical integration
           - Extract characteristic multipliers (eigenvalues of M)
           - Convert to characteristic exponents
           - Store damping and frequency values

        Returns:
        -------
        LTPStability
            Self with populated modal_solution list containing results at
            each operating point. Enables method chaining.

        Attributes Modified:
        -------------------
        modal_solution : list[ModalSolution]
            Populated with characteristic exponents at each RPM point.
            Length equals self.rotor_build.problem.number_points.

            Each entry contains:
            - OMEGA: Rotor speed (rad/s)
            - OMEGA_RPM: Rotor speed (RPM)
            - T: Period (s) = 2π/OMEGA
            - char_solution: Full solution dictionary with multipliers and exponents
            - damping: Real characteristic exponents (σ) - stability indicators
            - frequency: Imaginary characteristic exponents (ω) - oscillation frequencies

        Raises:
        ------
        RuntimeError
            If monodromy computation or characteristic exponent extraction
            fails at any operating point

        Notes:
        -----
        - This method can be computationally expensive for many points
        - Each point requires integration of n columns of the fundamental matrix
        - Each monodromy computation requires n ODE integrations over period T
        - For n DOF system, expect n complex conjugate pairs of exponents
        - Use continuation methods for more efficient tracking of specific modes
        - Method uses method chaining pattern for fluent interface
        - Integration tolerance controlled by base class defaults
          (DEFAULT_RTOL=1e-6, DEFAULT_ATOL=1e-8)

        Example:
        -------
        >>> rotor = RotorBuild.build_all()
        >>> rotor.problem.required_solver = "LTP"
        >>> rotor.problem.lower_rotor_RPM = 1000
        >>> rotor.problem.higher_rotor_RPM = 5000
        >>> rotor.problem.number_points = 20
        >>> ltp = LTPStability(rotor)
        >>> ltp = ltp.ltp_full_range()
        >>> # Check stability at all points
        >>> max_damping = max([np.max(sol.damping) for sol in ltp.modal_solution])
        >>> is_stable = max_damping < 0
        >>> print(f"System stable: {is_stable}")
        """
        # Assign OMEGA range and periods
        self.assign_range_OMEGA()
        self.assign_period_T()

        # Compute characteristic solutions at each point
        for i in range(self.rotor_build.problem.number_points):
            try:
                char_solution = self.ltp_single_point(
                    self.modal_solution[i].OMEGA,
                    self.modal_solution[i].T
                )

                # Store results
                self.modal_solution[i].char_solution = char_solution
                self.modal_solution[i].damping = char_solution['real_char_exp']
                self.modal_solution[i].frequency = char_solution['imag_char_exp']

            except (ValueError, RuntimeError, TypeError) as e:
                raise RuntimeError(
                    f"LTP analysis failed at point {i+1}/{self.rotor_build.problem.number_points}, "
                    f"OMEGA={self.modal_solution[i].OMEGA:.2f} rad/s "
                    f"({self.modal_solution[i].OMEGA_RPM:.2f} RPM), "
                    f"T={self.modal_solution[i].T:.6f}s: {e}"
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
