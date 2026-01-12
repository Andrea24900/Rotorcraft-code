"""HD (Harmonic Decomposition) stability analysis using Lopez and Prasad method.

This module implements Harmonic Decomposition (HD) stability analysis for rotor
systems with time-periodic coefficients. The method expands the state-space system 
in a lifted one with the contributions of the various harmonics.

Theoretical Background:
----------------------
For a linear time-periodic system: dx/dt = A(t)x, where A(t+T) = A(t)

The Harmonic Decomposition method:
1. Expands x(t) using Fourier series: x(t) = x₀ + Σ[xₙcos(nωt) + xₙsin(nωt)]
   or complex Foutier series x(t) = Σx_n * exp(inomegat)
2. Multiplies the resulting equations for sinusoidal or complex exponential terms
3. Integrates over one period
4. Arranges the resulting terms to build the A_HD matrix, which follows LTI stability

The HD matrix size grows with number of harmonics:
    size(A_HD) = (1 + 2N) × n, where N = number of harmonics, n = state dimension

Advantages:
- Algebraic eigenvalue problem (no ODE integration like monodromy)
- Can be faster than LTP for moderate number of harmonics

Disadvantages:
- Matrix size grows rapidly with harmonics: O((2N+1)×n)²
- Truncation error from finite harmonics

Converted from MATLAB HD_stability.m classdef.
"""

import numpy as np
from typing import Dict
from .base import StabilityAnalysis
from .ltp_stability import LTPStability


class HDStability(LTPStability):
    """Harmonic Decomposition (HD) stability analysis.

    This class extends LTP stability to use the Harmonic Decomposition (HD) method
    for analyzing periodic systems. Instead of computing the monodromy matrix through
    numerical integration (as in LTP), HD expands the periodic system
    using Fourier series and solves an augmented algebraic eigenvalue problem.

    The method is based on the Harmonic Decomposition by Lopez and Prasad.

    Key Features:
    ------------
    - Frequency domain representation using Fourier expansion
    - Algebraic eigenvalue problem (no ODE integration required)
    - Configurable number of harmonics for accuracy vs. cost trade-off
    - Both real and complex formulations supported
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
        - eigensolution: Dict with eigenvectors and eigenvalues of expanded system
        - damping: Real part of eigenvalues (σ)
        - frequency: Imaginary part of eigenvalues (ω)

    Methods:
    -------
    hd_single_point(OMEGA, T)
        Compute eigenvalues at a single operating point using HD method
    hd_full_range()
        Compute eigenvalues over full RPM range using HD method

    Notes:
    -----
    - Inherits from LTPStability but uses different analysis method
    - HD matrix size: (1 + 2N) × n where N = harmonics, n = state size
    - More harmonics → better accuracy but larger matrices
    - Eigenvalues are direct stability indicators (λ = σ + jω)
    - Negative real part (σ < 0) indicates stable mode
    - Can be faster than LTP monodromy for small N
    - No ODE integration required (advantage over LTP)

    Example:
    -------
    >>> rotor = RotorBuild.build_all()
    >>> rotor.problem.required_solver = "HD"
    >>> rotor.problem.number_harmonics = 3
    >>> rotor.problem.time_samples = 100
    >>> hd = HDStability(rotor)
    >>> hd = hd.hd_full_range()
    >>> # Check stability
    >>> max_damping = max([np.max(sol.damping) for sol in hd.modal_solution])
    >>> is_stable = max_damping < 0
    """

    def __init__(self, rotor_build) -> None:
        """Initialize HD stability analysis.

        Args:
            rotor_build: RotorBuild object containing rotor configuration,
                        state matrix functions, and problem definition including:
                        - number_harmonics: Number of harmonics to use in expansion
                        - time_samples: Number of time samples for Fourier coefficients
                        - hd_use_complex: Whether to use complex formulation
        """
        super().__init__(rotor_build)
    
    def hd_single_point(self, OMEGA: float, T: float) -> Dict[str, np.ndarray]:
        """Compute eigenvalues using Harmonic Decomposition method.

        This method constructs the Harmonic Decomposition matrix by:
        1. Sampling the time-periodic state matrix A(t) over one period
        2. Computing Fourier coefficients via FFT
        3. Assembling the augmented HD matrix A_HD
        4. Solving the algebraic eigenvalue problem: A_HD × X = λ × X

        The eigenvalues λ directly represent system stability (no conversion needed).

        Algorithm:
        ---------
        1. Create time-dependent state matrix function A(t, OMEGA)
        2. Generate time vector with specified number of samples over period T
        3. Call hd_computer to construct augmented matrix A_HD
        4. Compute eigenvalues and eigenvectors of A_HD
        5. Extract damping (real part) and frequency (imaginary part)

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
            - 'eigenvectors': Complex array of shape (n_hd, n_hd) where
              n_hd = (1 + 2×number_harmonics) × state_dimension
            - 'eigenvalues': Complex array of shape (n_hd,) containing
              eigenvalues λ = σ + jω where σ is damping, ω is frequency

        Raises:
        ------
        TypeError
            If OMEGA or T are not numeric
        ValueError
            If OMEGA <= 0, T <= 0, or problem configuration is invalid
        AttributeError
            If required problem attributes (number_harmonics, time_samples) are missing
        RuntimeError
            If HD matrix computation or eigenvalue extraction fails

        Notes:
        -----
        - Eigenvalues λ directly represent stability: σ < 0 → stable
        - Unlike LTP, no conversion from multipliers to exponents needed
        - HD matrix size: (1 + 2N) × n where N = number_harmonics
        - More harmonics → better accuracy but larger eigenvalue problem
        - Truncation error decreases with more harmonics
        - Requires time_samples for Fourier coefficient accuracy
        - Uses inherited hd_computer from StabilityAnalysis base class

        Problem Configuration:
        --------------------
        Requires rotor_build.problem to have:
        - number_harmonics: int >= 1, number of harmonics in Fourier expansion
        - time_samples: int >= 10, number of time points for sampling A(t)
        - hd_use_complex: bool, whether to use complex formulation

        Example:
        -------
        >>> hd = HDStability(rotor_build)
        >>> rotor_build.problem.number_harmonics = 3
        >>> rotor_build.problem.time_samples = 100
        >>> OMEGA = 25.0  # rad/s
        >>> T = 2 * np.pi / OMEGA
        >>> eigensol = hd.hd_single_point(OMEGA, T)
        >>> eigenvalues = eigensol['eigenvalues']
        >>> damping = np.real(eigenvalues)
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

        This method performs a sweep of operating points across the specified
        RPM range, computing eigenvalues at each point using the Harmonic
        Decomposition method. The results are stored in the modal_solution attribute.

        Process:
        -------
        1. Generate array of OMEGA values from problem definition using
           assign_range_OMEGA() method (inherited from base class)
        2. Compute corresponding periods T = 2π/OMEGA for each point
        3. For each operating point:
           - Construct HD matrix A_HD using Fourier expansion
           - Solve algebraic eigenvalue problem
           - Extract eigenvalues (stability indicators)
           - Store damping and frequency values

        Returns:
        -------
        HDStability
            Self with populated modal_solution list containing results at
            each operating point. Enables method chaining.

        Attributes Modified:
        -------------------
        modal_solution : list[ModalSolution]
            Populated with HD eigenvalue solutions at each RPM point.
            Length equals self.rotor_build.problem.number_points.

            Each entry contains:
            - OMEGA: Rotor speed (rad/s)
            - OMEGA_RPM: Rotor speed (RPM)
            - T: Period (s) = 2π/OMEGA
            - eigensolution: Dict with eigenvectors and eigenvalues of HD system
            - damping: Real part of eigenvalues (σ) - stability indicators
            - frequency: Imaginary part of eigenvalues (ω) - oscillation frequencies

        Raises:
        ------
        RuntimeError
            If HD matrix computation or eigenvalue extraction fails at
            any operating point
        AttributeError
            If problem configuration is missing required attributes

        Notes:
        -----
        - This method can be computationally expensive for many harmonics
        - HD matrix size: (1 + 2N) × n where N = harmonics, n = state dimension
        - Larger N → better accuracy but larger eigenvalue problems
        - No ODE integration required (advantage over LTP monodromy)
        - Eigenvalues directly represent stability (no conversion needed)
        - Method uses method chaining pattern for fluent interface
        - Can be faster than LTP for small N, slower for large N
        - HD has truncation error, LTP has integration error

        Example:
        -------
        >>> rotor = RotorBuild.build_all()
        >>> rotor.problem.required_solver = "HD"
        >>> rotor.problem.lower_rotor_RPM = 1000
        >>> rotor.problem.higher_rotor_RPM = 5000
        >>> rotor.problem.number_points = 20
        >>> rotor.problem.number_harmonics = 3
        >>> rotor.problem.time_samples = 100
        >>> hd = HDStability(rotor)
        >>> hd = hd.hd_full_range()
        >>> # Check stability at all points
        >>> max_damping = max([np.max(sol.damping) for sol in hd.modal_solution])
        >>> is_stable = max_damping < 0
        >>> print(f"System stable: {is_stable}")
        """
        # Assign OMEGA range and periods
        self.assign_range_OMEGA()
        self.assign_period_T()

        # Compute eigenvalues at each point
        for i in range(self.rotor_build.problem.number_points):
            try:
                eigensolution = self.hd_single_point(
                    self.modal_solution[i].OMEGA,
                    self.modal_solution[i].T
                )

                # Store results
                self.modal_solution[i].damping = np.real(eigensolution['eigenvalues'])
                self.modal_solution[i].frequency = np.imag(eigensolution['eigenvalues'])
                self.modal_solution[i].eigensolution = eigensolution

            except (ValueError, RuntimeError, TypeError, AttributeError) as e:
                raise RuntimeError(
                    f"HD analysis failed at point {i+1}/{self.rotor_build.problem.number_points}, "
                    f"OMEGA={self.modal_solution[i].OMEGA:.2f} rad/s "
                    f"({self.modal_solution[i].OMEGA_RPM:.2f} RPM), "
                    f"T={self.modal_solution[i].T:.6f}s: {e}"
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
