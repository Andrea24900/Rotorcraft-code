"""Base stability analysis class.

Converted from MATLAB stability_analysis.m classdef.
"""

import numpy as np
from typing import Callable, List, Dict, Any
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field


@dataclass
class ModalSolution:
    """Container for modal solution at a single operating point."""
    
    OMEGA: float = 0.0  # Rotor speed (rad/s)
    OMEGA_RPM: float = 0.0  # Rotor speed (RPM)
    T: float = 0.0  # Period (s)
    damping: np.ndarray = field(default_factory=lambda: np.array([]))
    frequency: np.ndarray = field(default_factory=lambda: np.array([]))
    eigensolution: Dict[str, Any] = field(default_factory=dict)
    char_solution: Dict[str, Any] = field(default_factory=dict)


class StabilityAnalysis:
    """Base class for stability analysis.
    
    This class contains main functions for various stability analyses
    including LTI, LTP, and HD methods.
    
    Attributes:
        rotor_build: Rotor build object with state matrices
        modal_solution: List of modal solutions at different operating points
    """
    
    def __init__(self, rotor_build):
        """Initialize stability analysis.
        
        Args:
            rotor_build: RotorBuild object
        """
        self.rotor_build = rotor_build
        self.modal_solution: List[ModalSolution] = []
    
    @staticmethod
    def run_stability(rotor_build):
        """Run full stability analysis according to problem definition.

        This is the main entry point that selects the appropriate analysis
        method based on the problem configuration.

        Args:
            rotor_build: RotorBuild object

        Returns:
            StabilityAnalysis subclass instance with computed results
        """
        from .lti_stability import LTIStability
        from .ltp_stability import LTPStability
        from .hd_stability import HDStability
        from .continuation_analysis import ContinuationAnalysis

        solver_type = rotor_build.problem.required_solver
        use_continuation = rotor_build.problem.continuation == "YES"

        if use_continuation:
            # Use continuation analysis for all solver types
            obj = ContinuationAnalysis(rotor_build)
            obj = obj.continuation()
        else:
            # Use standard analysis methods
            if solver_type == "LTI":
                obj = LTIStability(rotor_build)
                obj = obj.eigen_full_range()

            elif solver_type == "LTP":
                obj = LTPStability(rotor_build)
                obj = obj.LTP_full_range()

            elif solver_type == "HD":
                obj = HDStability(rotor_build)
                obj = obj.HD_full_range()
            else:
                raise ValueError(f"Unknown solver type: {solver_type}")

        return obj
    
    @staticmethod
    def monodromy_computer(
        A_handle_time: Callable[[float], np.ndarray],
        period_T: float
    ) -> np.ndarray:
        """Compute monodromy matrix for Linear Time-Periodic (LTP) systems.

        The monodromy matrix M characterizes the stability of time-periodic systems.
        It is defined as the fundamental solution matrix evaluated at one period:
        M = Φ(T), where dΦ/dt = A(t)Φ and Φ(0) = I.

        Algorithm:
        ---------
        For a system dx/dt = A(t)x with period T:
        1. Initialize fundamental matrix: Φ(0) = I (identity matrix)
        2. Integrate each column j: dΦⱼ/dt = A(t)Φⱼ from t=0 to t=T
        3. Assemble columns into monodromy matrix: M = [Φ₁(T), Φ₂(T), ..., Φₙ(T)]
        4. Eigenvalues of M are characteristic multipliers μ
        5. Stability: |μ| < 1 → stable, |μ| > 1 → unstable

        Integration Method:
        ------------------
        Uses scipy's DOP853 (Dormand-Prince 8th order) integrator with:
        - Relative tolerance: 1e-6
        - Absolute tolerance: 1e-8
        - Suitable for long-period integrations with high accuracy

        Args:
        ----
        A_handle_time : Callable[[float], np.ndarray]
            Function handle A(t) that returns the n×n state matrix at time t
        period_T : float
            Period of oscillation (s), typically T = 2π/OMEGA for rotor systems

        Returns:
        -------
        np.ndarray
            n×n monodromy matrix M whose eigenvalues are characteristic multipliers

        Notes:
        -----
        - Computational cost: O(n) integrations, each solving n ODEs
        - Total complexity: O(n²) for each monodromy computation
        - For large systems or many operating points, this can be expensive
        - The method is numerically stable for well-conditioned systems
        - Fixed from previous implementation that incorrectly reshaped vectors

        Example:
        -------
        >>> A = lambda t: np.array([[0, 1], [-1 - 0.1*np.sin(t), -0.1]])
        >>> T = 2 * np.pi  # Period of time-varying term
        >>> M = StabilityAnalysis.monodromy_computer(A, T)
        >>> multipliers = np.linalg.eigvals(M)
        >>> is_stable = np.all(np.abs(multipliers) < 1)
        """
        n = A_handle_time(0).shape[0]
        phi_0_matrix = np.eye(n)
        monodromy_matrix_M = np.zeros((n, n))

        def odefun_monodromy(t, phi_vec, A_func):
            """ODE function for fundamental solution matrix."""
            # phi_vec is a single column vector
            dphi_dt = A_func(t) @ phi_vec
            return dphi_dt

        # Integrate each column of the fundamental solution matrix
        for j in range(n):
            phi_0 = phi_0_matrix[:, j]

            sol = solve_ivp(
                lambda t, phi: odefun_monodromy(t, phi, A_handle_time),
                [0, period_T],
                phi_0,
                method='DOP853',
                rtol=1e-6,
                atol=1e-8
            )

            # Extract final value
            phi_T = sol.y[:, -1]
            monodromy_matrix_M[:, j] = phi_T

        return monodromy_matrix_M
    
    @staticmethod
    def HD_computer(
        A_handle_time: Callable[[float], np.ndarray],
        time: np.ndarray,
        number_harmonics: int,
        OMEGA: float
    ) -> np.ndarray:
        """Compute the state matrix for HD analysis.

        This method constructs the expanded HD matrix by computing
        Fourier coefficients of the time-periodic state matrix A(t).

        Args:
            A_handle_time: Function handle A(t)
            time: Time vector for sampling
            number_harmonics: Number of harmonics to include
            OMEGA: Fundamental frequency

        Returns:
            matrix A_HD
        """
        state_number = A_handle_time(0).shape[0]
        size_A_HD = (1 + 2 * number_harmonics) * state_number
        A_HD = np.zeros((size_A_HD, size_A_HD))

        # Compute time realizations of A(t)
        time_realisation_A = np.zeros((state_number, state_number, len(time)))
        for i, t in enumerate(time):
            time_realisation_A[:, :, i] = A_handle_time(t)

        # Compute FFT coefficients for each element
        # Storage: A_coeff[:,:,0] = A_0, A_coeff[:,:,2*i] = A_ic, A_coeff[:,:,2*i+1] = A_is
        A_coeff = np.zeros((state_number, state_number, 2 * number_harmonics + 1))

        for j in range(state_number):
            for k in range(state_number):
                # FFT of A[j,k](t)
                full_transform = np.fft.fft(time_realisation_A[j, k, :]) / len(time)

                # Extract coefficients: A_0, A_1c, A_1s, A_2c, A_2s, ...
                A_coeff[j, k, 0] = full_transform[0].real  # A_0

                index_coeff = 1
                for n in range(1, number_harmonics + 1):
                    A_coeff[j, k, index_coeff] = 2 * full_transform[n].real  # A_nc
                    A_coeff[j, k, index_coeff + 1] = -2 * full_transform[n].imag  # A_ns
                    index_coeff += 2

        # Assignment of first diagonal block (H0M)
        A_HD[0:state_number, 0:state_number] = StabilityAnalysis._H0M_assign(A_coeff)

        # Assignment of first row
        col_A_HD = state_number
        i = 1
        while col_A_HD < size_A_HD:
            A_HD[0:state_number, col_A_HD:col_A_HD + state_number] = \
                StabilityAnalysis._H0MiC_assign(A_coeff, i)
            A_HD[0:state_number, col_A_HD + state_number:col_A_HD + 2 * state_number] = \
                StabilityAnalysis._H0MiS_assign(A_coeff, i)
            i += 1
            col_A_HD += 2 * state_number

        # Assignment of first column
        row_A_HD = state_number
        i = 1
        while row_A_HD < size_A_HD:
            A_HD[row_A_HD:row_A_HD + state_number, 0:state_number] = \
                StabilityAnalysis._HiCM_assign(A_coeff, i)
            A_HD[row_A_HD + state_number:row_A_HD + 2 * state_number, 0:state_number] = \
                StabilityAnalysis._HiSM_assign(A_coeff, i)
            i += 1
            row_A_HD += 2 * state_number

        # Assignment of diagonal terms
        col_A_HD = state_number
        row_A_HD = state_number
        i = 1

        while row_A_HD < size_A_HD:
            # Upper left corner element: HiCMiC
            A_HD[row_A_HD:row_A_HD + state_number, col_A_HD:col_A_HD + state_number] = \
                StabilityAnalysis._HiCMjC_assign(A_coeff, i, i)

            # Upper right corner element: HiCMiS - i*OMEGA
            A_HD[row_A_HD:row_A_HD + state_number, col_A_HD + state_number:col_A_HD + 2 * state_number] = \
                StabilityAnalysis._HiCMjS_assign(A_coeff, i, i) - np.eye(state_number) * i * OMEGA

            # Lower left corner element: HiSMiC + i*OMEGA
            A_HD[row_A_HD + state_number:row_A_HD + 2 * state_number, col_A_HD:col_A_HD + state_number] = \
                StabilityAnalysis._HiSMjC_assign(A_coeff, i, i) + np.eye(state_number) * i * OMEGA

            # Lower right corner element: HiSMiS
            A_HD[row_A_HD + state_number:row_A_HD + 2 * state_number,
                 col_A_HD + state_number:col_A_HD + 2 * state_number] = \
                StabilityAnalysis._HiSMjS_assign(A_coeff, i, i)

            i += 1
            row_A_HD += 2 * state_number
            col_A_HD += 2 * state_number

        # Assignment of extra-diagonal terms (matrix scanned by columns, then by rows)
        col_A_HD = 3 * state_number  # Skip first diagonal block
        row_A_HD = state_number
        i = 1  # Index for "block rows"
        j = 2  # Index for "block columns"

        while row_A_HD < size_A_HD:
            while col_A_HD < size_A_HD:
                if i != j:
                    # Upper left corner element: HiCMjC
                    A_HD[row_A_HD:row_A_HD + state_number, col_A_HD:col_A_HD + state_number] = \
                        StabilityAnalysis._HiCMjC_assign(A_coeff, i, j)

                    # Upper right corner element: HiCMjS
                    A_HD[row_A_HD:row_A_HD + state_number,
                         col_A_HD + state_number:col_A_HD + 2 * state_number] = \
                        StabilityAnalysis._HiCMjS_assign(A_coeff, i, j)

                    # Lower left corner element: HiSMjC
                    A_HD[row_A_HD + state_number:row_A_HD + 2 * state_number,
                         col_A_HD:col_A_HD + state_number] = \
                        StabilityAnalysis._HiSMjC_assign(A_coeff, i, j)

                    # Lower right corner element: HiSMjS
                    A_HD[row_A_HD + state_number:row_A_HD + 2 * state_number,
                         col_A_HD + state_number:col_A_HD + 2 * state_number] = \
                        StabilityAnalysis._HiSMjS_assign(A_coeff, i, j)

                j += 1
                col_A_HD += 2 * state_number

            # Restart from first block column in the next block row
            col_A_HD = state_number
            j = 1  # Restart from first set of coefficients in the next row
            i += 1
            row_A_HD += 2 * state_number

        return A_HD

    @staticmethod
    def HD_computer_complex(
        A_handle_time: Callable[[float], np.ndarray],
        time: np.ndarray,
        number_harmonics: int,
        OMEGA: float
    ) -> np.ndarray:
        """Compute HD matrix using complex formulation.

        This is an alternative formulation using complex exponentials
        instead of separate sine/cosine terms.

        Args:
            A_handle_time: Function handle A(t)
            time: Time vector for sampling
            number_harmonics: Number of harmonics to include
            OMEGA: Fundamental frequency

        Returns:
            Complex HD matrix A_HD
        """
        state_number = A_handle_time(0).shape[0]
        size_A_HD = (1 + 2 * number_harmonics) * state_number
        A_HD = np.zeros((size_A_HD, size_A_HD), dtype=complex)

        # Compute time realizations of A(t)
        time_realisation_A = np.zeros((state_number, state_number, len(time)))
        for i, t in enumerate(time):
            time_realisation_A[:, :, i] = A_handle_time(t)

        # Compute FFT coefficients for each element
        A_coeff = np.zeros((state_number, state_number, len(time)), dtype=complex)

        for j in range(state_number):
            for k in range(state_number):
                A_coeff[j, k, :] = np.fft.fft(time_realisation_A[j, k, :]) / len(time)

        # k is used to cycle the harmonics from -N to N (columns)
        # j is used to cycle the blocks in the A matrix (rows)
        k = -number_harmonics
        j = -number_harmonics

        for block_row in range(2 * number_harmonics + 1):
            for block_column in range(2 * number_harmonics + 1):
                m = j - k

                if m < 0:
                    A_m = A_coeff[:, :, m]  # Negative index wraps around in Python
                elif m == 0:
                    A_m = A_coeff[:, :, 0] - 1j * k * OMEGA * np.eye(state_number)
                else:  # m > 0
                    A_m = A_coeff[:, :, m]

                A_HD[block_row * state_number:(block_row + 1) * state_number,
                     block_column * state_number:(block_column + 1) * state_number] = A_m

                # Cycle the column index
                k += 1

            # Cycle the row index
            j += 1
            # Reset the column index
            k = -number_harmonics

        return A_HD

    # H-matrix assignment helper functions
    @staticmethod
    def _H0M_assign(M_coeff: np.ndarray) -> np.ndarray:
        """Assign H_0_M term (constant part)."""
        return M_coeff[:, :, 0]

    @staticmethod
    def _H0MiC_assign(M_coeff: np.ndarray, i: int) -> np.ndarray:
        """Assign H_0_M_iC term (cosine coefficient for first row).

        Args:
            M_coeff: Coefficient array where [:,:,0]=A0, [:,:,1]=A1c, [:,:,2]=A1s, etc.
            i: Harmonic number (1-based: 1, 2, 3...)
        """
        return M_coeff[:, :, 2 * i - 1] / 2

    @staticmethod
    def _H0MiS_assign(M_coeff: np.ndarray, i: int) -> np.ndarray:
        """Assign H_0_M_iS term (sine coefficient for first row).

        Args:
            M_coeff: Coefficient array where [:,:,0]=A0, [:,:,1]=A1c, [:,:,2]=A1s, etc.
            i: Harmonic number (1-based: 1, 2, 3...)
        """
        return M_coeff[:, :, 2 * i] / 2

    @staticmethod
    def _HiCM_assign(M_coeff: np.ndarray, i: int) -> np.ndarray:
        """Assign H_iC_M term (cosine coefficient for first column).

        Args:
            M_coeff: Coefficient array where [:,:,0]=A0, [:,:,1]=A1c, [:,:,2]=A1s, etc.
            i: Harmonic number (1-based: 1, 2, 3...)
        """
        return M_coeff[:, :, 2 * i - 1]

    @staticmethod
    def _HiSM_assign(M_coeff: np.ndarray, i: int) -> np.ndarray:
        """Assign H_iS_M term (sine coefficient for first column).

        Args:
            M_coeff: Coefficient array where [:,:,0]=A0, [:,:,1]=A1c, [:,:,2]=A1s, etc.
            i: Harmonic number (1-based: 1, 2, 3...)
        """
        return M_coeff[:, :, 2 * i]

    @staticmethod
    def _HiCMjC_assign(M_coeff: np.ndarray, i: int, j: int) -> np.ndarray:
        """Assign H_iC_M_jC term (cosine-cosine interaction).

        Args:
            M_coeff: Coefficient array where [:,:,0]=A0, [:,:,1]=A1c, [:,:,2]=A1s, etc.
            i, j: Harmonic numbers (1-based)
        """
        if i == j:
            k = i + j
            # Safe access: check if k harmonic exists
            if 2 * k <= M_coeff.shape[2]:
                return M_coeff[:, :, 0] + M_coeff[:, :, 2 * k - 1] / 2
            else:
                return M_coeff[:, :, 0]
        elif i > j:
            k = i + j
            l = i - j
            result = M_coeff[:, :, 0] * 0  # Initialize with correct shape
            if 2 * l - 1 < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * l - 1] / 2
            if 2 * k - 1 < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * k - 1] / 2
            return result
        else:  # j > i
            k = i + j
            m = j - i
            result = M_coeff[:, :, 0] * 0
            if 2 * m - 1 < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * m - 1] / 2
            if 2 * k - 1 < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * k - 1] / 2
            return result

    @staticmethod
    def _HiCMjS_assign(M_coeff: np.ndarray, i: int, j: int) -> np.ndarray:
        """Assign H_iC_M_jS term (cosine-sine interaction).

        Args:
            M_coeff: Coefficient array where [:,:,0]=A0, [:,:,1]=A1c, [:,:,2]=A1s, etc.
            i, j: Harmonic numbers (1-based)
        """
        k = i + j
        result = M_coeff[:, :, 0] * 0

        if i == j:
            if 2 * k < M_coeff.shape[2]:
                return M_coeff[:, :, 2 * k] / 2
            return result
        elif i > j:
            l = i - j
            if 2 * k < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * k] / 2
            if 2 * l < M_coeff.shape[2]:
                result -= M_coeff[:, :, 2 * l] / 2
            return result
        else:  # j > i
            m = j - i
            if 2 * m < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * m] / 2
            if 2 * k < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * k] / 2
            return result

    @staticmethod
    def _HiSMjC_assign(M_coeff: np.ndarray, i: int, j: int) -> np.ndarray:
        """Assign H_iS_M_jC term (sine-cosine interaction).

        Args:
            M_coeff: Coefficient array where [:,:,0]=A0, [:,:,1]=A1c, [:,:,2]=A1s, etc.
            i, j: Harmonic numbers (1-based)
        """
        k = i + j
        result = M_coeff[:, :, 0] * 0

        if i == j:
            if 2 * k < M_coeff.shape[2]:
                return M_coeff[:, :, 2 * k] / 2
            return result
        elif i > j:
            l = i - j
            if 2 * k < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * k] / 2
            if 2 * l < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * l] / 2
            return result
        else:  # j > i
            m = j - i
            if 2 * k < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * k] / 2
            if 2 * m < M_coeff.shape[2]:
                result -= M_coeff[:, :, 2 * m] / 2
            return result

    @staticmethod
    def _HiSMjS_assign(M_coeff: np.ndarray, i: int, j: int) -> np.ndarray:
        """Assign H_iS_M_jS term (sine-sine interaction).

        Args:
            M_coeff: Coefficient array where [:,:,0]=A0, [:,:,1]=A1c, [:,:,2]=A1s, etc.
            i, j: Harmonic numbers (1-based)
        """
        if i == j:
            k = i + j
            if 2 * k - 1 < M_coeff.shape[2]:
                return M_coeff[:, :, 0] - M_coeff[:, :, 2 * k - 1] / 2
            else:
                return M_coeff[:, :, 0]
        elif i > j:
            k = i + j
            l = i - j
            result = M_coeff[:, :, 0] * 0
            if 2 * l - 1 < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * l - 1] / 2
            if 2 * k - 1 < M_coeff.shape[2]:
                result -= M_coeff[:, :, 2 * k - 1] / 2
            return result
        else:  # j > i
            k = i + j
            m = j - i
            result = M_coeff[:, :, 0] * 0
            if 2 * m - 1 < M_coeff.shape[2]:
                result += M_coeff[:, :, 2 * m - 1] / 2
            if 2 * k - 1 < M_coeff.shape[2]:
                result -= M_coeff[:, :, 2 * k - 1] / 2
            return result

    def assign_range_OMEGA(self) -> 'StabilityAnalysis':
        """Assign range of rotor speeds to analyze.
        
        Returns:
            Self (for method chaining)
        """
        rotor_OMEGA_vector = np.linspace(
            self.rotor_build.problem.lower_rotor_RPM / 60 * 2 * np.pi,
            self.rotor_build.problem.higher_rotor_RPM / 60 * 2 * np.pi,
            self.rotor_build.problem.number_points
        )
        
        rotor_OMEGA_vector_RPM = np.linspace(
            self.rotor_build.problem.lower_rotor_RPM,
            self.rotor_build.problem.higher_rotor_RPM,
            self.rotor_build.problem.number_points
        )
        
        # Reverse if solving right to left
        if self.rotor_build.problem.solution_direction == "R2L":
            rotor_OMEGA_vector = rotor_OMEGA_vector[::-1]
            rotor_OMEGA_vector_RPM = rotor_OMEGA_vector_RPM[::-1]
        
        # Initialize modal solution list
        self.modal_solution = []
        for i in range(self.rotor_build.problem.number_points):
            modal_sol = ModalSolution()
            modal_sol.OMEGA = rotor_OMEGA_vector[i]
            modal_sol.OMEGA_RPM = rotor_OMEGA_vector_RPM[i]
            self.modal_solution.append(modal_sol)
        
        return self
    
    def assign_period_T(self) -> 'StabilityAnalysis':
        """Assign period T = 2π/Ω for each operating point.
        
        Returns:
            Self (for method chaining)
        """
        for modal_sol in self.modal_solution:
            modal_sol.T = 2 * np.pi / modal_sol.OMEGA
        
        return self


if __name__ == "__main__":
    # Test stability analysis base class
    from ..rotor_build import RotorBuild
    
    print("Building rotor system...")
    rotor = RotorBuild.build_all()
    
    print("Initializing stability analysis...")
    stability = StabilityAnalysis(rotor)
    
    print("Assigning OMEGA range...")
    stability = stability.assign_range_OMEGA()
    print(f"Number of points: {len(stability.modal_solution)}")
    print(f"First OMEGA: {stability.modal_solution[0].OMEGA:.2f} rad/s")
    print(f"Last OMEGA: {stability.modal_solution[-1].OMEGA:.2f} rad/s")
