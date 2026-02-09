"""This is Bayer's sorting-free Hill's method, which uses Koopman embeddings to find the characteristic multipliers."""

import numpy as np
from scipy.linalg import expm
from typing import Dict
from .base import StabilityAnalysis


class SortingFreeStability(StabilityAnalysis):
    """Sorting-free Hill's method.

    Alternative to the integration of the monodromy matrix for time-periodic systems. Expands A(Ω,t) using Fourier series
    and solves an augmented algebraic eigenvalue problem (no ODE integration). Uses a sorting-free approach to find the characteristic multipliers from the
    eigenvalues of Hill's matrix (or HD state matrix), which are called candidates.

    This class contains more methods to perform this same task for comparison purposes:
    - choice of candidates through minimum imaginary part criterion.
    - choice of candidates through symmetry of eigenvectors.

    When applicable, it is best to compare the results with the direct integration of the monodromy matrix. 

    HD matrix size: (1 + 2N) × n where N = number_harmonics, n = state dimension.
    Eigenvalues eta = lambda + j*omega directly indicate stability (lambda < 0 stable).
    Characteristic multipliers Lambda = exp(diag(eta)*T).


    Returns from hill_single_point:
        {'A_HD': ndarray}

    Requires problem attributes: number_harmonics, time_samples, sorting_method, C_choice.
    Forces the use of complex HD through: use_hd_complex = True
    """

    def __init__(self, rotor_build) -> None:
        """Initialize sorting free stability analysis."""
        super().__init__(rotor_build)
    
    def hill_single_point(self, OMEGA: float) -> Dict[str, np.ndarray]:
        """COmpute Hill's matrix (or complex HD matrix) with the complex harmonic decomposition

        Args:
            OMEGA: Rotor speed (rad/s). Must be positive.

        Returns:
            A_HD: complex HD matrix
            number_states: the number of states in the original problem

        Raises:
            AttributeError: If required problem attributes missing.
            RuntimeError: If computation fails.
        """

        number_harmonics = self.rotor_build.problem.number_harmonics
        time_samples = self.rotor_build.problem.time_samples

        # Validate configuration values
        if not isinstance(number_harmonics, (int, np.integer)) or number_harmonics < 0:
            raise ValueError(
                f"number_harmonics must be integer >= 0, got {number_harmonics}"
            )

        if not isinstance(time_samples, (int, np.integer)) or time_samples < 2*number_harmonics+1:
            raise ValueError(
                f"time_samples must be integer >= {2*number_harmonics}, got {time_samples}"
            )

        # Compute period from OMEGA
        T = 2 * np.pi / OMEGA

        # Create time-dependent state matrix handle A_handle(OMEGA, t)
        try:
            A_handle = lambda Om, t: self.rotor_build.state_matrix_function(t, Om)
        except Exception as e:
            raise RuntimeError(
                f"Failed to create state matrix handle at OMEGA={OMEGA}: {e}"
            ) from e

        # Compute Harmonic Decomposition matrix (aka Hill's matrix)
        try:
            A_HD = StabilityAnalysis.hd_computer(
                A_handle,
                OMEGA,
                number_harmonics,
                time_samples,
                use_complex = True
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

        number_states = A_handle(1.0,0).shape[0]

        return A_HD,number_states
    
    
    def compute_W_matrix(self,size,number_harmonics) -> np.ndarray:
        """ Computes the constant W matrix according to the size of the problem.

        Args: 
                size: the size of the problem, i.e. of the HD matrix
                number_harmonics: the number of harmonics in the decomposition
        Returns: 
                W: vertical stack of identity matrices.
        """
        identity = np.eye(size)
        W = np.zeros((size*(2*number_harmonics+1),size))
        for i in range(1,2*number_harmonics+2):
            W[size*(i-1):size*i,:]=identity

        return W
    
    def compute_C_matrix(self,size,number_harmonics) -> np.ndarray:
        """ Computes the constant W matrix according to the size of the problem.

        Args: 
                size: the size of the problem, i.e. of the HD matrix
                number_harmonics: the number of harmonics in the decomposition
        Returns: 
                C: horizontal stack of identity matrices.
        """
        identity = np.eye(size)
        C = np.zeros((size,size*(2*number_harmonics+1)))
        C[:,size*number_harmonics:size*(number_harmonics+1)] = identity
        return C

    def char_mult_single_point(self,OMEGA: float) -> Dict[str,np.ndarray]:
        
        #should wrap the other functions
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
        
        number_harmonics = self.rotor_build.problem.number_harmonics

        A_HD,number_states = self.hill_single_point(OMEGA)

        W = self.compute_W_matrix(number_states,number_harmonics)
        
        C = self.compute_C_matrix(number_states,number_harmonics)
        
        expmat = expm(A_HD*T)

        monodromy = C@expmat@W

         # Compute eigenvalues and eigenvectors of expanded system
        try:
            eigenvalues, eigenvectors = np.linalg.eig(monodromy)
        except np.linalg.LinAlgError as e:
            raise RuntimeError(
                f"Eigenvalue computation failed for monodromy matrix at "
                f"OMEGA={OMEGA:.2f} rad/s: {e}"
            ) from e
        
        magnitude_squared = np.real(eigenvalues)**2 + np.imag(eigenvalues)**2

        # Avoid log of zero or negative values
        magnitude_squared = np.maximum(magnitude_squared, 1e-300)

        real_char_exp = (1 / (2 * T)) * np.log(magnitude_squared)

        # Imaginary part: ω = (1/T) * atan2(Im(μ), Re(μ))
        imag_char_exp = (1 / T) * np.arctan2(
            np.imag(eigenvalues),
            np.real(eigenvalues)
        )

        eigensolution = {
            'eigenvectors': eigenvectors,
            'char multipliers': eigenvalues,
            'real_char_exp': real_char_exp,
            'imag_char_exp': imag_char_exp

        }
        return eigensolution


    def char_mult_full_range(self) -> 'SortingFreeStability':
        """Compute characteristic multipliers over full RPM range using sorting-free method.

        Populates modal_solution with char multipliers and eigensolution
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
                eigensolution = self.char_mult_single_point(self.modal_solution[i].OMEGA)
                
                # Store results
                self.modal_solution[i].char_solution=eigensolution['char multipliers']
                self.modal_solution[i].eigensolution = eigensolution
                self.modal_solution[i].damping = eigensolution['real_char_exp']
                self.modal_solution[i].frequency = eigensolution['imag_char_exp']

            except (ValueError, RuntimeError, TypeError, AttributeError) as e:
                raise RuntimeError(
                    f"Sorting-free analysis failed at point {i+1}/{self.rotor_build.problem.number_points}, "
                    f"OMEGA={self.modal_solution[i].OMEGA:.2f} rad/s "
                    f"({self.modal_solution[i].OMEGA_RPM:.2f} RPM): {e}"
                ) from e

        return self


if __name__ == "__main__":
    # Test HD stability analysis
    from ..rotor_build import RotorBuild
    from ..config.problem_definition import create_default_problem
    from ..stability_analysis.ltp_stability import LTPStability
    
    
    print("Building rotor system...")
    problem = create_default_problem()

    problem.rotor_characteristics.damper_activation = "ODI"
    problem.sorting_method = "IMPART"
    problem.number_harmonics = 1
    problem.C_choice = "NAIVE"
    rotor = RotorBuild(problem)
    
    # Override to use desired solvers and sortings
    
    
    print("Running sorting-free stability analysis...")

    exact = LTPStability(rotor)
    exact = exact.ltp_full_range()

    sortfree = SortingFreeStability(rotor)
    sortfree = sortfree.char_mult_full_range()
    
    print(f"Number of solutions: {len(sortfree.modal_solution)}")
    print(f"First point OMEGA: {sortfree.modal_solution[0].OMEGA_RPM:.2f} RPM")
    if len(sortfree.modal_solution) > 0:
        print(f"Number of characteristic multipliers: {len(sortfree.modal_solution[0].damping)}")
        print(f"Sample characteristic multiplier at first point:")
        print(f"  Characteristic multipliers: {sortfree.modal_solution[0].char_solution}")
        print(f"Sample of exact characteristic multiplier at first point:")
        print(f"   Exact Characteristic multipliers: {exact.modal_solution[0].char_solution}")
