"""Modal participation analysis.

Converted from MATLAB modal_participation_analysis.m classdef.

This module implements modal participation analysis for rotor dynamics.
It computes how each harmonic component contributes to each mode, providing
insight into the frequency content of the system's dynamic behavior.

Two methods are supported:
- LTP (Linear Time-Periodic): Uses Floquet theory with ODE integration and FFT
- HD (Harmonic Balance): Extracts harmonics directly from expanded eigenvectors
"""

import numpy as np
from typing import List
from dataclasses import dataclass, field
from scipy.integrate import solve_ivp
from scipy.linalg import expm


@dataclass
class ModalParticipation:
    """Container for modal participation at a single operating point."""
    
    OMEGA: float = 0.0
    OMEGA_RPM: float = 0.0
    phi_jkn: np.ndarray = field(default_factory=lambda: np.array([]))


class ModalParticipationAnalysis:
    """Modal participation analysis class.
    
    This class computes modal participation factors showing how each
    harmonic contributes to each mode and state.
    
    Attributes:
        modal_participation: List of modal participation results
        modal_solution: Modal solutions from stability analysis
        rotor_build: Rotor build object
    """
    
    def __init__(self, modes):
        """Initialize modal participation analysis.
        
        Args:
            modes: StabilityAnalysis object with computed modes
        """
        self.modal_solution = modes.modal_solution
        self.rotor_build = modes.rotor_build
        self.modal_participation: List[ModalParticipation] = []
        
        # Run analysis
        self._modal_part_all()
    
    def _modal_part_all(self):
        """Compute modal participation for all operating points."""
        for i in range(len(self.modal_solution)):
            if self.rotor_build.problem.required_solver == "LTP":
                phi_jkn_mat = self._mod_part_single_ltp(i)
            elif self.rotor_build.problem.required_solver == "HD":
                phi_jkn_mat = self._mod_part_single_hd(i)
            else:
                # LTI doesn't have modal participation analysis
                continue

            # Store results
            mod_part = ModalParticipation()
            mod_part.OMEGA = self.modal_solution[i].OMEGA
            mod_part.OMEGA_RPM = self.modal_solution[i].OMEGA_RPM
            mod_part.phi_jkn = phi_jkn_mat
            self.modal_participation.append(mod_part)
    
    def _mod_part_single_ltp(self, index: int) -> np.ndarray:
        """Compute modal participation for single LTP point.

        Args:
            index: Index in modal_solution

        Returns:
            Modal participation matrix phi_jkn (j: state, k: mode, n: harmonic)
        """
        # Get the period T
        T = 2 * np.pi / self.modal_solution[index].OMEGA

        # Get state matrix function handle
        def A(t):
            return self.rotor_build.state_matrix_function(t, self.modal_solution[index].OMEGA)

        # Get size of A
        size_A = A(0).shape[0]

        # Retrieve characteristic exponents (complex eigenvalues)
        eta_M = (self.modal_solution[index].damping +
                 1j * self.modal_solution[index].frequency)

        # Create exponential matrix function
        def exp_mat(t):
            return expm(-np.diag(eta_M) * t)

        # Compute phi(t) for a set of time instants
        phi_0_matrix = np.eye(size_A)

        # Set ODE tolerances based on period
        if T > 1:
            rtol, atol = 1e-6, 1e-8
        else:
            rtol, atol = 1e-5, 1e-7

        time_phi = np.linspace(0, T, 500)
        phi_collection = np.zeros((size_A, size_A, len(time_phi)))

        # Integrate each column of phi_0_matrix
        for j in range(size_A):
            # Solve ODE for this initial condition
            sol = solve_ivp(
                lambda t, phi: self._odefun(t, phi, A),
                [0, T],
                phi_0_matrix[:, j],
                t_eval=time_phi,
                method='DOP853',
                rtol=rtol,
                atol=atol
            )

            # Store interpolated results
            for i in range(size_A):
                phi_collection[i, j, :] = sol.y[i, :]

        # Compute V(t) = phi(t) * eigenvectors * exp(-eta*t)
        V_collection = np.zeros((size_A, size_A, len(time_phi)), dtype=complex)
        eigenvectors = self.modal_solution[index].eigensolution['eigenvectors']

        for i in range(len(time_phi)):
            V_collection[:, :, i] = (phi_collection[:, :, i] @
                                     eigenvectors @
                                     exp_mat(time_phi[i]))

        # Compute Fourier series coefficients of V(t)
        # MATLAB uses number_harmonics*2 for FFT extraction
        N_harmonics = self.rotor_build.problem.number_harmonics * 2
        c_jk_matrix = np.zeros((size_A, size_A, 2*N_harmonics+1), dtype=complex)

        for j in range(size_A):  # j-th state
            for k in range(size_A):  # k-th mode
                c_jk = self._complex_coeff(V_collection[j, k, :], N_harmonics)
                c_jk_matrix[j, k, :] = c_jk

        # Normalization and computation of modal participation factors
        c_jk_matrix = np.abs(c_jk_matrix)

        phi_jkn_mat = np.zeros((size_A, size_A, 2*N_harmonics+1))

        for j in range(size_A):
            for k in range(size_A):
                # Normalize by sum across harmonics
                sum_harmonics = np.sum(c_jk_matrix[j, k, :])
                if sum_harmonics > 0:
                    phi_jkn_mat[j, k, :] = c_jk_matrix[j, k, :] / sum_harmonics

        return phi_jkn_mat
    
    def _mod_part_single_hd(self, index: int) -> np.ndarray:
        """Compute modal participation for single HD point.

        For Harmonic Balance method, eigenvectors are already expanded in harmonics.
        We extract the harmonic components and compute participation factors.

        Args:
            index: Index in modal_solution

        Returns:
            Modal participation matrix phi_jkn (j: state, k: mode, n: harmonic)
        """
        # Get eigenvector matrix (expanded with harmonics)
        eigenvector_matrix = self.modal_solution[index].eigensolution['eigenvectors']

        # Calculate block size (size of original state space)
        N_harmonics_full = 2 * self.rotor_build.problem.number_harmonics + 1
        block_size = eigenvector_matrix.shape[0] // N_harmonics_full

        n_modes = eigenvector_matrix.shape[1]

        # Initialize coefficient matrix
        c_jk_matrix = np.zeros((block_size, n_modes, N_harmonics_full), dtype=complex)

        # Process each eigenvector (mode)
        for k in range(n_modes):
            # Extract X_0 (constant component)
            X_0_k = eigenvector_matrix[0:block_size, k]

            # Extract harmonic components X_nc and X_ns
            X_nc_k = np.zeros((block_size, self.rotor_build.problem.number_harmonics), dtype=complex)
            X_ns_k = np.zeros((block_size, self.rotor_build.problem.number_harmonics), dtype=complex)

            j = block_size
            for n in range(self.rotor_build.problem.number_harmonics):
                X_nc_k[:, n] = eigenvector_matrix[j:j+block_size, k]
                X_ns_k[:, n] = eigenvector_matrix[j+block_size:j+2*block_size, k]
                j += 2 * block_size

            # Convert to complex exponential form
            # C_k^(+n) = (X_nc - i*X_ns) / 2  (positive harmonics)
            # C_k^(-n) = (X_nc + i*X_ns) / 2  (negative harmonics)
            C_k_nplus = 0.5 * (X_nc_k - 1j * X_ns_k)
            C_k_nminus = 0.5 * (X_nc_k + 1j * X_ns_k)

            # Flip for negative harmonics ordering
            C_k_nminus = np.fliplr(C_k_nminus)

            # Fill coefficient matrix
            # Negative harmonics (n = -N, ..., -1)
            for n in range(self.rotor_build.problem.number_harmonics):
                c_jk_matrix[:, k, n] = C_k_nminus[:, n]

            # Zero harmonic (n = 0)
            c_jk_matrix[:, k, self.rotor_build.problem.number_harmonics] = X_0_k

            # Positive harmonics (n = 1, ..., N)
            for n in range(self.rotor_build.problem.number_harmonics):
                c_jk_matrix[:, k, self.rotor_build.problem.number_harmonics + 1 + n] = C_k_nplus[:, n]

        # Normalization and computation of modal participation factors
        c_jk_matrix = np.abs(c_jk_matrix)

        phi_jkn_mat = np.zeros((block_size, n_modes, N_harmonics_full))

        for j in range(block_size):
            for k in range(n_modes):
                # Normalize by sum across harmonics
                sum_harmonics = np.sum(c_jk_matrix[j, k, :])
                if sum_harmonics > 0:
                    phi_jkn_mat[j, k, :] = c_jk_matrix[j, k, :] / sum_harmonics

        return phi_jkn_mat

    @staticmethod
    def _odefun(t: float, phi: np.ndarray, state_matrix_A) -> np.ndarray:
        """ODE function for integrating phi.

        Args:
            t: Time
            phi: State vector phi
            state_matrix_A: Function handle returning state matrix A(t)

        Returns:
            Time derivative dphi/dt = A(t) * phi
        """
        return state_matrix_A(t) @ phi

    @staticmethod
    def _complex_coeff(V_jk_t: np.ndarray, N_harmonics: int) -> np.ndarray:
        """Compute complex Fourier coefficients using FFT.

        This function uses FFT to determine the complex coefficients of the
        Fourier series for V(t).

        Args:
            V_jk_t: Time series data (complex array)
            N_harmonics: Number of harmonics to extract on each side

        Returns:
            Complex Fourier coefficients (2*N_harmonics+1 values)
        """
        # Compute FFT and center it
        all_coeff = np.fft.fftshift(np.fft.fft(V_jk_t) / len(V_jk_t))

        # Get middle index
        mid_index = len(V_jk_t) // 2

        # Extract required N_harmonics coefficients
        complex_coeffs = all_coeff[mid_index - N_harmonics:mid_index + N_harmonics + 1]

        return complex_coeffs


if __name__ == "__main__":
    # Test modal participation analysis
    from ..rotor_build import RotorBuild
    from ..stability_analysis.base import StabilityAnalysis
    
    print("Building rotor system...")
    rotor = RotorBuild.build_all()
    
    print("Running stability analysis...")
    modes = StabilityAnalysis.run_stability(rotor)
    
    print("Computing modal participation...")
    modal_part = ModalParticipationAnalysis(modes)
    
    print(f"Number of participation results: {len(modal_part.modal_participation)}")
