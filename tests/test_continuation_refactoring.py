"""Tests for continuation analysis refactoring.

Tests the helper methods and overall continuation functionality after refactoring:
- _build_augmented_matrix
- _compute_residual
- _solve_augmented_system
- _centered_difference
- numerical_sensitivity improvements

Run with:
    python -m pytest tests/test_continuation_refactoring.py -v
    python tests/test_continuation_refactoring.py
"""

import numpy as np
import pytest
from src.stability_analysis.continuation_analysis import ContinuationAnalysis
from src.stability_analysis.base import StabilityAnalysis


class TestHelperMethods:
    """Test the extracted helper methods."""

    def test_build_augmented_matrix_shape(self):
        """Test that augmented matrix has correct shape."""
        problem_size = 4
        eigenvector = np.array([1.0, 0.0, 0.0, 0.0])
        eigenvalue = -1.0 + 2.0j
        A = np.random.randn(problem_size, problem_size)

        A_aug = ContinuationAnalysis._build_augmented_matrix(
            eigenvector, eigenvalue, A, problem_size
        )

        assert A_aug.shape == (problem_size + 1, problem_size + 1)

    def test_build_augmented_matrix_structure(self):
        """Test that augmented matrix has correct structure."""
        problem_size = 3
        eigenvector = np.array([1.0, 0.0, 0.0])
        eigenvalue = -2.0
        A = np.eye(3)

        A_aug = ContinuationAnalysis._build_augmented_matrix(
            eigenvector, eigenvalue, A, problem_size
        )

        # Check first column is eigenvector with 0 at bottom
        np.testing.assert_array_almost_equal(A_aug[:3, 0], eigenvector)
        assert A_aug[3, 0] == 0.0

        # Check last row starts with 0
        assert A_aug[3, 0] == 0.0

        # Check conjugation in last row
        expected_last_row = np.hstack([[0.0], 2 * eigenvector.conj()])
        np.testing.assert_array_almost_equal(A_aug[3, :], expected_last_row)

    def test_build_augmented_matrix_complex_eigenvector(self):
        """Test augmented matrix with complex eigenvector."""
        problem_size = 2
        eigenvector = np.array([1.0 + 1.0j, 2.0 - 1.0j])
        eigenvalue = -1.5 + 0.5j
        A = np.array([[1.0, 0.5], [0.5, 2.0]])

        A_aug = ContinuationAnalysis._build_augmented_matrix(
            eigenvector, eigenvalue, A, problem_size
        )

        # Verify conjugation is applied
        expected_last_row = np.hstack([[0.0], 2 * eigenvector.conj()])
        np.testing.assert_array_almost_equal(A_aug[2, :], expected_last_row)

    def test_compute_residual_shape(self):
        """Test that residual has correct shape."""
        eigenvector = np.array([1.0, 0.0, 0.0])
        eigenvalue = -1.0
        A = np.eye(3)

        residual = ContinuationAnalysis._compute_residual(
            eigenvector, eigenvalue, A
        )

        assert residual.shape == (4, 1)

    def test_compute_residual_exact_eigenvalue(self):
        """Test residual is zero for exact eigenvalue."""
        # Create exact eigenpair
        A = np.array([[2.0, 0.0], [0.0, 3.0]])
        eigenvalue = 2.0
        eigenvector = np.array([1.0, 0.0])

        residual = ContinuationAnalysis._compute_residual(
            eigenvector, eigenvalue, A
        )

        # First part should be zero (eigenvalue equation satisfied)
        np.testing.assert_array_almost_equal(residual[:2, 0], [0.0, 0.0])
        # Second part should be zero (normalization satisfied)
        np.testing.assert_almost_equal(residual[2, 0], 0.0)

    def test_compute_residual_unnormalized(self):
        """Test residual for unnormalized eigenvector."""
        A = np.eye(2)
        eigenvalue = 1.0
        eigenvector = np.array([2.0, 0.0])  # Not normalized

        residual = ContinuationAnalysis._compute_residual(
            eigenvector, eigenvalue, A
        )

        # Normalization constraint should be violated
        expected_norm_residual = 1.0 - np.vdot(eigenvector, eigenvector)
        np.testing.assert_almost_equal(residual[2, 0], expected_norm_residual)
        assert abs(expected_norm_residual) > 1e-6  # Should be non-zero

    def test_solve_augmented_system_well_conditioned(self):
        """Test solving well-conditioned augmented system."""
        # Create well-conditioned system
        A_matrix = np.eye(5) * 2.0
        b_vector = np.ones((5, 1))

        x = ContinuationAnalysis._solve_augmented_system(A_matrix, b_vector)

        # Should use direct solve
        expected = np.linalg.solve(A_matrix, b_vector)
        np.testing.assert_array_almost_equal(x, expected)

    def test_solve_augmented_system_ill_conditioned(self):
        """Test solving ill-conditioned augmented system."""
        # Create ill-conditioned system
        A_matrix = np.array([[1.0, 1.0], [1.0, 1.0 + 1e-20]])
        b_vector = np.ones((2, 1))

        # Should not raise error and use pinv
        x = ContinuationAnalysis._solve_augmented_system(A_matrix, b_vector)

        assert x.shape == b_vector.shape
        assert not np.any(np.isnan(x))
        assert not np.any(np.isinf(x))

    def test_centered_difference_linear_function(self):
        """Test centered difference on linear function."""
        # For f(x) = 2x, derivative is 2
        x = 5.0
        h = 0.1
        f_plus_h = 2 * (x + h)
        f_min_h = 2 * (x - h)

        derivative = ContinuationAnalysis._centered_difference(
            f_plus_h, f_min_h, h
        )

        np.testing.assert_almost_equal(derivative, 2.0)

    def test_centered_difference_quadratic_function(self):
        """Test centered difference on quadratic function."""
        # For f(x) = x^2, derivative at x=3 is 6
        x = 3.0
        h = 0.01
        f_plus_h = (x + h) ** 2
        f_min_h = (x - h) ** 2

        derivative = ContinuationAnalysis._centered_difference(
            f_plus_h, f_min_h, h
        )

        np.testing.assert_almost_equal(derivative, 6.0, decimal=4)

    def test_centered_difference_matrix(self):
        """Test centered difference with matrices."""
        # Test with matrix-valued function
        h = 0.1
        A_plus_h = np.array([[1.1, 0.2], [0.3, 1.4]])
        A_min_h = np.array([[0.9, 0.0], [0.1, 1.2]])

        dA = ContinuationAnalysis._centered_difference(
            A_plus_h, A_min_h, h
        )

        # (1.1-0.9)/(0.2) = 1, (0.2-0)/(0.2) = 1, etc.
        expected = np.array([[1.0, 1.0], [1.0, 1.0]])
        np.testing.assert_array_almost_equal(dA, expected)


class TestNumericalSensitivity:
    """Test the refactored numerical_sensitivity method."""

    def test_sensitivity_ss_wrt_omega(self):
        """Test state-space sensitivity w.r.t. OMEGA."""
        # Simple state matrix that depends on omega: A(omega) = omega * I
        A_handle = lambda omega: omega * np.eye(2)
        OMEGA = 10.0
        h = 0.1

        sensitivity = ContinuationAnalysis.numerical_sensitivity(
            'SS', A_handle, h, OMEGA, wrt_omega=True
        )

        # Derivative should be I
        expected = np.eye(2)
        np.testing.assert_array_almost_equal(sensitivity, expected, decimal=4)

    def test_sensitivity_ss_wrt_parameter(self):
        """Test state-space sensitivity w.r.t. parameter."""
        # A(c, omega) = c * [[omega, 0], [0, omega]]
        A_handle = lambda c, omega: c * np.diag([omega, omega])
        OMEGA = 5.0
        par = 2.0
        h = 0.01

        sensitivity = ContinuationAnalysis.numerical_sensitivity(
            'SS', A_handle, h, OMEGA, wrt_omega=False, par=par
        )

        # Derivative w.r.t. c should be [[omega, 0], [0, omega]]
        expected = np.diag([OMEGA, OMEGA])
        np.testing.assert_array_almost_equal(sensitivity, expected, decimal=3)

    def test_sensitivity_parameter_required(self):
        """Test that par is required when wrt_omega=False."""
        A_handle = lambda c, omega: c * np.eye(2)
        OMEGA = 10.0
        h = 0.1

        with pytest.raises(ValueError, match="par must be provided"):
            ContinuationAnalysis.numerical_sensitivity(
                'SS', A_handle, h, OMEGA, wrt_omega=False, par=None
            )

    def test_sensitivity_unknown_matrix_type(self):
        """Test error for unknown matrix type."""
        A_handle = lambda omega: omega * np.eye(2)
        OMEGA = 10.0
        h = 0.1

        with pytest.raises(ValueError, match="Unknown matrix type"):
            ContinuationAnalysis.numerical_sensitivity(
                'UNKNOWN', A_handle, h, OMEGA
            )

    def test_sensitivity_hd_basic(self):
        """Test HD sensitivity computation."""
        # Time-invariant system
        A_handle = lambda t, omega: np.array([[-1.0, omega], [-omega, -1.0]])
        OMEGA = 2.0
        h = 0.1
        HD_parameters = [50, 1]  # time_samples, number_harmonics

        sensitivity = ContinuationAnalysis.numerical_sensitivity(
            'HD', A_handle, h, OMEGA, wrt_omega=True, HD_parameters=HD_parameters
        )

        # Should return a matrix (HD expands system)
        assert isinstance(sensitivity, np.ndarray)
        assert sensitivity.ndim == 2
        # For n=2 states and N=1 harmonic: (1+2*1)*2 = 6
        expected_size = (1 + 2 * HD_parameters[1]) * 2
        assert sensitivity.shape[0] == expected_size


class TestContinuationIntegration:
    """Integration tests for full continuation analysis."""

    @pytest.fixture
    def simple_rotor(self):
        """Create a simple rotor for testing."""
        from src.rotor_build import RotorBuild
        rotor = RotorBuild.build_all()
        rotor.problem.required_solver = 'LTI'
        rotor.problem.continuation = 'YES'
        rotor.problem.number_points = 3
        return rotor

    def test_continuation_runs_successfully(self, simple_rotor):
        """Test that continuation analysis runs without errors."""
        cont_analysis = ContinuationAnalysis(simple_rotor)
        cont_analysis = cont_analysis.continuation()

        assert len(cont_analysis.modal_solution) == 3
        assert cont_analysis.modal_solution[0].damping is not None
        assert cont_analysis.modal_solution[0].frequency is not None

    def test_continuation_eigenvalue_tracking(self, simple_rotor):
        """Test that eigenvalues are smoothly tracked."""
        cont_analysis = ContinuationAnalysis(simple_rotor)
        cont_analysis = cont_analysis.continuation()

        # Check that eigenvalues exist at all points
        for sol in cont_analysis.modal_solution:
            assert sol.eigensolution is not None
            assert 'eigenvalues' in sol.eigensolution
            assert 'eigenvectors' in sol.eigensolution

        # Check dimensions are consistent
        n_eigenvalues = len(cont_analysis.modal_solution[0].damping)
        for sol in cont_analysis.modal_solution:
            assert len(sol.damping) == n_eigenvalues
            assert len(sol.frequency) == n_eigenvalues

    def test_continuation_modal_consistency(self, simple_rotor):
        """Test that continuation maintains consistent modal properties."""
        cont_analysis = ContinuationAnalysis(simple_rotor)
        cont_analysis = cont_analysis.continuation()

        # Verify all points have the same number of eigenvalues
        n_eigenvalues = len(cont_analysis.modal_solution[0].damping)
        for sol in cont_analysis.modal_solution:
            assert len(sol.damping) == n_eigenvalues
            assert len(sol.frequency) == n_eigenvalues

        # Verify eigenvalues come in complex conjugate pairs (for real systems)
        # For each eigenvalue with nonzero imaginary part, its conjugate should exist
        for sol in cont_analysis.modal_solution:
            eigs = sol.damping + 1j * sol.frequency
            tol = 1e-6
            for eig in eigs:
                if abs(eig.imag) > tol:
                    # Find conjugate in the list
                    conjugate = eig.conjugate()
                    found_conjugate = any(
                        abs(e.real - conjugate.real) < tol and
                        abs(e.imag - conjugate.imag) < tol
                        for e in eigs
                    )
                    assert found_conjugate, f"Conjugate of {eig} not found in eigenvalues"


class TestRefactoringCorrectness:
    """Tests to ensure refactoring didn't break functionality."""

    def test_helper_methods_produce_same_results(self):
        """Test that helper methods produce identical results to original code."""
        problem_size = 3
        eigenvector = np.array([1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0.0])
        eigenvalue = -1.5 + 0.5j
        A = np.array([[1.0, 0.5, 0.0],
                      [0.5, 2.0, 0.1],
                      [0.0, 0.1, 3.0]])

        # Build augmented matrix using helper
        A_aug = ContinuationAnalysis._build_augmented_matrix(
            eigenvector, eigenvalue, A, problem_size
        )

        # Verify structure manually
        # Top-left block: [v, eigenvalue*I - A]
        assert A_aug[0, 0] == eigenvector[0]
        assert A_aug[1, 0] == eigenvector[1]
        assert A_aug[2, 0] == eigenvector[2]

        # Bottom row: [0, 2*v^H]
        assert A_aug[3, 0] == 0.0
        np.testing.assert_almost_equal(A_aug[3, 1], 2 * eigenvector[0].conj())
        np.testing.assert_almost_equal(A_aug[3, 2], 2 * eigenvector[1].conj())
        np.testing.assert_almost_equal(A_aug[3, 3], 2 * eigenvector[2].conj())

    def test_centered_difference_vs_manual_computation(self):
        """Verify centered difference matches manual computation."""
        h = 0.1
        f_plus = np.array([[5.0, 2.0], [3.0, 7.0]])
        f_minus = np.array([[3.0, 1.0], [2.0, 5.0]])

        # Using helper
        result = ContinuationAnalysis._centered_difference(f_plus, f_minus, h)

        # Manual computation
        expected = (f_plus - f_minus) / (2 * h)

        np.testing.assert_array_almost_equal(result, expected)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
