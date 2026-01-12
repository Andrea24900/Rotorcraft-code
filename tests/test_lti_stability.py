"""Unit tests for stability_analysis.lti_stability module.

Tests the LTIStability class and its methods including:
- eigen_single_point (single operating point analysis)
- eigen_full_range (full RPM range analysis)
- Input validation
- Error handling
- Integration with real rotor systems

Run with:
    python -m pytest tests/test_lti_stability.py -v
    python tests/test_lti_stability.py
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pytest
from src.stability_analysis.lti_stability import LTIStability
from src.rotor_build import RotorBuild


# ============================================================================
# Test LTIStability Initialization
# ============================================================================

def test_lti_initialization():
    """Test LTIStability initialization."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    assert lti.rotor_build is rotor
    assert isinstance(lti.modal_solution, list)
    assert len(lti.modal_solution) == 0


def test_lti_inherits_from_base():
    """Test that LTIStability properly inherits from StabilityAnalysis."""
    from src.stability_analysis.base import StabilityAnalysis

    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    assert isinstance(lti, StabilityAnalysis)
    # Should have inherited methods
    assert hasattr(lti, 'assign_range_OMEGA')
    assert hasattr(lti, 'assign_period_T')


# ============================================================================
# Test eigen_single_point - Basic Functionality
# ============================================================================

def test_eigen_single_point_basic():
    """Test eigen_single_point with valid rotor system."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    OMEGA = 25.0  # rad/s
    eigensol = lti.eigen_single_point(OMEGA)

    # Check return type
    assert isinstance(eigensol, dict)
    assert 'eigenvalues' in eigensol
    assert 'eigenvectors' in eigensol

    # Check array types
    assert isinstance(eigensol['eigenvalues'], np.ndarray)
    assert isinstance(eigensol['eigenvectors'], np.ndarray)

    # Check dimensions
    n_states = 2 * rotor.mass_matrix.shape[0]  # State space dimension
    assert eigensol['eigenvalues'].shape == (n_states,)
    assert eigensol['eigenvectors'].shape == (n_states, n_states)


def test_eigen_single_point_eigenvalues_complex():
    """Test that eigenvalues can be complex."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    eigensol = lti.eigen_single_point(50.0)
    eigenvalues = eigensol['eigenvalues']

    # Oscillatory systems should have complex eigenvalues
    has_complex = np.any(np.imag(eigenvalues) != 0)
    assert has_complex, "Rotor system should have complex eigenvalues"


def test_eigen_single_point_eigenvectors_shape():
    """Test that eigenvectors have correct shape."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    eigensol = lti.eigen_single_point(100.0)
    eigenvectors = eigensol['eigenvectors']

    # Should be square matrix
    assert eigenvectors.shape[0] == eigenvectors.shape[1]


def test_eigen_single_point_different_speeds():
    """Test eigen_single_point at different rotor speeds."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    speeds = [10.0, 50.0, 100.0, 200.0]  # rad/s
    eigenvalue_sets = []

    for omega in speeds:
        eigensol = lti.eigen_single_point(omega)
        eigenvalue_sets.append(eigensol['eigenvalues'])

    # Eigenvalues should be different at different speeds
    for i in range(len(speeds) - 1):
        difference = np.linalg.norm(eigenvalue_sets[i] - eigenvalue_sets[i+1])
        assert difference > 1e-6, f"Eigenvalues should differ between {speeds[i]} and {speeds[i+1]} rad/s"


# ============================================================================
# Test eigen_single_point - Input Validation
# ============================================================================

def test_eigen_single_point_negative_omega():
    """Test that negative OMEGA raises ValueError."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    with pytest.raises(ValueError, match="OMEGA must be positive"):
        lti.eigen_single_point(-10.0)


def test_eigen_single_point_zero_omega():
    """Test that zero OMEGA raises ValueError."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    with pytest.raises(ValueError, match="OMEGA must be positive"):
        lti.eigen_single_point(0.0)


def test_eigen_single_point_non_numeric_omega():
    """Test that non-numeric OMEGA raises TypeError."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    with pytest.raises(TypeError, match="OMEGA must be numeric"):
        lti.eigen_single_point("not_a_number")

    with pytest.raises(TypeError, match="OMEGA must be numeric"):
        lti.eigen_single_point(None)

    with pytest.raises(TypeError, match="OMEGA must be numeric"):
        lti.eigen_single_point([25.0])


def test_eigen_single_point_numpy_scalar():
    """Test that numpy scalar types are accepted."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    # Should not raise exceptions
    eigensol1 = lti.eigen_single_point(np.float64(25.0))
    eigensol2 = lti.eigen_single_point(np.int32(25))

    assert isinstance(eigensol1, dict)
    assert isinstance(eigensol2, dict)


# ============================================================================
# Test eigen_full_range - Basic Functionality
# ============================================================================

def test_eigen_full_range_basic():
    """Test eigen_full_range with default rotor configuration."""
    rotor = RotorBuild.build_all()
    rotor.problem.required_solver = "LTI"

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    # Check that modal_solution is populated
    assert len(lti.modal_solution) == rotor.problem.number_points
    assert len(lti.modal_solution) > 0


def test_eigen_full_range_modal_solution_structure():
    """Test that modal_solution has correct structure."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 100.0
    rotor.problem.higher_rotor_RPM = 200.0
    rotor.problem.number_points = 5

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    # Check each modal solution
    for i, modal_sol in enumerate(lti.modal_solution):
        # Check required attributes
        assert hasattr(modal_sol, 'OMEGA')
        assert hasattr(modal_sol, 'OMEGA_RPM')
        assert hasattr(modal_sol, 'damping')
        assert hasattr(modal_sol, 'frequency')
        assert hasattr(modal_sol, 'eigensolution')

        # Check OMEGA is positive
        assert modal_sol.OMEGA > 0

        # Check arrays are populated
        assert len(modal_sol.damping) > 0
        assert len(modal_sol.frequency) > 0
        assert len(modal_sol.damping) == len(modal_sol.frequency)

        # Check eigensolution structure
        assert 'eigenvalues' in modal_sol.eigensolution
        assert 'eigenvectors' in modal_sol.eigensolution


def test_eigen_full_range_rpm_values():
    """Test that RPM values are correctly set."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 1000.0
    rotor.problem.higher_rotor_RPM = 2000.0
    rotor.problem.number_points = 11

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    # Extract RPM values
    rpm_values = [sol.OMEGA_RPM for sol in lti.modal_solution]

    # Check range (default is R2L, so descending)
    if rotor.problem.solution_direction == "R2L":
        assert np.isclose(rpm_values[0], 2000.0, atol=1.0)
        assert np.isclose(rpm_values[-1], 1000.0, atol=1.0)
    else:  # L2R
        assert np.isclose(rpm_values[0], 1000.0, atol=1.0)
        assert np.isclose(rpm_values[-1], 2000.0, atol=1.0)


def test_eigen_full_range_damping_frequency_relationship():
    """Test that damping and frequency are real and imaginary parts of eigenvalues."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 3

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    for modal_sol in lti.modal_solution:
        eigenvalues = modal_sol.eigensolution['eigenvalues']
        damping = modal_sol.damping
        frequency = modal_sol.frequency

        # Check that damping = real part
        assert np.allclose(damping, np.real(eigenvalues))

        # Check that frequency = imaginary part
        assert np.allclose(frequency, np.imag(eigenvalues))


def test_eigen_full_range_method_chaining():
    """Test that eigen_full_range returns self for method chaining."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 3

    lti = LTIStability(rotor)
    result = lti.eigen_full_range()

    assert result is lti, "Method should return self for chaining"


# ============================================================================
# Test eigen_full_range - Different Configurations
# ============================================================================

def test_eigen_full_range_single_point():
    """Test eigen_full_range with single operating point."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 1000.0
    rotor.problem.higher_rotor_RPM = 1000.0
    rotor.problem.number_points = 1

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    assert len(lti.modal_solution) == 1
    assert np.isclose(lti.modal_solution[0].OMEGA_RPM, 1000.0)


def test_eigen_full_range_many_points():
    """Test eigen_full_range with many operating points."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 100.0  # Changed from 0.0 to avoid OMEGA=0
    rotor.problem.higher_rotor_RPM = 5000.0
    rotor.problem.number_points = 50

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    assert len(lti.modal_solution) == 50

    # Check that all points were computed
    for modal_sol in lti.modal_solution:
        assert len(modal_sol.damping) > 0
        assert len(modal_sol.frequency) > 0


def test_eigen_full_range_l2r_direction():
    """Test eigen_full_range with L2R solution direction."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 100.0
    rotor.problem.higher_rotor_RPM = 200.0
    rotor.problem.number_points = 5
    rotor.problem.solution_direction = "L2R"

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    # Check ascending order
    rpm_values = [sol.OMEGA_RPM for sol in lti.modal_solution]
    assert np.isclose(rpm_values[0], 100.0, atol=1.0)
    assert np.isclose(rpm_values[-1], 200.0, atol=1.0)


def test_eigen_full_range_r2l_direction():
    """Test eigen_full_range with R2L solution direction."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 100.0
    rotor.problem.higher_rotor_RPM = 200.0
    rotor.problem.number_points = 5
    rotor.problem.solution_direction = "R2L"

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    # Check descending order
    rpm_values = [sol.OMEGA_RPM for sol in lti.modal_solution]
    assert np.isclose(rpm_values[0], 200.0, atol=1.0)
    assert np.isclose(rpm_values[-1], 100.0, atol=1.0)


# ============================================================================
# Test Stability Analysis Properties
# ============================================================================

def test_eigenvalue_properties():
    """Test that eigenvalues have expected properties."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 3

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    for modal_sol in lti.modal_solution:
        eigenvalues = modal_sol.eigensolution['eigenvalues']

        # Check that complex conjugate pairs exist for oscillatory modes
        # (this is a property of real state matrices)
        for i, lam in enumerate(eigenvalues):
            if np.abs(np.imag(lam)) > 1e-10:  # If significantly complex
                # Find conjugate pair
                conjugate = np.conj(lam)
                found_conjugate = np.any(np.abs(eigenvalues - conjugate) < 1e-8)
                assert found_conjugate, f"Complex eigenvalue {lam} should have conjugate pair"


def test_damping_interpretation():
    """Test interpretation of damping values."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 3

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    # For each operating point, check damping interpretation
    for modal_sol in lti.modal_solution:
        damping = modal_sol.damping

        # All values should be real (since they're real parts of eigenvalues)
        assert np.all(np.isreal(damping))

        # Damping < 0 indicates stable mode
        # Damping > 0 indicates unstable mode
        stable_modes = np.sum(damping < 0)
        unstable_modes = np.sum(damping > 0)

        # At least some modes should exist
        assert stable_modes + unstable_modes == len(damping)


def test_natural_frequency_computation():
    """Test that natural frequency can be computed from eigenvalues."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 2

    lti = LTIStability(rotor)
    lti = lti.eigen_full_range()

    for modal_sol in lti.modal_solution:
        damping = modal_sol.damping
        frequency = modal_sol.frequency

        # Natural frequency: ω_n = sqrt(σ² + ω²)
        natural_freq = np.sqrt(damping**2 + frequency**2)

        # Natural frequency should be non-negative
        assert np.all(natural_freq >= 0)

        # For oscillatory modes (frequency != 0), natural frequency > 0
        oscillatory_mask = np.abs(frequency) > 1e-10
        if np.any(oscillatory_mask):
            assert np.all(natural_freq[oscillatory_mask] > 0)


# ============================================================================
# Test Integration with run_stability
# ============================================================================

def test_run_stability_lti():
    """Test LTI stability through run_stability factory method."""
    from src.stability_analysis.base import StabilityAnalysis

    rotor = RotorBuild.build_all()
    rotor.problem.required_solver = "LTI"
    rotor.problem.continuation = "NO"
    rotor.problem.number_points = 5

    result = StabilityAnalysis.run_stability(rotor)

    # Should return LTIStability instance
    assert isinstance(result, LTIStability)

    # Should have computed results
    assert len(result.modal_solution) == 5

    # Check results are populated
    for modal_sol in result.modal_solution:
        assert len(modal_sol.damping) > 0


# ============================================================================
# Test Edge Cases
# ============================================================================

def test_eigen_single_point_very_small_omega():
    """Test eigen_single_point with very small but positive OMEGA."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    # Should work with very small positive values
    eigensol = lti.eigen_single_point(1e-3)
    assert isinstance(eigensol, dict)
    assert 'eigenvalues' in eigensol


def test_eigen_single_point_very_large_omega():
    """Test eigen_single_point with very large OMEGA."""
    rotor = RotorBuild.build_all()
    lti = LTIStability(rotor)

    # Should work with large values
    eigensol = lti.eigen_single_point(1000.0)
    assert isinstance(eigensol, dict)
    assert 'eigenvalues' in eigensol


def test_consistency_single_vs_range():
    """Test that eigen_single_point matches eigen_full_range at same OMEGA."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 1000.0
    rotor.problem.higher_rotor_RPM = 2000.0
    rotor.problem.number_points = 5

    lti = LTIStability(rotor)

    # Get full range results
    lti = lti.eigen_full_range()

    # Compare with single point computation
    for modal_sol in lti.modal_solution:
        single_result = lti.eigen_single_point(modal_sol.OMEGA)

        # Eigenvalues should match
        assert np.allclose(
            single_result['eigenvalues'],
            modal_sol.eigensolution['eigenvalues'],
            rtol=1e-12
        ), f"Single point and full range should give same eigenvalues at OMEGA={modal_sol.OMEGA}"


# ============================================================================
# Test Main Block Execution
# ============================================================================

def test_main_block_execution():
    """Test that __main__ block executes without errors."""
    # This ensures the example code in the module works
    import subprocess

    result = subprocess.run(
        [sys.executable, "-m", "src.stability_analysis.lti_stability"],
        capture_output=True,
        text=True,
        timeout=30
    )

    # Should complete successfully
    assert result.returncode == 0, f"Main block failed: {result.stderr}"

    # Should print expected output
    assert "Building rotor system" in result.stdout or "Building rotor system" in result.stderr
    assert "Running LTI stability analysis" in result.stdout or "Running LTI stability analysis" in result.stderr


# ============================================================================
# Run tests manually
# ============================================================================

if __name__ == "__main__":
    import traceback

    # List all test functions
    test_functions = [
        test_lti_initialization,
        test_lti_inherits_from_base,
        test_eigen_single_point_basic,
        test_eigen_single_point_eigenvalues_complex,
        test_eigen_single_point_eigenvectors_shape,
        test_eigen_single_point_different_speeds,
        test_eigen_single_point_negative_omega,
        test_eigen_single_point_zero_omega,
        test_eigen_single_point_non_numeric_omega,
        test_eigen_single_point_numpy_scalar,
        test_eigen_full_range_basic,
        test_eigen_full_range_modal_solution_structure,
        test_eigen_full_range_rpm_values,
        test_eigen_full_range_damping_frequency_relationship,
        test_eigen_full_range_method_chaining,
        test_eigen_full_range_single_point,
        test_eigen_full_range_many_points,
        test_eigen_full_range_l2r_direction,
        test_eigen_full_range_r2l_direction,
        test_eigenvalue_properties,
        test_damping_interpretation,
        test_natural_frequency_computation,
        test_run_stability_lti,
        test_eigen_single_point_very_small_omega,
        test_eigen_single_point_very_large_omega,
        test_consistency_single_vs_range,
        test_main_block_execution,
    ]

    passed = 0
    failed = 0

    for test in test_functions:
        try:
            test()
            print(f"PASS {test.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL {test.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"FAIL {test.__name__}: {type(e).__name__}: {e}")
            traceback.print_exc()
            failed += 1

    print(f"\n{'='*60}")
    print(f"Test Results: {passed} passed, {failed} failed")
    print(f"{'='*60}")
    sys.exit(0 if failed == 0 else 1)
