"""Unit tests for stability_analysis.hd_stability module.

Tests the HDStability class and its methods including:
- hd_single_point (single operating point with HD analysis)
- hd_full_range (full RPM range analysis)
- HD matrix construction and eigenvalue extraction
- Input validation
- Error handling
- Integration with real rotor systems

Run with:
    python -m pytest tests/test_hd_stability.py -v
    python tests/test_hd_stability.py
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pytest
from src.stability_analysis.hd_stability import HDStability
from src.rotor_build import RotorBuild


# ============================================================================
# Test HDStability Initialization
# ============================================================================

def test_hd_initialization():
    """Test HDStability initialization."""
    rotor = RotorBuild.build_all()
    hd = HDStability(rotor)

    assert hd.rotor_build is rotor
    assert isinstance(hd.modal_solution, list)
    assert len(hd.modal_solution) == 0


def test_hd_inherits_from_stability_analysis():
    """Test that HDStability properly inherits from StabilityAnalysis.

    Note: HDStability inherits directly from StabilityAnalysis, not LTPStability.
    HD and LTP are parallel approaches for time-periodic systems, not parent/child.
    """
    from src.stability_analysis.base import StabilityAnalysis

    rotor = RotorBuild.build_all()
    hd = HDStability(rotor)

    assert isinstance(hd, StabilityAnalysis)
    # Should have inherited methods from StabilityAnalysis
    assert hasattr(hd, 'assign_range_OMEGA')
    assert hasattr(hd, 'assign_period_T')
    assert hasattr(hd, 'hd_computer')  # Static method from base class


# ============================================================================
# Test hd_single_point - Basic Functionality
# ============================================================================

def test_hd_single_point_basic():
    """Test hd_single_point with valid rotor system."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_harmonics = 1
    rotor.problem.time_samples = 50
    rotor.problem.hd_use_complex = False

    hd = HDStability(rotor)

    OMEGA = 25.0  # rad/s
    eigensol = hd.hd_single_point(OMEGA)

    # Check return type
    assert isinstance(eigensol, dict)
    assert 'eigenvalues' in eigensol
    assert 'eigenvectors' in eigensol

    # Check array types
    assert isinstance(eigensol['eigenvalues'], np.ndarray)
    assert isinstance(eigensol['eigenvectors'], np.ndarray)

    # Check dimensions - HD expands the system
    n_states = 2 * rotor.mass_matrix.shape[0]
    n_harmonics = 1
    expected_size = (1 + 2 * n_harmonics) * n_states
    assert eigensol['eigenvalues'].shape == (expected_size,)
    assert eigensol['eigenvectors'].shape == (expected_size, expected_size)


def test_hd_single_point_multiple_harmonics():
    """Test hd_single_point with multiple harmonics."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_harmonics = 3
    rotor.problem.time_samples = 100
    rotor.problem.hd_use_complex = False

    hd = HDStability(rotor)

    OMEGA = 50.0
    eigensol = hd.hd_single_point(OMEGA)

    # Check expanded system size
    n_states = 2 * rotor.mass_matrix.shape[0]
    n_harmonics = 3
    expected_size = (1 + 2 * n_harmonics) * n_states
    assert eigensol['eigenvalues'].shape == (expected_size,)


def test_hd_single_point_complex_formulation():
    """Test hd_single_point with complex formulation."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_harmonics = 2
    rotor.problem.time_samples = 80
    rotor.problem.hd_use_complex = True

    hd = HDStability(rotor)

    OMEGA = 30.0
    eigensol = hd.hd_single_point(OMEGA)

    # Should return eigenvalues (may be complex)
    assert isinstance(eigensol['eigenvalues'], np.ndarray)
    assert isinstance(eigensol['eigenvectors'], np.ndarray)

    # Check dimensions match expected size
    n_states = 2 * rotor.mass_matrix.shape[0]
    n_harmonics = 2
    expected_size = (1 + 2 * n_harmonics) * n_states
    assert eigensol['eigenvalues'].shape == (expected_size,)


def test_hd_real_vs_complex_comparison():
    """Test that real and complex formulations give similar results."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_harmonics = 2
    rotor.problem.time_samples = 100

    OMEGA = 25.0

    # Test with real formulation
    rotor.problem.hd_use_complex = False
    hd_real = HDStability(rotor)
    eigensol_real = hd_real.hd_single_point(OMEGA)

    # Test with complex formulation
    rotor.problem.hd_use_complex = True
    hd_complex = HDStability(rotor)
    eigensol_complex = hd_complex.hd_single_point(OMEGA)

    # Both should produce same dimension results
    assert eigensol_real['eigenvalues'].shape == eigensol_complex['eigenvalues'].shape

    # Extract damping (real parts) - the stability indicators
    damping_real = np.real(eigensol_real['eigenvalues'])
    damping_complex = np.real(eigensol_complex['eigenvalues'])

    # Sort damping values for comparison
    damping_real_sorted = np.sort(damping_real)
    damping_complex_sorted = np.sort(damping_complex)

    # Damping values (stability indicators) should be very close
    # Both methods approximate the same system
    assert np.allclose(damping_real_sorted, damping_complex_sorted, rtol=1e-3, atol=1e-6)

    # Extract frequencies (imaginary parts) - should also be similar in magnitude
    freq_real = np.abs(np.imag(eigensol_real['eigenvalues']))
    freq_complex = np.abs(np.imag(eigensol_complex['eigenvalues']))

    # Sort frequencies for comparison
    freq_real_sorted = np.sort(freq_real)
    freq_complex_sorted = np.sort(freq_complex)

    # Frequencies should be similar (may have small differences due to truncation)
    assert np.allclose(freq_real_sorted, freq_complex_sorted, rtol=1e-2, atol=1e-4)


def test_hd_full_range_with_complex():
    """Test hd_full_range with complex formulation."""
    rotor = RotorBuild.build_all()
    rotor.problem.required_solver = "HD"
    rotor.problem.number_harmonics = 2
    rotor.problem.number_points = 3
    rotor.problem.hd_use_complex = True

    hd = HDStability(rotor)
    hd = hd.hd_full_range()

    # Check that all solutions were computed
    assert len(hd.modal_solution) == 3

    # Check each solution has results
    for modal_sol in hd.modal_solution:
        assert hasattr(modal_sol, 'damping')
        assert hasattr(modal_sol, 'frequency')
        assert len(modal_sol.damping) > 0
        assert len(modal_sol.frequency) > 0


# ============================================================================
# Test hd_single_point - Input Validation
# ============================================================================

def test_hd_single_point_negative_omega():
    """Test that negative OMEGA raises ValueError."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_harmonics = 1
    hd = HDStability(rotor)

    with pytest.raises(ValueError, match="OMEGA must be positive"):
        hd.hd_single_point(-10.0)


def test_hd_single_point_zero_omega():
    """Test that zero OMEGA raises ValueError."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_harmonics = 1
    hd = HDStability(rotor)

    with pytest.raises(ValueError, match="OMEGA must be positive"):
        hd.hd_single_point(0.0)


def test_hd_single_point_invalid_harmonics():
    """Test that invalid number_harmonics raises ValueError."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_harmonics = 0  # Invalid
    rotor.problem.time_samples = 50
    hd = HDStability(rotor)

    OMEGA = 25.0

    with pytest.raises(ValueError, match="number_harmonics must be integer >= 1"):
        hd.hd_single_point(OMEGA)


def test_hd_single_point_invalid_time_samples():
    """Test that invalid time_samples raises ValueError."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_harmonics = 1
    rotor.problem.time_samples = 5  # Too few
    hd = HDStability(rotor)

    OMEGA = 25.0

    with pytest.raises(ValueError, match="time_samples must be integer >= 10"):
        hd.hd_single_point(OMEGA)


# ============================================================================
# Test hd_full_range - Basic Functionality
# ============================================================================

def test_hd_full_range_basic():
    """Test hd_full_range with default rotor configuration."""
    rotor = RotorBuild.build_all()
    rotor.problem.required_solver = "HD"
    rotor.problem.number_harmonics = 1
    rotor.problem.number_points = 2  # Keep small for speed

    hd = HDStability(rotor)
    hd = hd.hd_full_range()

    # Check that modal_solution is populated
    assert len(hd.modal_solution) == 2


def test_hd_full_range_modal_solution_structure():
    """Test that modal_solution has correct structure."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 1000.0
    rotor.problem.higher_rotor_RPM = 2000.0
    rotor.problem.number_points = 3
    rotor.problem.number_harmonics = 1

    hd = HDStability(rotor)
    hd = hd.hd_full_range()

    # Check each modal solution
    for modal_sol in hd.modal_solution:
        # Check required attributes
        assert hasattr(modal_sol, 'OMEGA')
        assert hasattr(modal_sol, 'OMEGA_RPM')
        assert hasattr(modal_sol, 'T')
        assert hasattr(modal_sol, 'damping')
        assert hasattr(modal_sol, 'frequency')
        assert hasattr(modal_sol, 'eigensolution')

        # Check arrays are populated
        assert len(modal_sol.damping) > 0
        assert len(modal_sol.frequency) > 0

        # Check eigensolution structure
        assert 'eigenvalues' in modal_sol.eigensolution
        assert 'eigenvectors' in modal_sol.eigensolution


def test_hd_full_range_method_chaining():
    """Test that hd_full_range returns self for method chaining."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 2
    rotor.problem.number_harmonics = 1

    hd = HDStability(rotor)
    result = hd.hd_full_range()

    assert result is hd, "Method should return self for chaining"


# ============================================================================
# Test HD Properties
# ============================================================================

def test_hd_expanded_system_size():
    """Test that HD matrix has correct expanded size."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_harmonics = 2
    rotor.problem.time_samples = 50
    rotor.problem.number_points = 1

    hd = HDStability(rotor)
    hd = hd.hd_full_range()

    # Original state dimension
    n_original = 2 * rotor.mass_matrix.shape[0]
    # Number of harmonics
    N = 2
    # Expected expanded dimension
    expected_size = (1 + 2 * N) * n_original

    # Check eigenvalue count
    assert len(hd.modal_solution[0].damping) == expected_size


def test_hd_damping_frequency_relationship():
    """Test that damping and frequency are from eigenvalues."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 2
    rotor.problem.number_harmonics = 1

    hd = HDStability(rotor)
    hd = hd.hd_full_range()

    for modal_sol in hd.modal_solution:
        eigenvalues = modal_sol.eigensolution['eigenvalues']
        damping = modal_sol.damping
        frequency = modal_sol.frequency

        # Check relationships
        assert np.allclose(damping, np.real(eigenvalues))
        assert np.allclose(frequency, np.imag(eigenvalues))


# ============================================================================
# Test Integration with run_stability
# ============================================================================

def test_run_stability_hd():
    """Test HD stability through run_stability factory method."""
    from src.stability_analysis.base import StabilityAnalysis

    rotor = RotorBuild.build_all()
    rotor.problem.required_solver = "HD"
    rotor.problem.continuation = "NO"
    rotor.problem.number_points = 2
    rotor.problem.number_harmonics = 1

    result = StabilityAnalysis.run_stability(rotor)

    # Should return HDStability instance
    assert isinstance(result, HDStability)

    # Should have computed results
    assert len(result.modal_solution) == 2


# ============================================================================
# Test Main Block Execution
# ============================================================================

def test_main_block_execution():
    """Test that __main__ block executes without errors."""
    import subprocess

    result = subprocess.run(
        [sys.executable, "-m", "src.stability_analysis.hd_stability"],
        capture_output=True,
        text=True,
        timeout=60
    )

    # Should complete successfully
    assert result.returncode == 0, f"Main block failed: {result.stderr}"

    # Should print expected output
    assert "Building rotor system" in result.stdout or "Building rotor system" in result.stderr
    assert "Running HD stability analysis" in result.stdout or "Running HD stability analysis" in result.stderr


# ============================================================================
# Run tests manually
# ============================================================================

if __name__ == "__main__":
    import traceback

    # List all test functions
    test_functions = [
        test_hd_initialization,
        test_hd_inherits_from_ltp,
        test_hd_single_point_basic,
        test_hd_single_point_multiple_harmonics,
        test_hd_single_point_complex_formulation,
        test_hd_real_vs_complex_comparison,
        test_hd_full_range_with_complex,
        test_hd_single_point_negative_omega,
        test_hd_single_point_zero_omega,
        test_hd_single_point_negative_T,
        test_hd_single_point_invalid_harmonics,
        test_hd_single_point_invalid_time_samples,
        test_hd_full_range_basic,
        test_hd_full_range_modal_solution_structure,
        test_hd_full_range_method_chaining,
        test_hd_expanded_system_size,
        test_hd_damping_frequency_relationship,
        test_run_stability_hd,
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
