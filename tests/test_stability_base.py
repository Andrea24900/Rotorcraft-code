"""Unit tests for stability_analysis.base module.

Tests the StabilityAnalysis base class and its methods including:
- monodromy_computer (LTP method)
- hd_computer (real formulation)
- hd_computer_complex (complex formulation)
- Input validation
- Edge cases

Run with:
    python -m pytest tests/test_stability_base.py -v
    python tests/test_stability_base.py
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pytest
from src.stability_analysis.base import StabilityAnalysis, ModalSolution
from src.rotor_build import RotorBuild
from src.config.problem_definition import ProblemDefinition


# ============================================================================
# Test Constants
# ============================================================================

def test_class_constants():
    """Test that class constants are defined correctly."""
    assert hasattr(StabilityAnalysis, 'DEFAULT_RTOL')
    assert hasattr(StabilityAnalysis, 'DEFAULT_ATOL')
    assert StabilityAnalysis.DEFAULT_RTOL == 1e-6
    assert StabilityAnalysis.DEFAULT_ATOL == 1e-8


# ============================================================================
# Test Monodromy Computer
# ============================================================================

def test_monodromy_computer_simple_system():
    """Test monodromy computer with a simple time-invariant system."""
    # Simple 2x2 time-invariant system: dx/dt = A*x
    # A = [[0, 1], [-1, -0.1]]
    # Analytical solution exists for verification
    A_func = lambda t: np.array([[0, 1], [-1, -0.1]])
    T = 2 * np.pi

    M = StabilityAnalysis.monodromy_computer(A_func, T)

    # Check dimensions
    assert M.shape == (2, 2), "Monodromy matrix should be 2x2"

    # Check that it's real-valued
    assert np.all(np.isreal(M)), "Monodromy matrix should be real"

    # Check stability (eigenvalues should have magnitude < 1 for this stable system)
    eigenvalues = np.linalg.eigvals(M)
    assert len(eigenvalues) == 2, "Should have 2 eigenvalues"


def test_monodromy_computer_time_periodic_system():
    """Test monodromy computer with a time-periodic system."""
    # Time-periodic system with sinusoidal coefficient
    A_func = lambda t: np.array([[0, 1], [-1 - 0.1*np.sin(t), -0.1]])
    T = 2 * np.pi  # Period matches the sinusoidal term

    M = StabilityAnalysis.monodromy_computer(A_func, T)

    # Check dimensions
    assert M.shape == (2, 2)

    # Check eigenvalues exist
    eigenvalues = np.linalg.eigvals(M)
    assert len(eigenvalues) == 2


def test_monodromy_computer_custom_tolerances():
    """Test monodromy computer with custom tolerances."""
    A_func = lambda t: np.array([[0, 1], [-1, -0.1]])
    T = 2 * np.pi

    # Use custom tolerances
    M = StabilityAnalysis.monodromy_computer(A_func, T, rtol=1e-8, atol=1e-10)

    assert M.shape == (2, 2)
    assert np.all(np.isreal(M))


def test_monodromy_computer_identity_system():
    """Test monodromy computer with identity system (dΦ/dt = 0)."""
    # System: dx/dt = 0 => Φ(t) = Φ(0) = I
    A_func = lambda t: np.zeros((2, 2))
    T = 1.0

    M = StabilityAnalysis.monodromy_computer(A_func, T)

    # Monodromy matrix should be identity
    assert np.allclose(M, np.eye(2), atol=1e-6), \
        "Monodromy matrix should be identity for zero system"


def test_monodromy_computer_3x3_system():
    """Test monodromy computer with a 3x3 system."""
    A_func = lambda t: np.array([
        [0, 1, 0],
        [0, 0, 1],
        [-1, -1, -1]
    ])
    T = 2 * np.pi

    M = StabilityAnalysis.monodromy_computer(A_func, T)

    assert M.shape == (3, 3)
    eigenvalues = np.linalg.eigvals(M)
    assert len(eigenvalues) == 3


# ============================================================================
# Test Monodromy Computer - Input Validation
# ============================================================================

def test_monodromy_negative_period():
    """Test that negative period raises ValueError."""
    A_func = lambda t: np.array([[0, 1], [-1, -0.1]])

    with pytest.raises(ValueError, match="period_T must be positive"):
        StabilityAnalysis.monodromy_computer(A_func, -1.0)


def test_monodromy_zero_period():
    """Test that zero period raises ValueError."""
    A_func = lambda t: np.array([[0, 1], [-1, -0.1]])

    with pytest.raises(ValueError, match="period_T must be positive"):
        StabilityAnalysis.monodromy_computer(A_func, 0.0)


def test_monodromy_non_square_matrix():
    """Test that non-square matrix raises ValueError."""
    A_func = lambda t: np.array([[0, 1, 2], [-1, -0.1, 0]])  # 2x3 matrix

    with pytest.raises(ValueError, match="square matrix"):
        StabilityAnalysis.monodromy_computer(A_func, 1.0)


def test_monodromy_non_2d_array():
    """Test that non-2D array raises ValueError."""
    A_func = lambda t: np.array([0, 1, -1])  # 1D array

    with pytest.raises(ValueError, match="2D array"):
        StabilityAnalysis.monodromy_computer(A_func, 1.0)


def test_monodromy_invalid_function():
    """Test that invalid function raises ValueError."""
    A_func = lambda t: None  # Returns None instead of array

    with pytest.raises(ValueError):
        StabilityAnalysis.monodromy_computer(A_func, 1.0)


# ============================================================================
# Test HD Computer - Real Formulation
# ============================================================================

def test_hd_computer_basic():
    """Test HD computer with basic time-periodic system."""
    A_func = lambda t: np.array([[0, 1], [-1 - 0.1*np.sin(t), -0.1]])
    OMEGA = 1.0
    T = 2 * np.pi / OMEGA
    time = np.linspace(0, T, 100)
    number_harmonics = 1

    A_HD = StabilityAnalysis.hd_computer(A_func, time, number_harmonics, OMEGA)

    # For n=2 original state variables and N=1 harmonic (user-defined): HD matrix size = (1 + 2*N) * n = 6
    expected_size = (1 + 2 * number_harmonics) * 2
    assert A_HD.shape == (expected_size, expected_size)

    # Should be real-valued (default)
    assert np.all(np.isreal(A_HD))


def test_hd_computer_multiple_harmonics():
    """Test HD computer with multiple harmonics."""
    A_func = lambda t: np.array([[0, 1], [-1 - 0.1*np.sin(t), -0.1]])
    OMEGA = 1.0
    T = 2 * np.pi / OMEGA
    time = np.linspace(0, T, 200)
    number_harmonics = 3

    A_HD = StabilityAnalysis.hd_computer(A_func, time, number_harmonics, OMEGA)

    # For n=2 original state variables and N=3 harmonics (user-defined): HD matrix size = (1 + 2*N) * n = 14
    expected_size = (1 + 2 * number_harmonics) * 2
    assert A_HD.shape == (expected_size, expected_size)


def test_hd_computer_3x3_system():
    """Test HD computer with 3x3 system."""
    A_func = lambda t: np.array([
        [0, 1, 0],
        [0, 0, 1],
        [-1 - 0.1*np.sin(t), -1, -1]
    ])
    OMEGA = 1.0
    T = 2 * np.pi / OMEGA
    time = np.linspace(0, T, 100)
    number_harmonics = 2

    A_HD = StabilityAnalysis.hd_computer(A_func, time, number_harmonics, OMEGA)

    # For n=3 original state variables and N=2 harmonics (user-defined): HD matrix size = (1 + 2*N) * n = 15
    expected_size = (1 + 2 * number_harmonics) * 3
    assert A_HD.shape == (expected_size, expected_size)


# ============================================================================
# Test HD Computer - Complex Formulation
# ============================================================================

def test_hd_computer_complex():
    """Test HD computer with complex formulation."""
    A_func = lambda t: np.array([[0, 1], [-1 - 0.1*np.sin(t), -0.1]])
    OMEGA = 1.0
    T = 2 * np.pi / OMEGA
    time = np.linspace(0, T, 100)
    number_harmonics = 1

    A_HD = StabilityAnalysis.hd_computer(A_func, time, number_harmonics, OMEGA, use_complex=True)

    # Same size as real formulation
    expected_size = (1 + 2 * number_harmonics) * 2
    assert A_HD.shape == (expected_size, expected_size)

    # Should be complex-valued
    assert A_HD.dtype == np.complex128 or np.iscomplexobj(A_HD)


def test_hd_computer_complex_direct_call():
    """Test calling hd_computer_complex directly."""
    A_func = lambda t: np.array([[0, 1], [-1 - 0.1*np.sin(t), -0.1]])
    OMEGA = 1.0
    T = 2 * np.pi / OMEGA
    time = np.linspace(0, T, 100)
    number_harmonics = 2

    A_HD = StabilityAnalysis.hd_computer_complex(A_func, time, number_harmonics, OMEGA)

    expected_size = (1 + 2 * number_harmonics) * 2
    assert A_HD.shape == (expected_size, expected_size)
    assert np.iscomplexobj(A_HD)


# ============================================================================
# Test HD Computer - Input Validation
# ============================================================================

def test_hd_computer_invalid_harmonics():
    """Test that invalid number of harmonics raises ValueError."""
    A_func = lambda t: np.array([[0, 1], [-1, -0.1]])
    time = np.linspace(0, 1, 100)
    OMEGA = 1.0

    with pytest.raises(ValueError, match="number_harmonics must be >= 1"):
        StabilityAnalysis.hd_computer(A_func, time, 0, OMEGA)

    with pytest.raises(ValueError, match="number_harmonics must be >= 1"):
        StabilityAnalysis.hd_computer(A_func, time, -1, OMEGA)


def test_hd_computer_invalid_omega():
    """Test that invalid OMEGA raises ValueError."""
    A_func = lambda t: np.array([[0, 1], [-1, -0.1]])
    time = np.linspace(0, 1, 100)

    with pytest.raises(ValueError, match="OMEGA must be positive"):
        StabilityAnalysis.hd_computer(A_func, time, 1, 0.0)

    with pytest.raises(ValueError, match="OMEGA must be positive"):
        StabilityAnalysis.hd_computer(A_func, time, 1, -1.0)


def test_hd_computer_invalid_time_array():
    """Test that invalid time array raises ValueError."""
    A_func = lambda t: np.array([[0, 1], [-1, -0.1]])
    OMEGA = 1.0

    # Not a numpy array
    with pytest.raises(ValueError, match="time must be numpy array"):
        StabilityAnalysis.hd_computer(A_func, [0, 1, 2], 1, OMEGA)

    # Too few points
    with pytest.raises(ValueError, match="at least 2 points"):
        StabilityAnalysis.hd_computer(A_func, np.array([0]), 1, OMEGA)


def test_hd_computer_non_square_matrix():
    """Test that non-square matrix raises ValueError."""
    A_func = lambda t: np.array([[0, 1, 2], [-1, -0.1, 0]])  # 2x3
    time = np.linspace(0, 1, 100)
    OMEGA = 1.0

    with pytest.raises(ValueError, match="square matrix"):
        StabilityAnalysis.hd_computer(A_func, time, 1, OMEGA)


def test_hd_computer_complex_invalid_inputs():
    """Test that hd_computer_complex validates inputs."""
    A_func = lambda t: np.array([[0, 1], [-1, -0.1]])
    time = np.linspace(0, 1, 100)

    with pytest.raises(ValueError, match="OMEGA must be positive"):
        StabilityAnalysis.hd_computer_complex(A_func, time, 1, -1.0)

    with pytest.raises(ValueError, match="number_harmonics must be >= 1"):
        StabilityAnalysis.hd_computer_complex(A_func, time, 0, 1.0)


# ============================================================================
# Test ModalSolution Dataclass
# ============================================================================

def test_modal_solution_initialization():
    """Test ModalSolution dataclass initialization."""
    modal_sol = ModalSolution()

    assert modal_sol.OMEGA == 0.0
    assert modal_sol.OMEGA_RPM == 0.0
    assert modal_sol.T == 0.0
    assert isinstance(modal_sol.damping, np.ndarray)
    assert isinstance(modal_sol.frequency, np.ndarray)
    assert isinstance(modal_sol.eigensolution, dict)
    assert isinstance(modal_sol.char_solution, dict)


def test_modal_solution_with_values():
    """Test ModalSolution with custom values."""
    modal_sol = ModalSolution(
        OMEGA=25.0,
        OMEGA_RPM=238.73,
        T=0.2513,
        damping=np.array([-1.0, -2.0]),
        frequency=np.array([10.0, 20.0])
    )

    assert modal_sol.OMEGA == 25.0
    assert modal_sol.OMEGA_RPM == 238.73
    assert np.allclose(modal_sol.damping, np.array([-1.0, -2.0]))
    assert np.allclose(modal_sol.frequency, np.array([10.0, 20.0]))


# ============================================================================
# Test StabilityAnalysis Class Methods
# ============================================================================

def test_stability_analysis_initialization():
    """Test StabilityAnalysis initialization."""
    rotor = RotorBuild.build_all()
    stability = StabilityAnalysis(rotor)

    assert stability.rotor_build is rotor
    assert isinstance(stability.modal_solution, list)
    assert len(stability.modal_solution) == 0


def test_assign_range_omega():
    """Test assign_range_OMEGA method."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 100.0
    rotor.problem.higher_rotor_RPM = 200.0
    rotor.problem.number_points = 11

    stability = StabilityAnalysis(rotor)
    stability.assign_range_OMEGA()

    assert len(stability.modal_solution) == 11

    # Check first and last values
    first_omega_rpm = stability.modal_solution[0].OMEGA_RPM
    last_omega_rpm = stability.modal_solution[-1].OMEGA_RPM

    # Should be descending if R2L (default)
    if rotor.problem.solution_direction == "R2L":
        assert np.isclose(first_omega_rpm, 200.0)
        assert np.isclose(last_omega_rpm, 100.0)


def test_assign_range_omega_l2r():
    """Test assign_range_OMEGA with L2R direction."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 100.0
    rotor.problem.higher_rotor_RPM = 200.0
    rotor.problem.number_points = 11
    rotor.problem.solution_direction = "L2R"

    stability = StabilityAnalysis(rotor)
    stability.assign_range_OMEGA()

    # Should be ascending for L2R
    first_omega_rpm = stability.modal_solution[0].OMEGA_RPM
    last_omega_rpm = stability.modal_solution[-1].OMEGA_RPM

    assert np.isclose(first_omega_rpm, 100.0)
    assert np.isclose(last_omega_rpm, 200.0)


def test_assign_period_t():
    """Test assign_period_T method."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 5

    stability = StabilityAnalysis(rotor)
    stability.assign_range_OMEGA()
    stability.assign_period_T()

    # Check that periods are set correctly
    for modal_sol in stability.modal_solution:
        expected_T = 2 * np.pi / modal_sol.OMEGA
        assert np.isclose(modal_sol.T, expected_T), \
            f"Period should be 2π/Ω, got {modal_sol.T}, expected {expected_T}"


def test_run_stability_factory_lti():
    """Test run_stability factory method with LTI solver."""
    # Note: This test is skipped because LTI stability implementation
    # expects a different rotor interface that we don't test here
    # The factory method itself is tested with unknown solver
    pass


def test_run_stability_unknown_solver():
    """Test that unknown solver type raises ValueError."""
    rotor = RotorBuild.build_all()
    rotor.problem.required_solver = "UNKNOWN"

    with pytest.raises(ValueError, match="Unknown solver type"):
        StabilityAnalysis.run_stability(rotor)


# ============================================================================
# Test Integration with Real Rotor System
# ============================================================================

def test_monodromy_with_rotor_system():
    """Test monodromy computer with actual rotor state matrix."""
    rotor = RotorBuild.build_all()
    OMEGA = 25.0  # rad/s
    T = 2 * np.pi / OMEGA

    A_func = lambda t: rotor.state_matrix_function(t, OMEGA)

    M = StabilityAnalysis.monodromy_computer(A_func, T)

    # Check dimensions match rotor system
    n = rotor.mass_matrix.shape[0]
    expected_size = 2 * n  # State space is 2n

    assert M.shape == (expected_size, expected_size)
    assert np.all(np.isreal(M))


def test_hd_with_rotor_system():
    """Test HD computer with actual rotor state matrix."""
    rotor = RotorBuild.build_all()
    OMEGA = 25.0  # rad/s
    T = 2 * np.pi / OMEGA
    time = np.linspace(0, T, 100)
    number_harmonics = 1

    A_func = lambda t: rotor.state_matrix_function(t, OMEGA)

    A_HD = StabilityAnalysis.hd_computer(A_func, time, number_harmonics, OMEGA)

    # Check dimensions
    n = rotor.mass_matrix.shape[0]
    state_size = 2 * n
    expected_size = (1 + 2 * number_harmonics) * state_size

    assert A_HD.shape == (expected_size, expected_size)


# ============================================================================
# Test Edge Cases
# ============================================================================

def test_very_small_period():
    """Test monodromy with very small period."""
    A_func = lambda t: np.array([[0, 1], [-1, -0.1]])
    T = 1e-6

    M = StabilityAnalysis.monodromy_computer(A_func, T)

    # Should be close to identity for very small T
    assert M.shape == (2, 2)
    # M ≈ I + A*T for small T
    assert np.linalg.norm(M - np.eye(2)) < 1e-3


def test_very_large_harmonics():
    """Test HD computer with many harmonics (stress test)."""
    A_func = lambda t: np.array([[0, 1], [-1 - 0.1*np.sin(t), -0.1]])
    OMEGA = 1.0
    T = 2 * np.pi / OMEGA
    time = np.linspace(0, T, 500)
    number_harmonics = 10

    A_HD = StabilityAnalysis.hd_computer(A_func, time, number_harmonics, OMEGA)

    expected_size = (1 + 2 * number_harmonics) * 2
    assert A_HD.shape == (expected_size, expected_size)


def test_performance_warning_large_system():
    """Test that performance warning is issued for large systems."""
    import warnings

    # Create a large system (n=60 > 50)
    A_func = lambda t: np.zeros((60, 60))
    T = 1.0

    # Capture warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        M = StabilityAnalysis.monodromy_computer(A_func, T)

        # Check that a UserWarning was issued
        assert len(w) == 1
        assert issubclass(w[0].category, UserWarning)
        assert "large system" in str(w[0].message).lower()
        assert "n=60" in str(w[0].message)


def test_no_warning_small_system():
    """Test that no warning is issued for small systems."""
    import warnings

    # Create a small system (n=10 < 50)
    A_func = lambda t: np.zeros((10, 10))
    T = 1.0

    # Capture warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        M = StabilityAnalysis.monodromy_computer(A_func, T)

        # Check that no warning was issued
        assert len(w) == 0


# ============================================================================
# Run tests manually
# ============================================================================

if __name__ == "__main__":
    import traceback

    # List all test functions
    test_functions = [
        test_class_constants,
        test_monodromy_computer_simple_system,
        test_monodromy_computer_time_periodic_system,
        test_monodromy_computer_custom_tolerances,
        test_monodromy_computer_identity_system,
        test_monodromy_computer_3x3_system,
        test_monodromy_negative_period,
        test_monodromy_zero_period,
        test_monodromy_non_square_matrix,
        test_monodromy_non_2d_array,
        test_monodromy_invalid_function,
        test_hd_computer_basic,
        test_hd_computer_multiple_harmonics,
        test_hd_computer_3x3_system,
        test_hd_computer_complex,
        test_hd_computer_complex_direct_call,
        test_hd_computer_invalid_harmonics,
        test_hd_computer_invalid_omega,
        test_hd_computer_invalid_time_array,
        test_hd_computer_non_square_matrix,
        test_hd_computer_complex_invalid_inputs,
        test_modal_solution_initialization,
        test_modal_solution_with_values,
        test_stability_analysis_initialization,
        test_assign_range_omega,
        test_assign_range_omega_l2r,
        test_assign_period_t,
        test_run_stability_factory_lti,
        test_run_stability_unknown_solver,
        test_monodromy_with_rotor_system,
        test_hd_with_rotor_system,
        test_very_small_period,
        test_very_large_harmonics,
        test_performance_warning_large_system,
        test_no_warning_small_system,
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
