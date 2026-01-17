"""Unit tests for rotor_build module.

Tests the RotorBuild class for constructing state-space matrices.
Validates mass matrix, damping/stiffness matrices, and state matrix assembly.

Run with:
    python -m pytest tests/test_rotor_build.py -v
    python tests/test_rotor_build.py
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
from src.rotor_build import RotorBuild
from src.config.problem_definition import ProblemDefinition


def test_rotor_build_initialization():
    """Test RotorBuild initialization with auto_build."""
    problem = ProblemDefinition()
    rotor = RotorBuild(problem, auto_build=True)

    assert rotor.problem is problem, "Problem should be stored"
    assert rotor.state_matrix_function is not None, "State matrix function should be initialized"


def test_rotor_build_initialization_no_auto_build():
    """Test RotorBuild initialization without auto_build."""
    problem = ProblemDefinition()
    rotor = RotorBuild(problem, auto_build=False)

    assert rotor.problem is problem, "Problem should be stored"
    assert rotor.state_matrix_function is None, "State matrix function should not be initialized"


def test_build_all_factory_method():
    """Test the build_all factory method."""
    rotor = RotorBuild.build_all()

    assert rotor is not None, "RotorBuild object should be created"
    assert rotor.problem is not None, "Problem definition should exist"
    assert rotor.state_matrix_function is not None, "State matrix function should be initialized"

    # Check default problem values
    assert rotor.problem.number_blades == 4, "Default should be 4 blades"
    assert rotor.problem.damper_connection == "H2B", "Default should be H2B"


def test_state_matrix_function_exists():
    """Test that state matrix function is created."""
    rotor = RotorBuild.build_all()

    assert rotor.state_matrix_function is not None, \
        "State matrix function should be initialized"
    assert callable(rotor.state_matrix_function), \
        "State matrix function should be callable"


def test_state_matrix_dimensions():
    """Test state matrix dimensions."""
    rotor = RotorBuild.build_all()

    # Evaluate state matrix
    A = rotor.state_matrix_function(t=0.0, rotor_OMEGA=25.0)

    # State matrix should be square and even-dimensional (2n x 2n)
    assert A.shape[0] == A.shape[1], "State matrix should be square"
    assert A.shape[0] % 2 == 0, "State matrix dimension should be even (2n)"
    assert A.shape[0] > 0, "State matrix should have positive size"


def test_state_matrix_structure():
    """Test state matrix block structure."""
    rotor = RotorBuild.build_all()
    A = rotor.state_matrix_function(t=0.0, rotor_OMEGA=25.0)
    n = A.shape[0] // 2  # State matrix is 2n x 2n

    # Extract blocks
    # A = [[-M^-1 C, -M^-1 K],
    #      [   I,        0   ]]
    upper_left = A[:n, :n]
    upper_right = A[:n, n:]
    lower_left = A[n:, :n]
    lower_right = A[n:, n:]

    # Lower left should be identity
    assert np.allclose(lower_left, np.eye(n)), \
        "Lower left block should be identity matrix"

    # Lower right should be zero
    assert np.allclose(lower_right, np.zeros((n, n))), \
        "Lower right block should be zero matrix"

    # Upper blocks should not be all zero (sanity check)
    assert not np.allclose(upper_left, np.zeros((n, n))), \
        "Upper left block should not be all zeros"


def test_state_matrix_at_different_speeds():
    """Test state matrix evaluation at different rotor speeds."""
    rotor = RotorBuild.build_all()

    speeds = [10.0, 25.0, 50.0, 100.0]  # rad/s

    for omega in speeds:
        A = rotor.state_matrix_function(t=0.0, rotor_OMEGA=omega)

        # Should return valid matrix
        assert A is not None, f"State matrix should exist at omega={omega}"
        assert not np.any(np.isnan(A)), f"State matrix should not contain NaN at omega={omega}"
        assert not np.any(np.isinf(A)), f"State matrix should not contain Inf at omega={omega}"


def test_state_matrix_at_different_times():
    """Test state matrix evaluation at different times (for time-varying systems)."""
    rotor = RotorBuild.build_all()

    times = [0.0, 0.1, 0.5, 1.0]  # seconds
    omega = 25.0  # rad/s

    for t in times:
        A = rotor.state_matrix_function(t=t, rotor_OMEGA=omega)

        # Should return valid matrix
        assert A is not None, f"State matrix should exist at t={t}"
        assert not np.any(np.isnan(A)), f"State matrix should not contain NaN at t={t}"
        assert not np.any(np.isinf(A)), f"State matrix should not contain Inf at t={t}"


def test_different_damper_activation_modes():
    """Test different damper activation modes produce different matrices."""
    # Test ALL vs ODI activation modes with H2B (which is well-supported)
    problem_all = ProblemDefinition(
        damper_connection="H2B",
        damper_activation="ALL",
        number_blades=4
    )
    problem_odi = ProblemDefinition(
        damper_connection="H2B",
        damper_activation="ODI",
        number_blades=4
    )

    rotor_all = RotorBuild(problem_all, auto_build=True)
    rotor_odi = RotorBuild(problem_odi, auto_build=True)

    # Evaluate state matrices
    A_all = rotor_all.state_matrix_function(t=0.0, rotor_OMEGA=25.0)
    A_odi = rotor_odi.state_matrix_function(t=0.0, rotor_OMEGA=25.0)

    # Matrices should be different (different damper activation)
    assert not np.allclose(A_all, A_odi), \
        "ALL and ODI damper activation should produce different state matrices"


def test_state_matrix_numerical_properties():
    """Test numerical properties of state matrix."""
    rotor = RotorBuild.build_all()
    A = rotor.state_matrix_function(t=0.0, rotor_OMEGA=25.0)

    # Check condition number is reasonable (not too ill-conditioned)
    cond = np.linalg.cond(A)
    assert cond < 1e10, f"Condition number {cond:.2e} is too large"

    # Matrix should be real-valued
    assert np.all(np.isreal(A)), "State matrix should be real-valued"


if __name__ == "__main__":
    # Run tests manually without pytest
    import traceback

    tests = [
        test_rotor_build_initialization,
        test_rotor_build_initialization_no_auto_build,
        test_build_all_factory_method,
        test_state_matrix_function_exists,
        test_state_matrix_dimensions,
        test_state_matrix_structure,
        test_state_matrix_at_different_speeds,
        test_state_matrix_at_different_times,
        test_different_damper_activation_modes,
        test_state_matrix_numerical_properties,
    ]

    passed = 0
    failed = 0

    for test in tests:
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

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)
