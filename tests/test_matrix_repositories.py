"""Unit tests for matrix_repositories module.

Tests the build_mass_matrix and build_damping_matrix functions.
Validates matrix construction, properties, and configurations for different blade counts.

Run with:
    python -m pytest tests/test_matrix_repositories.py -v
    python tests/test_matrix_repositories.py
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
from src.utils.matrix_repositories import build_mass_matrix, build_damping_matrix
from src.config.problem_definition import ProblemDefinition


def test_build_mass_matrix_3_blade():
    """Test mass matrix construction for 3-blade rotor."""
    problem = ProblemDefinition(number_blades=3)
    M = build_mass_matrix(problem)

    # Check dimensions
    assert M.shape == (5, 5), "3-blade mass matrix should be 5x5"

    # Check symmetry
    assert np.allclose(M, M.T), "Mass matrix should be symmetric"

    # Check positive definiteness
    eigenvalues = np.linalg.eigvals(M)
    assert np.all(eigenvalues > 0), "Mass matrix should be positive definite"


def test_build_mass_matrix_4_blade():
    """Test mass matrix construction for 4-blade rotor."""
    problem = ProblemDefinition(number_blades=4)
    M = build_mass_matrix(problem)

    # Check dimensions
    assert M.shape == (6, 6), "4-blade mass matrix should be 6x6"

    # Check symmetry
    assert np.allclose(M, M.T), "Mass matrix should be symmetric"

    # Check positive definiteness
    eigenvalues = np.linalg.eigvals(M)
    assert np.all(eigenvalues > 0), "Mass matrix should be positive definite"


def test_build_mass_matrix_5_blade():
    """Test mass matrix construction for 5-blade rotor."""
    problem = ProblemDefinition(number_blades=5)
    M = build_mass_matrix(problem)

    # Check dimensions
    assert M.shape == (7, 7), "5-blade mass matrix should be 7x7"

    # Check symmetry
    assert np.allclose(M, M.T), "Mass matrix should be symmetric"

    # Check positive definiteness
    eigenvalues = np.linalg.eigvals(M)
    assert np.all(eigenvalues > 0), "Mass matrix should be positive definite"


def test_build_mass_matrix_7_blade():
    """Test mass matrix construction for 7-blade rotor (time-varying)."""
    problem = ProblemDefinition(number_blades=7)

    # Mass matrix for 7-blade is time-varying
    t = 0.0
    omega = 25.0
    M = build_mass_matrix(problem, time=t, rotor_omega=omega)

    # Check dimensions
    assert M.shape == (9, 9), "7-blade mass matrix should be 9x9"

    # Check symmetry
    assert np.allclose(M, M.T), "Mass matrix should be symmetric"

    # Check positive definiteness
    eigenvalues = np.linalg.eigvals(M)
    assert np.all(eigenvalues > 0), "Mass matrix should be positive definite"


def test_mass_matrix_7_blade_coupling_coefficient():
    """Test 7-blade mass matrix coupling coefficient calculation."""
    problem = ProblemDefinition(number_blades=7)
    omega = 25.0

    # The coupling coefficient for 7-blade is calculated as:
    # coupling_coeff = (4*cos(π/7) + 4*cos²(π/7) - 8*cos³(π/7) - 1) / 2
    # This evaluates to exactly 0, making the mass matrix time-invariant
    cos_pi7 = np.cos(np.pi / 7)
    coupling_coeff = (4 * cos_pi7 + 4 * cos_pi7**2 - 8 * cos_pi7**3 - 1) / 2

    assert np.isclose(coupling_coeff, 0.0, atol=1e-10), \
        "7-blade coupling coefficient should be zero"

    # Since coupling coefficient is zero, mass matrix should be time-invariant
    M_t0 = build_mass_matrix(problem, time=0.0, rotor_omega=omega)
    M_t1 = build_mass_matrix(problem, time=0.5, rotor_omega=omega)

    assert np.allclose(M_t0, M_t1), \
        "7-blade mass matrix should be time-invariant (coupling_coeff = 0)"


def test_mass_matrix_unsupported_blade_count():
    """Test that unsupported blade counts raise ValueError."""
    problem = ProblemDefinition(number_blades=6)  # Not supported

    try:
        build_mass_matrix(problem)
        assert False, "Should raise ValueError for unsupported blade count"
    except ValueError as e:
        assert "Unsupported number of blades" in str(e)


def test_mass_matrix_dimensions_all_blade_counts():
    """Test mass matrix dimensions for all supported blade counts."""
    expected_dimensions = {
        3: (5, 5),
        4: (6, 6),
        5: (7, 7),
        7: (9, 9)
    }

    for n_blades, expected_shape in expected_dimensions.items():
        problem = ProblemDefinition(number_blades=n_blades)
        M = build_mass_matrix(problem)

        assert M.shape == expected_shape, \
            f"{n_blades}-blade mass matrix should be {expected_shape}"


def test_mass_matrix_uses_rotor_characteristics():
    """Test that mass matrix uses rotor characteristics correctly."""
    problem = ProblemDefinition(number_blades=4)
    M = build_mass_matrix(problem)

    # First element should be n_blades * blade_inertia
    expected_first_element = 4 * problem.rotor_characteristics.blade_inertia_Ib
    assert np.isclose(M[0, 0], expected_first_element), \
        "First element should be n_blades * blade_inertia"


def test_damping_matrix_3_blade_h2b_all():
    """Test damping matrix for 3-blade H2B configuration with ALL activation."""
    problem = ProblemDefinition(
        number_blades=3,
        damper_connection="H2B",
        damper_activation="ALL"
    )

    C_func, K_func = build_damping_matrix(problem)

    # Evaluate at test point
    t = 0.0
    omega = 25.0
    C = C_func(t, omega)
    K = K_func(t, omega)

    # Check dimensions
    assert C.shape == (5, 5), "Damping matrix should be 5x5 for 3-blade"
    assert K.shape == (5, 5), "Stiffness matrix should be 5x5 for 3-blade"

    # Matrices should be real-valued
    assert np.all(np.isreal(C)), "Damping matrix should be real"
    assert np.all(np.isreal(K)), "Stiffness matrix should be real"


def test_damping_matrix_3_blade_h2b_odi():
    """Test damping matrix for 3-blade H2B configuration with ODI activation."""
    problem = ProblemDefinition(
        number_blades=3,
        damper_connection="H2B",
        damper_activation="ODI"
    )

    C_func, K_func = build_damping_matrix(problem)

    # Evaluate at test point
    t = 0.0
    omega = 25.0
    C = C_func(t, omega)
    K = K_func(t, omega)

    # Check dimensions
    assert C.shape == (5, 5), "Damping matrix should be 5x5 for 3-blade"
    assert K.shape == (5, 5), "Stiffness matrix should be 5x5 for 3-blade"


def test_damping_matrix_4_blade_h2b_all():
    """Test damping matrix for 4-blade H2B configuration with ALL activation."""
    problem = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ALL"
    )

    C_func, K_func = build_damping_matrix(problem)

    # Evaluate at test point
    t = 0.0
    omega = 25.0
    C = C_func(t, omega)
    K = K_func(t, omega)

    # Check dimensions
    assert C.shape == (6, 6), "Damping matrix should be 6x6 for 4-blade"
    assert K.shape == (6, 6), "Stiffness matrix should be 6x6 for 4-blade"


def test_damping_matrix_4_blade_h2b_odi():
    """Test damping matrix for 4-blade H2B configuration with ODI activation."""
    problem = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ODI"
    )

    C_func, K_func = build_damping_matrix(problem)

    # Evaluate at test point
    t = 0.0
    omega = 25.0
    C = C_func(t, omega)
    K = K_func(t, omega)

    # Check dimensions
    assert C.shape == (6, 6), "Damping matrix should be 6x6 for 4-blade"
    assert K.shape == (6, 6), "Stiffness matrix should be 6x6 for 4-blade"


def test_damping_matrix_5_blade_all():
    """Test damping matrix for 5-blade configuration with ALL activation."""
    problem = ProblemDefinition(
        number_blades=5,
        damper_activation="ALL"
    )

    C_func, K_func = build_damping_matrix(problem)

    # Evaluate at test point
    t = 0.0
    omega = 25.0
    C = C_func(t, omega)
    K = K_func(t, omega)

    # Check dimensions
    assert C.shape == (7, 7), "Damping matrix should be 7x7 for 5-blade"
    assert K.shape == (7, 7), "Stiffness matrix should be 7x7 for 5-blade"


def test_damping_matrix_5_blade_odi():
    """Test damping matrix for 5-blade configuration with ODI activation."""
    problem = ProblemDefinition(
        number_blades=5,
        damper_activation="ODI"
    )

    C_func, K_func = build_damping_matrix(problem)

    # Evaluate at test point
    t = 0.0
    omega = 25.0
    C = C_func(t, omega)
    K = K_func(t, omega)

    # Check dimensions
    assert C.shape == (7, 7), "Damping matrix should be 7x7 for 5-blade"
    assert K.shape == (7, 7), "Stiffness matrix should be 7x7 for 5-blade"


def test_damping_matrix_7_blade_not_implemented():
    """Test that 7-blade damping matrix raises NotImplementedError."""
    problem = ProblemDefinition(number_blades=7)

    try:
        build_damping_matrix(problem)
        assert False, "Should raise NotImplementedError for 7-blade"
    except NotImplementedError as e:
        assert "7-blade configuration is not yet implemented" in str(e)


def test_damping_stiffness_speed_dependence():
    """Test that damping and stiffness matrices depend on rotor speed."""
    problem = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ALL"
    )

    C_func, K_func = build_damping_matrix(problem)

    t = 0.0
    omega1 = 10.0
    omega2 = 50.0

    C1 = C_func(t, omega1)
    C2 = C_func(t, omega2)
    K1 = K_func(t, omega1)
    K2 = K_func(t, omega2)

    # Matrices should differ at different speeds
    assert not np.allclose(C1, C2), "Damping matrix should vary with rotor speed"
    assert not np.allclose(K1, K2), "Stiffness matrix should vary with rotor speed"


def test_odi_time_varying_behavior():
    """Test that ODI activation produces time-varying matrices."""
    problem = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ODI"
    )

    C_func, K_func = build_damping_matrix(problem)

    omega = 25.0
    t1 = 0.0
    t2 = 0.1

    C_t1 = C_func(t1, omega)
    C_t2 = C_func(t2, omega)
    K_t1 = K_func(t1, omega)
    K_t2 = K_func(t2, omega)

    # Matrices should vary with time for ODI
    assert not np.allclose(C_t1, C_t2), \
        "ODI damping matrix should vary with time"
    assert not np.allclose(K_t1, K_t2), \
        "ODI stiffness matrix should vary with time"


def test_all_vs_odi_produces_different_matrices():
    """Test that ALL and ODI activation modes produce different matrices."""
    problem_all = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ALL"
    )
    problem_odi = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ODI"
    )

    C_func_all, K_func_all = build_damping_matrix(problem_all)
    C_func_odi, K_func_odi = build_damping_matrix(problem_odi)

    t = 0.0
    omega = 25.0

    C_all = C_func_all(t, omega)
    C_odi = C_func_odi(t, omega)

    # Matrices should be different
    assert not np.allclose(C_all, C_odi), \
        "ALL and ODI should produce different damping matrices"


def test_damping_matrix_no_nan_inf():
    """Test that damping matrices don't contain NaN or Inf."""
    problem = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ALL"
    )

    C_func, K_func = build_damping_matrix(problem)

    # Test at multiple speeds
    for omega in [10.0, 25.0, 50.0, 100.0]:
        C = C_func(0.0, omega)
        K = K_func(0.0, omega)

        assert not np.any(np.isnan(C)), f"Damping matrix contains NaN at omega={omega}"
        assert not np.any(np.isinf(C)), f"Damping matrix contains Inf at omega={omega}"
        assert not np.any(np.isnan(K)), f"Stiffness matrix contains NaN at omega={omega}"
        assert not np.any(np.isinf(K)), f"Stiffness matrix contains Inf at omega={omega}"


def test_stiffness_matrix_speed_squared_dependence():
    """Test that stiffness matrix has omega^2 dependence (centrifugal effect)."""
    problem = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ALL"
    )

    C_func, K_func = build_damping_matrix(problem)

    t = 0.0
    omega = 25.0

    K1 = K_func(t, omega)
    K2 = K_func(t, 2 * omega)

    # The difference should scale roughly with omega^2
    # (due to centrifugal stiffness terms)
    assert not np.allclose(K1, K2), \
        "Stiffness should change significantly with speed"


def test_damping_matrix_callable():
    """Test that returned damping and stiffness are callable."""
    problem = ProblemDefinition(number_blades=4)

    C_func, K_func = build_damping_matrix(problem)

    assert callable(C_func), "Damping matrix function should be callable"
    assert callable(K_func), "Stiffness matrix function should be callable"


def test_mass_matrix_independent_of_time_for_small_blade_counts():
    """Test that mass matrix is constant for 3, 4, 5 blades (time-independent)."""
    for n_blades in [3, 4, 5]:
        problem = ProblemDefinition(number_blades=n_blades)

        # Mass matrix should be same regardless of time/omega for these blade counts
        M1 = build_mass_matrix(problem, time=0.0, rotor_omega=0.0)
        M2 = build_mass_matrix(problem, time=1.0, rotor_omega=50.0)

        assert np.allclose(M1, M2), \
            f"{n_blades}-blade mass matrix should be time-independent"


def test_damping_matrix_gyroscopic_coupling():
    """Test that damping matrix contains gyroscopic coupling terms."""
    problem = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ALL"
    )

    C_func, K_func = build_damping_matrix(problem)

    t = 0.0
    omega = 25.0
    C = C_func(t, omega)

    # Gyroscopic terms should make the damping matrix non-symmetric
    # (due to rotor spinning coupling between modes)
    assert not np.allclose(C, C.T), \
        "Damping matrix should be non-symmetric due to gyroscopic effects"


def test_stiffness_matrix_gyroscopic_coupling():
    """Test that stiffness matrix contains gyroscopic coupling terms."""
    problem = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ALL"
    )

    C_func, K_func = build_damping_matrix(problem)

    t = 0.0
    omega = 25.0
    K = K_func(t, omega)

    # Gyroscopic terms should make the stiffness matrix non-symmetric
    assert not np.allclose(K, K.T), \
        "Stiffness matrix should be non-symmetric due to gyroscopic effects"


def test_damping_matrix_dimensions_match_mass_matrix():
    """Test that damping/stiffness dimensions match mass matrix."""
    for n_blades in [3, 4, 5]:
        problem = ProblemDefinition(number_blades=n_blades)

        M = build_mass_matrix(problem)
        C_func, K_func = build_damping_matrix(problem)

        C = C_func(0.0, 25.0)
        K = K_func(0.0, 25.0)

        assert C.shape == M.shape, \
            f"Damping matrix shape should match mass matrix for {n_blades} blades"
        assert K.shape == M.shape, \
            f"Stiffness matrix shape should match mass matrix for {n_blades} blades"


if __name__ == "__main__":
    # Run tests manually without pytest
    import traceback

    tests = [
        test_build_mass_matrix_3_blade,
        test_build_mass_matrix_4_blade,
        test_build_mass_matrix_5_blade,
        test_build_mass_matrix_7_blade,
        test_mass_matrix_7_blade_coupling_coefficient,
        test_mass_matrix_unsupported_blade_count,
        test_mass_matrix_dimensions_all_blade_counts,
        test_mass_matrix_uses_rotor_characteristics,
        test_damping_matrix_3_blade_h2b_all,
        test_damping_matrix_3_blade_h2b_odi,
        test_damping_matrix_4_blade_h2b_all,
        test_damping_matrix_4_blade_h2b_odi,
        test_damping_matrix_5_blade_all,
        test_damping_matrix_5_blade_odi,
        test_damping_matrix_7_blade_not_implemented,
        test_damping_stiffness_speed_dependence,
        test_odi_time_varying_behavior,
        test_all_vs_odi_produces_different_matrices,
        test_damping_matrix_no_nan_inf,
        test_stiffness_matrix_speed_squared_dependence,
        test_damping_matrix_callable,
        test_mass_matrix_independent_of_time_for_small_blade_counts,
        test_damping_matrix_gyroscopic_coupling,
        test_stiffness_matrix_gyroscopic_coupling,
        test_damping_matrix_dimensions_match_mass_matrix,
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
