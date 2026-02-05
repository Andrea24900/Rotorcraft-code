"""Unit tests for problem_definition module.

Tests the ProblemDefinition and RotorCharacteristics configuration classes.
Validates default values, H2B/B2B configurations, and geometric calculations.

Run with:
    pytest tests/test_problem_definition.py -v
    python tests/test_problem_definition.py
"""

import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

import numpy as np
from config.problem_definition import ProblemDefinition, RotorCharacteristics, create_default_problem


def test_rotor_characteristics_default_values():
    """Test that RotorCharacteristics default values are set correctly."""
    rotor_char = RotorCharacteristics()

    assert rotor_char.blade_mass_Mb == 94.9, "Default blade mass should be 94.9 kg"
    assert rotor_char.blade_static_Sb == 289.1, "Default blade static moment should be 289.1 kg*m"
    assert rotor_char.blade_inertia_Ib == 1085.0, "Default blade inertia should be 1085.0 kg*m^2"
    assert rotor_char.hub_mass_Mx == 8027.0, "Default hub mass X should be 8027.0 kg"
    assert rotor_char.hub_mass_My == 3284.0, "Default hub mass Y should be 3284.0 kg"
    assert rotor_char.lag_hinge_offset_e == 0.3048, "Default lag hinge offset should be 0.3048 m"


def test_rotor_characteristics_custom_values():
    """Test RotorCharacteristics custom initialization."""
    rotor_char = RotorCharacteristics(
        blade_mass_Mb=100.0,
        hub_mass_Mx=8000.0
    )

    assert rotor_char.blade_mass_Mb == 100.0, "Custom blade mass should be set"
    assert rotor_char.hub_mass_Mx == 8000.0, "Custom hub mass should be set"
    # Other values should still be defaults
    assert rotor_char.blade_static_Sb == 289.1, "Non-specified values should remain default"


def test_problem_definition_default_creation():
    """Test default problem creation."""
    problem = create_default_problem()

    assert problem.rotor_characteristics.number_blades == 4, "Default should be 4 blades"
    assert problem.rotor_characteristics.damper_connection == "H2B", "Default damper connection should be H2B"
    assert problem.required_solver == "LTI", "Default solver should be LTI"
    assert problem.lower_rotor_RPM == 20.0, "Default lower RPM should be 20.0"
    assert problem.higher_rotor_RPM == 400.0, "Default higher RPM should be 400.0"
    assert problem.number_points == 200, "Default number of points should be 200"


def test_h2b_configuration():
    """Test Hub-to-Blade damper configuration."""
    problem = ProblemDefinition(rotor_characteristics=RotorCharacteristics(damper_connection="H2B"))

    assert problem.rotor_characteristics.nominal_damping_Cd == 4067.0, \
        "H2B damping should be 4067.0 Nms/rad"
    assert problem.rotor_characteristics.nominal_stiffness_Kd == 0.0, \
        "H2B stiffness should be 0.0 Nm/rad"


def test_b2b_configuration():
    """Test Blade-to-Blade damper configuration."""
    problem = ProblemDefinition(
        rotor_characteristics=RotorCharacteristics(damper_connection="B2B", number_blades=4)
    )

    # Check damper values
    assert problem.rotor_characteristics.nominal_damping_Cd == 36500.0, \
        "B2B damping should be 36500.0 Ns/m"

    # Check geometric calculations
    expected_delta_psi = np.pi / 2  # For 4 blades
    assert np.isclose(problem.rotor_characteristics.delta_psi, expected_delta_psi), \
        f"delta_psi should be approximately {expected_delta_psi}"

    assert problem.rotor_characteristics.phi > 0, "phi should be positive"
    assert problem.rotor_characteristics.hinge_to_hinge_c > 0, \
        "hinge-to-hinge distance should be positive"
    assert problem.rotor_characteristics.length_undef_l0 > 0, \
        "undeformed length should be positive"
    assert problem.rotor_characteristics.length_eq_lE > 0, \
        "equilibrium length should be positive"


def test_b2b_geometric_calculations():
    """Test B2B geometric parameter calculations for different blade counts."""
    for n_blades in [3, 4, 5, 7]:
        problem = ProblemDefinition(
            rotor_characteristics=RotorCharacteristics(damper_connection="B2B", number_blades=n_blades)
        )

        expected_delta_psi = 2 * np.pi / n_blades
        assert np.isclose(problem.rotor_characteristics.delta_psi, expected_delta_psi), \
            f"delta_psi for {n_blades} blades should be {expected_delta_psi}"

        # All geometric values should be calculated and positive
        assert problem.rotor_characteristics.phi > 0
        assert problem.rotor_characteristics.hinge_to_hinge_c > 0
        assert problem.rotor_characteristics.length_undef_l0 > 0


def test_number_blades_validation():
    """Test that valid blade numbers are accepted."""
    valid_blades = [3, 4, 5, 7]

    for n_blades in valid_blades:
        problem = ProblemDefinition(rotor_characteristics=RotorCharacteristics(number_blades=n_blades))
        assert problem.rotor_characteristics.number_blades == n_blades, \
            f"Should accept {n_blades} blades"


def test_solver_types():
    """Test different solver configurations."""
    solvers = ["LTI", "LTP", "HD"]

    for solver in solvers:
        problem = ProblemDefinition(required_solver=solver)
        assert problem.required_solver == solver, \
            f"Solver should be set to {solver}"


def test_damper_activation_modes():
    """Test different damper activation configurations."""
    activation_modes = ["ALL", "ODI", "2DI ADJ", "2DI OPP", "3DI"]

    for mode in activation_modes:
        problem = ProblemDefinition(rotor_characteristics=RotorCharacteristics(damper_activation=mode))
        assert problem.rotor_characteristics.damper_activation == mode, \
            f"Damper activation should be set to {mode}"


def test_rpm_range_configuration():
    """Test RPM range configuration."""
    problem = ProblemDefinition(
        lower_rotor_RPM=50.0,
        higher_rotor_RPM=500.0,
        number_points=100
    )

    assert problem.lower_rotor_RPM == 50.0, "Lower RPM should be 50.0"
    assert problem.higher_rotor_RPM == 500.0, "Higher RPM should be 500.0"
    assert problem.number_points == 100, "Number of points should be 100"


def test_harmonics_setting():
    """Test number of harmonics for HD analysis."""
    problem = ProblemDefinition(
        required_solver="HD",
        number_harmonics=3
    )

    assert problem.number_harmonics == 3, "Number of harmonics should be 3"


def test_continuation_parameters():
    """Test continuation method parameters."""
    problem = ProblemDefinition(
        continuation="YES",
        continuation_tolerance=1e-3,
        continuation_max_iter=200,
        step_h=1e-5
    )

    assert problem.continuation == "YES", "Continuation should be enabled"
    assert problem.continuation_tolerance == 1e-3, "Tolerance should be 1e-3"
    assert problem.continuation_max_iter == 200, "Max iterations should be 200"
    assert problem.step_h == 1e-5, "Step size should be 1e-5"


def test_rotor_characteristics_embedded():
    """Test that rotor characteristics are properly embedded in problem definition."""
    problem = ProblemDefinition()

    assert isinstance(problem.rotor_characteristics, RotorCharacteristics), \
        "Should contain RotorCharacteristics instance"
    assert problem.rotor_characteristics.blade_mass_Mb > 0, \
        "Blade mass should be positive"
    assert problem.rotor_characteristics.hub_mass_Mx > 0, \
        "Hub mass should be positive"


if __name__ == "__main__":
    # Run tests manually without pytest
    import traceback

    tests = [
        test_rotor_characteristics_default_values,
        test_rotor_characteristics_custom_values,
        test_problem_definition_default_creation,
        test_h2b_configuration,
        test_b2b_configuration,
        test_b2b_geometric_calculations,
        test_number_blades_validation,
        test_solver_types,
        test_damper_activation_modes,
        test_rpm_range_configuration,
        test_harmonics_setting,
        test_continuation_parameters,
        test_rotor_characteristics_embedded,
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
