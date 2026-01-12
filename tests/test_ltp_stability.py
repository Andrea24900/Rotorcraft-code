"""Unit tests for stability_analysis.ltp_stability module.

Tests the LTPStability class and its methods including:
- ltp_single_point (single operating point with Floquet analysis)
- ltp_full_range (full RPM range analysis)
- Characteristic multipliers and exponents
- Input validation
- Error handling
- Integration with real rotor systems

Run with:
    python -m pytest tests/test_ltp_stability.py -v
    python tests/test_ltp_stability.py
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pytest
from src.stability_analysis.ltp_stability import LTPStability
from src.rotor_build import RotorBuild


# ============================================================================
# Test LTPStability Initialization
# ============================================================================

def test_ltp_initialization():
    """Test LTPStability initialization."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    assert ltp.rotor_build is rotor
    assert isinstance(ltp.modal_solution, list)
    assert len(ltp.modal_solution) == 0


def test_ltp_inherits_from_base():
    """Test that LTPStability properly inherits from StabilityAnalysis."""
    from src.stability_analysis.base import StabilityAnalysis

    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    assert isinstance(ltp, StabilityAnalysis)
    # Should have inherited methods
    assert hasattr(ltp, 'assign_range_OMEGA')
    assert hasattr(ltp, 'assign_period_T')
    assert hasattr(ltp, 'monodromy_computer')


# ============================================================================
# Test ltp_single_point - Basic Functionality
# ============================================================================

def test_ltp_single_point_basic():
    """Test ltp_single_point with valid rotor system."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    OMEGA = 25.0  # rad/s
    T = 2 * np.pi / OMEGA
    char_sol = ltp.ltp_single_point(OMEGA, T)

    # Check return type
    assert isinstance(char_sol, dict)
    assert 'char_multipliers' in char_sol
    assert 'real_char_exp' in char_sol
    assert 'imag_char_exp' in char_sol

    # Check array types
    assert isinstance(char_sol['char_multipliers'], np.ndarray)
    assert isinstance(char_sol['real_char_exp'], np.ndarray)
    assert isinstance(char_sol['imag_char_exp'], np.ndarray)

    # Check dimensions
    n_states = 2 * rotor.mass_matrix.shape[0]
    assert char_sol['char_multipliers'].shape == (n_states,)
    assert char_sol['real_char_exp'].shape == (n_states,)
    assert char_sol['imag_char_exp'].shape == (n_states,)


def test_ltp_single_point_multipliers_complex():
    """Test that characteristic multipliers can be complex."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    OMEGA = 50.0
    T = 2 * np.pi / OMEGA
    char_sol = ltp.ltp_single_point(OMEGA, T)

    multipliers = char_sol['char_multipliers']

    # Multipliers should generally be complex for oscillatory systems
    has_complex = np.any(np.imag(multipliers) != 0)
    assert has_complex or np.all(np.abs(multipliers) > 0), \
        "Characteristic multipliers should exist"


def test_ltp_single_point_exponents_real():
    """Test that characteristic exponents have real damping component."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    OMEGA = 100.0
    T = 2 * np.pi / OMEGA
    char_sol = ltp.ltp_single_point(OMEGA, T)

    real_exp = char_sol['real_char_exp']
    imag_exp = char_sol['imag_char_exp']

    # Real exponents should be real-valued
    assert np.all(np.isreal(real_exp))
    assert np.all(np.isreal(imag_exp))


def test_ltp_single_point_different_speeds():
    """Test ltp_single_point at different rotor speeds."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    speeds = [10.0, 50.0, 100.0]  # rad/s
    exponent_sets = []

    for omega in speeds:
        T = 2 * np.pi / omega
        char_sol = ltp.ltp_single_point(omega, T)
        exponent_sets.append(char_sol['real_char_exp'])

    # Exponents should be different at different speeds
    for i in range(len(speeds) - 1):
        difference = np.linalg.norm(exponent_sets[i] - exponent_sets[i+1])
        assert difference > 1e-6, \
            f"Characteristic exponents should differ between {speeds[i]} and {speeds[i+1]} rad/s"


# ============================================================================
# Test ltp_single_point - Input Validation
# ============================================================================

def test_ltp_single_point_negative_omega():
    """Test that negative OMEGA raises ValueError."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    with pytest.raises(ValueError, match="OMEGA must be positive"):
        ltp.ltp_single_point(-10.0, 1.0)


def test_ltp_single_point_zero_omega():
    """Test that zero OMEGA raises ValueError."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    with pytest.raises(ValueError, match="OMEGA must be positive"):
        ltp.ltp_single_point(0.0, 1.0)


def test_ltp_single_point_negative_T():
    """Test that negative T raises ValueError."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    with pytest.raises(ValueError, match="T .* must be positive"):
        ltp.ltp_single_point(25.0, -1.0)


def test_ltp_single_point_zero_T():
    """Test that zero T raises ValueError."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    with pytest.raises(ValueError, match="T .* must be positive"):
        ltp.ltp_single_point(25.0, 0.0)


def test_ltp_single_point_non_numeric_omega():
    """Test that non-numeric OMEGA raises TypeError."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    T = 0.1

    with pytest.raises(TypeError, match="OMEGA must be numeric"):
        ltp.ltp_single_point("not_a_number", T)

    with pytest.raises(TypeError, match="OMEGA must be numeric"):
        ltp.ltp_single_point(None, T)


def test_ltp_single_point_non_numeric_T():
    """Test that non-numeric T raises TypeError."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    OMEGA = 25.0

    with pytest.raises(TypeError, match="T must be numeric"):
        ltp.ltp_single_point(OMEGA, "not_a_number")

    with pytest.raises(TypeError, match="T must be numeric"):
        ltp.ltp_single_point(OMEGA, None)


def test_ltp_single_point_numpy_scalars():
    """Test that numpy scalar types are accepted."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    OMEGA = np.float64(25.0)
    T = np.float32(2 * np.pi / 25.0)

    # Should not raise exceptions
    char_sol = ltp.ltp_single_point(OMEGA, T)
    assert isinstance(char_sol, dict)


def test_ltp_single_point_inconsistent_omega_T():
    """Test warning when T doesn't match 2π/OMEGA."""
    import warnings

    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    OMEGA = 25.0
    T_wrong = 1.0  # Should be ~0.251 for OMEGA=25

    # Should issue a warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        char_sol = ltp.ltp_single_point(OMEGA, T_wrong)

        # Check that a warning was issued
        assert len(w) >= 1
        assert any("does not match expected value" in str(warning.message) for warning in w)


# ============================================================================
# Test ltp_full_range - Basic Functionality
# ============================================================================

def test_ltp_full_range_basic():
    """Test ltp_full_range with default rotor configuration."""
    rotor = RotorBuild.build_all()
    rotor.problem.required_solver = "LTP"
    rotor.problem.number_points = 3  # Keep small for speed

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    # Check that modal_solution is populated
    assert len(ltp.modal_solution) == 3
    assert len(ltp.modal_solution) > 0


def test_ltp_full_range_modal_solution_structure():
    """Test that modal_solution has correct structure."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 1000.0
    rotor.problem.higher_rotor_RPM = 2000.0
    rotor.problem.number_points = 3

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    # Check each modal solution
    for i, modal_sol in enumerate(ltp.modal_solution):
        # Check required attributes
        assert hasattr(modal_sol, 'OMEGA')
        assert hasattr(modal_sol, 'OMEGA_RPM')
        assert hasattr(modal_sol, 'T')
        assert hasattr(modal_sol, 'damping')
        assert hasattr(modal_sol, 'frequency')
        assert hasattr(modal_sol, 'char_solution')

        # Check OMEGA and T are positive
        assert modal_sol.OMEGA > 0
        assert modal_sol.T > 0

        # Check arrays are populated
        assert len(modal_sol.damping) > 0
        assert len(modal_sol.frequency) > 0
        assert len(modal_sol.damping) == len(modal_sol.frequency)

        # Check char_solution structure
        assert 'char_multipliers' in modal_sol.char_solution
        assert 'real_char_exp' in modal_sol.char_solution
        assert 'imag_char_exp' in modal_sol.char_solution


def test_ltp_full_range_rpm_values():
    """Test that RPM values are correctly set."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 1000.0
    rotor.problem.higher_rotor_RPM = 2000.0
    rotor.problem.number_points = 5

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    # Extract RPM values
    rpm_values = [sol.OMEGA_RPM for sol in ltp.modal_solution]

    # Check range (default is R2L, so descending)
    if rotor.problem.solution_direction == "R2L":
        assert np.isclose(rpm_values[0], 2000.0, atol=1.0)
        assert np.isclose(rpm_values[-1], 1000.0, atol=1.0)
    else:  # L2R
        assert np.isclose(rpm_values[0], 1000.0, atol=1.0)
        assert np.isclose(rpm_values[-1], 2000.0, atol=1.0)


def test_ltp_full_range_period_calculation():
    """Test that periods are correctly calculated."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 3

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    # Check that T = 2π/OMEGA at each point
    for modal_sol in ltp.modal_solution:
        expected_T = 2 * np.pi / modal_sol.OMEGA
        assert np.isclose(modal_sol.T, expected_T, rtol=1e-12), \
            f"Period should be 2π/OMEGA, got {modal_sol.T}, expected {expected_T}"


def test_ltp_full_range_exponent_relationship():
    """Test that damping and frequency are from characteristic exponents."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 2

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    for modal_sol in ltp.modal_solution:
        char_solution = modal_sol.char_solution
        damping = modal_sol.damping
        frequency = modal_sol.frequency

        # Check that damping = real characteristic exponents
        assert np.allclose(damping, char_solution['real_char_exp'])

        # Check that frequency = imaginary characteristic exponents
        assert np.allclose(frequency, char_solution['imag_char_exp'])


def test_ltp_full_range_method_chaining():
    """Test that ltp_full_range returns self for method chaining."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 2

    ltp = LTPStability(rotor)
    result = ltp.ltp_full_range()

    assert result is ltp, "Method should return self for chaining"


# ============================================================================
# Test ltp_full_range - Different Configurations
# ============================================================================

def test_ltp_full_range_single_point():
    """Test ltp_full_range with single operating point."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 1000.0
    rotor.problem.higher_rotor_RPM = 1000.0
    rotor.problem.number_points = 1

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    assert len(ltp.modal_solution) == 1
    assert np.isclose(ltp.modal_solution[0].OMEGA_RPM, 1000.0)


def test_ltp_full_range_l2r_direction():
    """Test ltp_full_range with L2R solution direction."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 1000.0
    rotor.problem.higher_rotor_RPM = 2000.0
    rotor.problem.number_points = 3
    rotor.problem.solution_direction = "L2R"

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    # Check ascending order
    rpm_values = [sol.OMEGA_RPM for sol in ltp.modal_solution]
    assert np.isclose(rpm_values[0], 1000.0, atol=1.0)
    assert np.isclose(rpm_values[-1], 2000.0, atol=1.0)


def test_ltp_full_range_r2l_direction():
    """Test ltp_full_range with R2L solution direction."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 1000.0
    rotor.problem.higher_rotor_RPM = 2000.0
    rotor.problem.number_points = 3
    rotor.problem.solution_direction = "R2L"

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    # Check descending order
    rpm_values = [sol.OMEGA_RPM for sol in ltp.modal_solution]
    assert np.isclose(rpm_values[0], 2000.0, atol=1.0)
    assert np.isclose(rpm_values[-1], 1000.0, atol=1.0)


# ============================================================================
# Test Floquet Theory Properties
# ============================================================================

def test_characteristic_multiplier_magnitude():
    """Test that characteristic multipliers have expected magnitude properties."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 2

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    for modal_sol in ltp.modal_solution:
        multipliers = modal_sol.char_solution['char_multipliers']

        # All multipliers should have positive magnitude
        magnitudes = np.abs(multipliers)
        assert np.all(magnitudes > 0), "All characteristic multipliers should have positive magnitude"


def test_damping_from_multipliers():
    """Test relationship between multipliers and damping."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 2

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    for modal_sol in ltp.modal_solution:
        multipliers = modal_sol.char_solution['char_multipliers']
        damping = modal_sol.damping
        T = modal_sol.T

        # Verify relationship: σ = (1/2T) * ln(|μ|²)
        expected_damping = (1 / (2 * T)) * np.log(np.abs(multipliers)**2)

        assert np.allclose(damping, expected_damping, rtol=1e-10), \
            "Damping should match formula: σ = (1/2T) * ln(|μ|²)"


def test_frequency_from_multipliers():
    """Test relationship between multipliers and frequency."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 2

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    for modal_sol in ltp.modal_solution:
        multipliers = modal_sol.char_solution['char_multipliers']
        frequency = modal_sol.frequency
        T = modal_sol.T

        # Verify relationship: ω = (1/T) * atan2(Im(μ), Re(μ))
        expected_frequency = (1 / T) * np.arctan2(
            np.imag(multipliers),
            np.real(multipliers)
        )

        assert np.allclose(frequency, expected_frequency, rtol=1e-10), \
            "Frequency should match formula: ω = (1/T) * atan2(Im(μ), Re(μ))"


def test_stability_interpretation():
    """Test interpretation of stability from characteristic exponents."""
    rotor = RotorBuild.build_all()
    rotor.problem.number_points = 2

    ltp = LTPStability(rotor)
    ltp = ltp.ltp_full_range()

    for modal_sol in ltp.modal_solution:
        damping = modal_sol.damping

        # All values should be real
        assert np.all(np.isreal(damping))

        # Damping < 0 indicates stable mode
        # Damping > 0 indicates unstable mode
        stable_modes = np.sum(damping < 0)
        unstable_modes = np.sum(damping > 0)
        neutral_modes = np.sum(np.abs(damping) < 1e-10)

        # At least some modes should exist
        assert stable_modes + unstable_modes + neutral_modes == len(damping)


# ============================================================================
# Test Integration with run_stability
# ============================================================================

def test_run_stability_ltp():
    """Test LTP stability through run_stability factory method."""
    from src.stability_analysis.base import StabilityAnalysis

    rotor = RotorBuild.build_all()
    rotor.problem.required_solver = "LTP"
    rotor.problem.continuation = "NO"
    rotor.problem.number_points = 2

    result = StabilityAnalysis.run_stability(rotor)

    # Should return LTPStability instance
    assert isinstance(result, LTPStability)

    # Should have computed results
    assert len(result.modal_solution) == 2

    # Check results are populated
    for modal_sol in result.modal_solution:
        assert len(modal_sol.damping) > 0


# ============================================================================
# Test Edge Cases
# ============================================================================

def test_ltp_single_point_very_small_omega():
    """Test ltp_single_point with very small but positive OMEGA."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    OMEGA = 1e-2  # Very small
    T = 2 * np.pi / OMEGA

    # Should work with very small positive values
    char_sol = ltp.ltp_single_point(OMEGA, T)
    assert isinstance(char_sol, dict)
    assert 'char_multipliers' in char_sol


def test_ltp_single_point_very_large_omega():
    """Test ltp_single_point with very large OMEGA."""
    rotor = RotorBuild.build_all()
    ltp = LTPStability(rotor)

    OMEGA = 500.0  # Large value
    T = 2 * np.pi / OMEGA

    # Should work with large values
    char_sol = ltp.ltp_single_point(OMEGA, T)
    assert isinstance(char_sol, dict)
    assert 'char_multipliers' in char_sol


def test_consistency_single_vs_range():
    """Test that ltp_single_point matches ltp_full_range at same OMEGA."""
    rotor = RotorBuild.build_all()
    rotor.problem.lower_rotor_RPM = 1000.0
    rotor.problem.higher_rotor_RPM = 2000.0
    rotor.problem.number_points = 3

    ltp = LTPStability(rotor)

    # Get full range results
    ltp = ltp.ltp_full_range()

    # Compare with single point computation
    for modal_sol in ltp.modal_solution:
        single_result = ltp.ltp_single_point(modal_sol.OMEGA, modal_sol.T)

        # Characteristic multipliers should match
        assert np.allclose(
            single_result['char_multipliers'],
            modal_sol.char_solution['char_multipliers'],
            rtol=1e-10
        ), f"Single point and full range should give same multipliers at OMEGA={modal_sol.OMEGA}"

        # Exponents should match
        assert np.allclose(
            single_result['real_char_exp'],
            modal_sol.damping,
            rtol=1e-10
        ), "Real exponents should match"


# ============================================================================
# Test Main Block Execution
# ============================================================================

def test_main_block_execution():
    """Test that __main__ block executes without errors."""
    import subprocess

    result = subprocess.run(
        [sys.executable, "-m", "src.stability_analysis.ltp_stability"],
        capture_output=True,
        text=True,
        timeout=60
    )

    # Should complete successfully
    assert result.returncode == 0, f"Main block failed: {result.stderr}"

    # Should print expected output
    assert "Building rotor system" in result.stdout or "Building rotor system" in result.stderr
    assert "Running LTP stability analysis" in result.stdout or "Running LTP stability analysis" in result.stderr


# ============================================================================
# Run tests manually
# ============================================================================

if __name__ == "__main__":
    import traceback

    # List all test functions
    test_functions = [
        test_ltp_initialization,
        test_ltp_inherits_from_base,
        test_ltp_single_point_basic,
        test_ltp_single_point_multipliers_complex,
        test_ltp_single_point_exponents_real,
        test_ltp_single_point_different_speeds,
        test_ltp_single_point_negative_omega,
        test_ltp_single_point_zero_omega,
        test_ltp_single_point_negative_T,
        test_ltp_single_point_zero_T,
        test_ltp_single_point_non_numeric_omega,
        test_ltp_single_point_non_numeric_T,
        test_ltp_single_point_numpy_scalars,
        test_ltp_single_point_inconsistent_omega_T,
        test_ltp_full_range_basic,
        test_ltp_full_range_modal_solution_structure,
        test_ltp_full_range_rpm_values,
        test_ltp_full_range_period_calculation,
        test_ltp_full_range_exponent_relationship,
        test_ltp_full_range_method_chaining,
        test_ltp_full_range_single_point,
        test_ltp_full_range_l2r_direction,
        test_ltp_full_range_r2l_direction,
        test_characteristic_multiplier_magnitude,
        test_damping_from_multipliers,
        test_frequency_from_multipliers,
        test_stability_interpretation,
        test_run_stability_ltp,
        test_ltp_single_point_very_small_omega,
        test_ltp_single_point_very_large_omega,
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
