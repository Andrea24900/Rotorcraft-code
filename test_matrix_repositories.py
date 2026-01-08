"""Test script for matrix repository functions.

This script tests the mass, damping, and stiffness matrix functions
for different rotor configurations.
"""

import numpy as np
from src.config.problem_definition import ProblemDefinition
from src.utils.matrix_repositories import build_mass_matrix, build_damping_matrix


def print_matrix_info(matrix, name):
    """Print information about a matrix."""
    print(f"\n{name}:")
    print(f"  Shape: {matrix.shape}")
    print(f"  Symmetric: {np.allclose(matrix, matrix.T)}")
    print(f"  Max value: {np.max(np.abs(matrix)):.2e}")
    print(f"  Min value: {np.min(np.abs(matrix[matrix != 0])):.2e}")
    print(f"  Sample values (first 3x3):\n{matrix[:3, :3]}")


def test_3_blade_configuration():
    """Test 3-blade rotor configuration."""
    print("\n" + "="*70)
    print("Testing 3-BLADE CONFIGURATION")
    print("="*70)

    # Test with ALL activation (H2B)
    print("\n--- 3-Blade, H2B/ALL ---")
    problem = ProblemDefinition(
        number_blades=3,
        damper_connection="H2B",
        damper_activation="ALL"
    )

    # Test mass matrix
    time = 0.0
    rotor_omega = 25.0  # rad/s (about 240 RPM)
    M = build_mass_matrix(problem, time, rotor_omega)
    print_matrix_info(M, "Mass Matrix M")

    # Test damping and stiffness matrices
    C_func, K_func = build_damping_matrix(problem)
    C = C_func(time, rotor_omega)
    K = K_func(time, rotor_omega)
    print_matrix_info(C, "Damping Matrix C")
    print_matrix_info(K, "Stiffness Matrix K")

    # Test with ODI activation (time-varying)
    print("\n--- 3-Blade, H2B/ODI ---")
    problem.damper_activation = "ODI"
    problem._initialize_rotor_characteristics()

    C_func, K_func = build_damping_matrix(problem)

    # Test at different time points
    times = [0.0, 0.1, 0.2]
    for t in times:
        C = C_func(t, rotor_omega)
        K = K_func(t, rotor_omega)
        print(f"\nTime t={t:.1f}s:")
        print(f"  C[0,1] = {C[0,1]:.2f}, C[1,2] = {C[1,2]:.2f}")
        print(f"  K[0,1] = {K[0,1]:.2f}, K[1,2] = {K[1,2]:.2f}")


def test_4_blade_configuration():
    """Test 4-blade rotor configuration."""
    print("\n" + "="*70)
    print("Testing 4-BLADE CONFIGURATION")
    print("="*70)

    # Test with ALL activation
    print("\n--- 4-Blade, H2B/ALL ---")
    problem = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ALL"
    )

    time = 0.0
    rotor_omega = 30.0  # rad/s
    M = build_mass_matrix(problem, time, rotor_omega)
    print_matrix_info(M, "Mass Matrix M")

    C_func, K_func = build_damping_matrix(problem)
    C = C_func(time, rotor_omega)
    K = K_func(time, rotor_omega)
    print_matrix_info(C, "Damping Matrix C")
    print_matrix_info(K, "Stiffness Matrix K")

    # Test with ODI activation
    print("\n--- 4-Blade, H2B/ODI ---")
    problem.damper_activation = "ODI"
    problem._initialize_rotor_characteristics()

    C_func, K_func = build_damping_matrix(problem)
    C = C_func(0.0, rotor_omega)
    K = K_func(0.0, rotor_omega)
    print_matrix_info(C, "Damping Matrix C (t=0)")
    print_matrix_info(K, "Stiffness Matrix K (t=0)")


def test_5_blade_configuration():
    """Test 5-blade rotor configuration."""
    print("\n" + "="*70)
    print("Testing 5-BLADE CONFIGURATION")
    print("="*70)

    # Test with ALL activation
    print("\n--- 5-Blade, ALL ---")
    problem = ProblemDefinition(
        number_blades=5,
        damper_activation="ALL"
    )

    time = 0.0
    rotor_omega = 20.0  # rad/s
    M = build_mass_matrix(problem, time, rotor_omega)
    print_matrix_info(M, "Mass Matrix M")

    C_func, K_func = build_damping_matrix(problem)
    C = C_func(time, rotor_omega)
    K = K_func(time, rotor_omega)
    print_matrix_info(C, "Damping Matrix C")
    print_matrix_info(K, "Stiffness Matrix K")

    # Test with ODI activation
    print("\n--- 5-Blade, ODI ---")
    problem.damper_activation = "ODI"
    problem._initialize_rotor_characteristics()

    C_func, K_func = build_damping_matrix(problem)
    C = C_func(0.0, rotor_omega)
    K = K_func(0.0, rotor_omega)
    print_matrix_info(C, "Damping Matrix C (t=0)")
    print_matrix_info(K, "Stiffness Matrix K (t=0)")


def test_7_blade_mass_matrix():
    """Test 7-blade mass matrix (time-varying)."""
    print("\n" + "="*70)
    print("Testing 7-BLADE MASS MATRIX (Time-Varying)")
    print("="*70)

    problem = ProblemDefinition(number_blades=7)
    rotor_omega = 25.0  # rad/s

    # Test at different times to show time variation
    times = [0.0, 0.05, 0.1, 0.15]
    print("\nTime-varying elements:")
    for t in times:
        M = build_mass_matrix(problem, t, rotor_omega)
        print(f"\nTime t={t:.2f}s:")
        print(f"  M[1,2] (inertia coupling at 2*omega) = {M[1,2]:.4f}")
        print(f"  M[3,4] (inertia coupling at 4*omega) = {M[3,4]:.4f}")
        print(f"  M[5,6] (inertia coupling at 6*omega) = {M[5,6]:.4f}")


def test_rotor_speed_variation():
    """Test how matrices change with rotor speed."""
    print("\n" + "="*70)
    print("Testing ROTOR SPEED VARIATION (4-Blade, H2B/ALL)")
    print("="*70)

    problem = ProblemDefinition(
        number_blades=4,
        damper_connection="H2B",
        damper_activation="ALL"
    )

    C_func, K_func = build_damping_matrix(problem)

    # Test at different rotor speeds
    rotor_speeds = [10.0, 20.0, 30.0, 40.0]  # rad/s

    print("\nGyroscopic coupling (C[1,2]) vs rotor speed:")
    for omega in rotor_speeds:
        C = C_func(0.0, omega)
        rpm = omega * 60 / (2 * np.pi)
        print(f"  omega = {omega:4.1f} rad/s ({rpm:5.1f} RPM): C[1,2] = {C[1,2]:8.1f} Nms/rad")

    print("\nStiffness K[1,1] (centrifugal + nominal) vs rotor speed:")
    for omega in rotor_speeds:
        K = K_func(0.0, omega)
        rpm = omega * 60 / (2 * np.pi)
        print(f"  omega = {omega:4.1f} rad/s ({rpm:5.1f} RPM): K[1,1] = {K[1,1]:10.1f} Nm/rad")


def test_matrix_properties():
    """Test physical properties of the matrices."""
    print("\n" + "="*70)
    print("Testing MATRIX PHYSICAL PROPERTIES")
    print("="*70)

    problem = ProblemDefinition(number_blades=4)
    rotor_omega = 25.0

    # Mass matrix should be symmetric and positive definite
    M = build_mass_matrix(problem, 0.0, rotor_omega)
    print("\nMass Matrix Properties:")
    print(f"  Symmetric: {np.allclose(M, M.T)}")
    eigenvalues = np.linalg.eigvals(M)
    print(f"  Positive definite: {np.all(eigenvalues > 0)}")
    print(f"  Eigenvalues: {eigenvalues}")

    # Damping matrix should have positive elements on diagonal
    C_func, K_func = build_damping_matrix(problem)
    C = C_func(0.0, rotor_omega)
    print("\nDamping Matrix Properties:")
    print(f"  Diagonal elements > 0: {np.all(np.diag(C) > 0)}")
    print(f"  Diagonal: {np.diag(C)}")

    # Stiffness matrix properties
    K = K_func(0.0, rotor_omega)
    print("\nStiffness Matrix Properties:")
    print(f"  Diagonal elements > 0: {np.all(np.diag(K) > 0)}")
    K_eigenvalues = np.linalg.eigvals(K)
    print(f"  Positive definite: {np.all(K_eigenvalues > 0)}")


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print("MATRIX REPOSITORY TEST SUITE")
    print("="*70)

    try:
        test_3_blade_configuration()
        test_4_blade_configuration()
        test_5_blade_configuration()
        test_7_blade_mass_matrix()
        test_rotor_speed_variation()
        test_matrix_properties()

        print("\n" + "="*70)
        print("ALL TESTS COMPLETED SUCCESSFULLY!")
        print("="*70 + "\n")

    except Exception as e:
        print(f"\n\nERROR: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
