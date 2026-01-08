"""Matrix repository functions for mass, damping, and stiffness matrices.

Converted from MATLAB mass_matrix_repository.m and damping_matrix_repository.m.
"""

import numpy as np
from typing import Callable, Tuple
from ..config.problem_definition import ProblemDefinition


def build_mass_matrix(
    problem: ProblemDefinition,
    time: float = 0.0,
    rotor_omega: float = 0.0
) -> np.ndarray:
    """Build mass matrix for the rotor system.

    Args:
        problem: Problem definition containing rotor characteristics
        time: Current time (required for 7-blade configuration)
        rotor_omega: Rotor angular velocity (required for 7-blade configuration)

    Returns:
        Mass matrix M

    Note:
        Matrix dimensions depend on number of blades:
        - 3 blades: 5x5 matrix
        - 4 blades: 6x6 matrix
        - 5 blades: 7x7 matrix
        - 7 blades: 9x9 matrix (time-varying)
    """
    num_blades = problem.number_blades
    rotor_chars = problem.rotor_characteristics

    # Extract rotor characteristics for clarity
    blade_inertia = rotor_chars.blade_inertia_Ib
    blade_mass = rotor_chars.blade_mass_Mb
    blade_static_moment = rotor_chars.blade_static_Sb
    hub_mass_x = rotor_chars.hub_mass_Mx
    hub_mass_y = rotor_chars.hub_mass_My

    # Total hub mass including blades
    total_hub_mass_x = hub_mass_x + blade_mass * num_blades
    total_hub_mass_y = hub_mass_y + blade_mass * num_blades

    if num_blades == 3:
        # 5x5 mass matrix for 3-blade rotor
        mass_matrix = np.array([
            [3 * blade_inertia, 0, 0, 0, 0],
            [0, 1.5 * blade_inertia, 0, 0, -1.5 * blade_static_moment],
            [0, 0, 1.5 * blade_inertia, 1.5 * blade_static_moment, 0],
            [0, 0, 1.5 * blade_static_moment, total_hub_mass_x, 0],
            [0, -1.5 * blade_static_moment, 0, 0, total_hub_mass_y]
        ])

    elif num_blades == 4:
        # 6x6 mass matrix for 4-blade rotor
        mass_matrix = np.array([
            [4 * blade_inertia, 0, 0, 0, 0, 0],
            [0, 2 * blade_inertia, 0, 0, 0, -2 * blade_static_moment],
            [0, 0, 2 * blade_inertia, 0, 2 * blade_static_moment, 0],
            [0, 0, 0, 4 * blade_inertia, 0, 0],
            [0, 0, 2 * blade_static_moment, 0, total_hub_mass_x, 0],
            [0, -2 * blade_static_moment, 0, 0, 0, total_hub_mass_y]
        ])

    elif num_blades == 5:
        # 7x7 mass matrix for 5-blade rotor
        mass_matrix = np.array([
            [5 * blade_inertia, 0, 0, 0, 0, 0, 0],
            [0, 2.5 * blade_inertia, 0, 0, 0, 0, -2.5 * blade_static_moment],
            [0, 0, 2.5 * blade_inertia, 0, 0, 2.5 * blade_static_moment, 0],
            [0, 0, 0, 2.5 * blade_inertia, 0, 0, 0],
            [0, 0, 0, 0, 2.5 * blade_inertia, 0, 0],
            [0, 0, 2.5 * blade_static_moment, 0, 0, total_hub_mass_x, 0],
            [0, -2.5 * blade_static_moment, 0, 0, 0, 0, total_hub_mass_y]
        ])

    elif num_blades == 7:
        # 9x9 mass matrix for 7-blade rotor (time-varying)
        # Time-varying coupling coefficient
        cos_pi7 = np.cos(np.pi / 7)
        coupling_coeff = (4 * cos_pi7 + 4 * cos_pi7**2 - 8 * cos_pi7**3 - 1) / 2

        # Time-varying terms
        sin_2wt = np.sin(2 * rotor_omega * time)
        sin_4wt = np.sin(4 * rotor_omega * time)
        sin_6wt = np.sin(6 * rotor_omega * time)

        # Off-diagonal inertia coupling terms
        inertia_coupling_2w = blade_inertia * sin_2wt * coupling_coeff
        inertia_coupling_4w = blade_inertia * sin_4wt * coupling_coeff
        inertia_coupling_6w = blade_inertia * sin_6wt * coupling_coeff

        # Off-diagonal static moment coupling terms
        static_coupling_2w = blade_static_moment * sin_2wt * coupling_coeff

        mass_matrix = np.array([
            [7 * blade_inertia, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 3.5 * blade_inertia, inertia_coupling_2w, 0, 0, 0, 0, static_coupling_2w, -3.5 * blade_static_moment],
            [0, inertia_coupling_2w, 3.5 * blade_inertia, 0, 0, 0, 0, 3.5 * blade_static_moment, -static_coupling_2w],
            [0, 0, 0, 3.5 * blade_inertia, inertia_coupling_4w, 0, 0, 0, 0],
            [0, 0, 0, inertia_coupling_4w, 3.5 * blade_inertia, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 3.5 * blade_inertia, inertia_coupling_6w, 0, 0],
            [0, 0, 0, 0, 0, inertia_coupling_6w, 3.5 * blade_inertia, 0, 0],
            [0, static_coupling_2w, 3.5 * blade_static_moment, 0, 0, 0, 0, 7 * blade_mass + hub_mass_x, 0],
            [0, -3.5 * blade_static_moment, -static_coupling_2w, 0, 0, 0, 0, 0, 7 * blade_mass + hub_mass_y]
        ])

    else:
        raise ValueError(f"Unsupported number of blades: {num_blades}. Supported values: 3, 4, 5, 7")

    return mass_matrix


def build_damping_matrix(
    problem: ProblemDefinition
) -> Tuple[Callable[[float, float], np.ndarray], Callable[[float, float], np.ndarray]]:
    """Build damping and stiffness matrices as function handles.

    Args:
        problem: Problem definition containing rotor characteristics

    Returns:
        Tuple of (damping_matrix_C, stiffness_matrix_K) as function handles
        Both functions take (t, rotor_omega) as arguments and return matrices

    Note:
        Converted from MATLAB damping_matrix_repository.m.
        Currently supports blade counts 3, 4, and 5.
        Case 7 (7-blade) is not yet implemented due to complexity.
    """
    num_blades = problem.number_blades
    rotor_chars = problem.rotor_characteristics

    # Extract commonly used parameters for clarity
    blade_inertia = rotor_chars.blade_inertia_Ib
    blade_static_moment = rotor_chars.blade_static_Sb
    lag_hinge_offset = rotor_chars.lag_hinge_offset_e
    nominal_damping = rotor_chars.nominal_damping_Cd
    nominal_stiffness = rotor_chars.nominal_stiffness_Kd
    hub_damping_x = rotor_chars.hub_damping_Cx
    hub_damping_y = rotor_chars.hub_damping_Cy
    hub_stiffness_x = rotor_chars.hub_stiffness_Kx
    hub_stiffness_y = rotor_chars.hub_stiffness_Ky

    if num_blades == 3:
        return _build_matrices_3_blade(
            blade_inertia, blade_static_moment, lag_hinge_offset,
            nominal_damping, nominal_stiffness,
            hub_damping_x, hub_damping_y, hub_stiffness_x, hub_stiffness_y,
            problem.damper_connection, problem.damper_activation
        )
    elif num_blades == 4:
        return _build_matrices_4_blade(
            blade_inertia, blade_static_moment, lag_hinge_offset,
            nominal_damping, nominal_stiffness,
            hub_damping_x, hub_damping_y, hub_stiffness_x, hub_stiffness_y,
            problem.damper_connection, problem.damper_activation
        )
    elif num_blades == 5:
        return _build_matrices_5_blade(
            blade_inertia, blade_static_moment, lag_hinge_offset,
            nominal_damping, nominal_stiffness,
            hub_damping_x, hub_damping_y, hub_stiffness_x, hub_stiffness_y,
            problem.damper_activation
        )
    elif num_blades == 7:
        raise NotImplementedError(
            "7-blade configuration is not yet implemented due to its complexity. "
            "Please use 3, 4, or 5 blades."
        )
    else:
        raise ValueError(f"Unsupported number of blades: {num_blades}. Supported values: 3, 4, 5")


def _build_matrices_3_blade(
    blade_inertia: float, blade_static_moment: float, lag_hinge_offset: float,
    nominal_damping: float, nominal_stiffness: float,
    hub_damping_x: float, hub_damping_y: float,
    hub_stiffness_x: float, hub_stiffness_y: float,
    damper_connection: str, damper_activation: str
) -> Tuple[Callable[[float, float], np.ndarray], Callable[[float, float], np.ndarray]]:
    """Build damping and stiffness matrices for 3-blade rotor."""

    if damper_connection == "H2B":
        if damper_activation == "ALL":
            def damping_matrix_C(time: float, rotor_omega: float) -> np.ndarray:
                return np.array([
                    [3 * nominal_damping, 0, 0, 0, 0],
                    [0, 1.5 * nominal_damping, 3 * blade_inertia * rotor_omega, 0, 0],
                    [0, -3 * blade_inertia * rotor_omega, 1.5 * nominal_damping, 0, 0],
                    [0, 0, 0, hub_damping_x, 0],
                    [0, 0, 0, 0, hub_damping_y]
                ])

            def stiffness_matrix_K(time: float, rotor_omega: float) -> np.ndarray:
                centrifugal_stiffness = blade_static_moment * lag_hinge_offset * rotor_omega**2
                gyroscopic_term = 1.5 * nominal_damping * rotor_omega
                dynamic_stiffness = 1.5 * nominal_stiffness - 1.5 * blade_inertia * rotor_omega**2 + \
                                    1.5 * centrifugal_stiffness

                return np.array([
                    [3 * centrifugal_stiffness + 3 * nominal_stiffness, 0, 0, 0, 0],
                    [0, dynamic_stiffness, gyroscopic_term, 0, 0],
                    [0, -gyroscopic_term, dynamic_stiffness, 0, 0],
                    [0, 0, 0, hub_stiffness_x, 0],
                    [0, 0, 0, 0, hub_stiffness_y]
                ])

        elif damper_activation == "ODI":
            def damping_matrix_C(time: float, rotor_omega: float) -> np.ndarray:
                wt = rotor_omega * time
                cos_wt = np.cos(wt)
                sin_wt = np.sin(wt)
                cos_2wt = np.cos(2 * wt)
                sin_2wt = np.sin(2 * wt)

                return np.array([
                    [2 * nominal_damping, -nominal_damping * cos_wt, -nominal_damping * sin_wt, 0, 0],
                    [-nominal_damping * cos_wt, -0.5 * nominal_damping * (cos_2wt - 2),
                     3 * blade_inertia * rotor_omega - 0.5 * nominal_damping * sin_2wt, 0, 0],
                    [-nominal_damping * sin_wt, -3 * blade_inertia * rotor_omega - 0.5 * nominal_damping * sin_2wt,
                     0.5 * nominal_damping * (cos_2wt + 2), 0, 0],
                    [0, 0, 0, hub_damping_x, 0],
                    [0, 0, 0, 0, hub_damping_y]
                ])

            def stiffness_matrix_K(time: float, rotor_omega: float) -> np.ndarray:
                wt = rotor_omega * time
                cos_wt = np.cos(wt)
                sin_wt = np.sin(wt)
                cos_2wt = np.cos(2 * wt)
                sin_2wt = np.sin(2 * wt)

                centrifugal_stiffness = blade_static_moment * lag_hinge_offset * rotor_omega**2
                base_stiffness = 1.5 * nominal_stiffness - 1.5 * blade_inertia * rotor_omega**2 + \
                                 1.5 * centrifugal_stiffness

                return np.array([
                    [3 * centrifugal_stiffness + 3 * nominal_stiffness,
                     nominal_damping * rotor_omega * sin_wt,
                     -nominal_damping * rotor_omega * cos_wt, 0, 0],
                    [0, base_stiffness + 0.5 * nominal_damping * rotor_omega * sin_2wt,
                     -0.5 * nominal_damping * rotor_omega * (cos_2wt - 2), 0, 0],
                    [0, -0.5 * nominal_damping * rotor_omega * (cos_2wt + 2),
                     base_stiffness - 0.5 * nominal_damping * rotor_omega * sin_2wt, 0, 0],
                    [0, 0, 0, hub_stiffness_x, 0],
                    [0, 0, 0, 0, hub_stiffness_y]
                ])
        else:
            raise ValueError(f"Unsupported damper_activation for 3-blade H2B: {damper_activation}")
    else:
        raise ValueError(f"Unsupported damper_connection for 3-blade: {damper_connection}")

    return damping_matrix_C, stiffness_matrix_K


def _build_matrices_4_blade(
    blade_inertia: float, blade_static_moment: float, lag_hinge_offset: float,
    nominal_damping: float, nominal_stiffness: float,
    hub_damping_x: float, hub_damping_y: float,
    hub_stiffness_x: float, hub_stiffness_y: float,
    damper_connection: str, damper_activation: str
) -> Tuple[Callable[[float, float], np.ndarray], Callable[[float, float], np.ndarray]]:
    """Build damping and stiffness matrices for 4-blade rotor."""

    if damper_connection == "H2B":
        if damper_activation == "ALL":
            def damping_matrix_C(time: float, rotor_omega: float) -> np.ndarray:
                return np.array([
                    [4 * nominal_damping, 0, 0, 0, 0, 0],
                    [0, 2 * nominal_damping, 4 * blade_inertia * rotor_omega, 0, 0, 0],
                    [0, -4 * blade_inertia * rotor_omega, 2 * nominal_damping, 0, 0, 0],
                    [0, 0, 0, 4 * nominal_damping, 0, 0],
                    [0, 0, 0, 0, hub_damping_x, 0],
                    [0, 0, 0, 0, 0, hub_damping_y]
                ])

            def stiffness_matrix_K(time: float, rotor_omega: float) -> np.ndarray:
                centrifugal_stiffness = blade_static_moment * lag_hinge_offset * rotor_omega**2
                gyroscopic_term = 2 * nominal_damping * rotor_omega
                dynamic_stiffness = 2 * nominal_stiffness - 2 * blade_inertia * rotor_omega**2 + \
                                    2 * centrifugal_stiffness

                return np.array([
                    [4 * centrifugal_stiffness + 4 * nominal_stiffness, 0, 0, 0, 0, 0],
                    [0, dynamic_stiffness, gyroscopic_term, 0, 0, 0],
                    [0, -gyroscopic_term, dynamic_stiffness, 0, 0, 0],
                    [0, 0, 0, 4 * centrifugal_stiffness + 4 * nominal_stiffness, 0, 0],
                    [0, 0, 0, 0, hub_stiffness_x, 0],
                    [0, 0, 0, 0, 0, hub_stiffness_y]
                ])

        elif damper_activation == "ODI":
            def damping_matrix_C(time: float, rotor_omega: float) -> np.ndarray:
                wt = rotor_omega * time
                cos_wt = np.cos(wt)
                sin_wt = np.sin(wt)
                cos_2wt = np.cos(2 * wt)
                sin_2wt = np.sin(2 * wt)

                return np.array([
                    [3 * nominal_damping, -nominal_damping * cos_wt, -nominal_damping * sin_wt,
                     nominal_damping, 0, 0],
                    [-nominal_damping * cos_wt, nominal_damping * (sin_wt**2 + 1),
                     4 * blade_inertia * rotor_omega - 0.5 * nominal_damping * sin_2wt,
                     nominal_damping * cos_wt, 0, 0],
                    [-nominal_damping * sin_wt, -4 * blade_inertia * rotor_omega - 0.5 * nominal_damping * sin_2wt,
                     0.5 * nominal_damping * (cos_2wt + 3), nominal_damping * sin_wt, 0, 0],
                    [nominal_damping, nominal_damping * cos_wt, nominal_damping * sin_wt,
                     3 * nominal_damping, 0, 0],
                    [0, 0, 0, 0, hub_damping_x, 0],
                    [0, 0, 0, 0, 0, hub_damping_y]
                ])

            def stiffness_matrix_K(time: float, rotor_omega: float) -> np.ndarray:
                wt = rotor_omega * time
                cos_wt = np.cos(wt)
                sin_wt = np.sin(wt)
                cos_2wt = np.cos(2 * wt)
                sin_2wt = np.sin(2 * wt)

                centrifugal_stiffness = blade_static_moment * lag_hinge_offset * rotor_omega**2
                base_stiffness = 2 * nominal_stiffness - 2 * blade_inertia * rotor_omega**2 + \
                                 2 * centrifugal_stiffness

                return np.array([
                    [4 * centrifugal_stiffness + 4 * nominal_stiffness,
                     nominal_damping * rotor_omega * sin_wt,
                     -nominal_damping * rotor_omega * cos_wt, 0, 0, 0],
                    [0, base_stiffness + 0.5 * nominal_damping * rotor_omega * sin_2wt,
                     -0.5 * nominal_damping * rotor_omega * (cos_2wt - 3), 0, 0, 0],
                    [0, -0.5 * nominal_damping * rotor_omega * (cos_2wt + 3),
                     base_stiffness - 0.5 * nominal_damping * rotor_omega * sin_2wt, 0, 0, 0],
                    [0, -nominal_damping * rotor_omega * sin_wt,
                     nominal_damping * rotor_omega * cos_wt,
                     4 * centrifugal_stiffness + 4 * nominal_stiffness, 0, 0],
                    [0, 0, 0, 0, hub_stiffness_x, 0],
                    [0, 0, 0, 0, 0, hub_stiffness_y]
                ])
        else:
            raise ValueError(f"Unsupported damper_activation for 4-blade H2B: {damper_activation}")
    else:
        raise ValueError(f"Unsupported damper_connection for 4-blade: {damper_connection}")

    return damping_matrix_C, stiffness_matrix_K


def _build_matrices_5_blade(
    blade_inertia: float, blade_static_moment: float, lag_hinge_offset: float,
    nominal_damping: float, nominal_stiffness: float,
    hub_damping_x: float, hub_damping_y: float,
    hub_stiffness_x: float, hub_stiffness_y: float,
    damper_activation: str
) -> Tuple[Callable[[float, float], np.ndarray], Callable[[float, float], np.ndarray]]:
    """Build damping and stiffness matrices for 5-blade rotor."""

    if damper_activation == "ALL":
        def damping_matrix_C(time: float, rotor_omega: float) -> np.ndarray:
            return np.array([
                [5 * nominal_damping, 0, 0, 0, 0, 0, 0],
                [0, 2.5 * nominal_damping, 5 * blade_inertia * rotor_omega, 0, 0, 0, 0],
                [0, -5 * blade_inertia * rotor_omega, 2.5 * nominal_damping, 0, 0, 0, 0],
                [0, 0, 0, 2.5 * nominal_damping, 10 * blade_inertia * rotor_omega, 0, 0],
                [0, 0, 0, -10 * blade_inertia * rotor_omega, 2.5 * nominal_damping, 0, 0],
                [0, 0, 0, 0, 0, hub_damping_x, 0],
                [0, 0, 0, 0, 0, 0, hub_damping_y]
            ])

        def stiffness_matrix_K(time: float, rotor_omega: float) -> np.ndarray:
            centrifugal_stiffness = blade_static_moment * lag_hinge_offset * rotor_omega**2
            gyroscopic_term_1 = 2.5 * nominal_damping * rotor_omega
            gyroscopic_term_2 = 5 * nominal_damping * rotor_omega

            dynamic_stiffness_1 = 2.5 * nominal_stiffness - 2.5 * blade_inertia * rotor_omega**2 + \
                                  2.5 * centrifugal_stiffness
            dynamic_stiffness_2 = 2.5 * nominal_stiffness - 10 * blade_inertia * rotor_omega**2 + \
                                  2.5 * centrifugal_stiffness

            return np.array([
                [5 * centrifugal_stiffness + 5 * nominal_stiffness, 0, 0, 0, 0, 0, 0],
                [0, dynamic_stiffness_1, gyroscopic_term_1, 0, 0, 0, 0],
                [0, -gyroscopic_term_1, dynamic_stiffness_1, 0, 0, 0, 0],
                [0, 0, 0, dynamic_stiffness_2, gyroscopic_term_2, 0, 0],
                [0, 0, 0, -gyroscopic_term_2, dynamic_stiffness_2, 0, 0],
                [0, 0, 0, 0, 0, hub_stiffness_x, 0],
                [0, 0, 0, 0, 0, 0, hub_stiffness_y]
            ])

    elif damper_activation == "ODI":
        def damping_matrix_C(time: float, rotor_omega: float) -> np.ndarray:
            wt = rotor_omega * time
            cos_wt = np.cos(wt)
            sin_wt = np.sin(wt)
            cos_2wt = np.cos(2 * wt)
            sin_2wt = np.sin(2 * wt)
            cos_4wt = np.cos(4 * wt)
            sin_4wt = np.sin(4 * wt)

            return np.array([
                [4 * nominal_damping, -nominal_damping * cos_wt, -nominal_damping * sin_wt,
                 -nominal_damping * cos_2wt, -nominal_damping * sin_2wt, 0, 0],
                [-nominal_damping * cos_wt, -0.5 * nominal_damping * (cos_2wt - 4),
                 5 * blade_inertia * rotor_omega - 0.5 * nominal_damping * sin_2wt,
                 -nominal_damping * cos_wt * (2 * cos_wt**2 - 1),
                 2 * nominal_damping * sin_wt * (sin_wt**2 - 1), 0, 0],
                [-nominal_damping * sin_wt, -0.5 * nominal_damping * sin_2wt - 5 * blade_inertia * rotor_omega,
                 0.5 * nominal_damping * (cos_2wt + 4),
                 nominal_damping * sin_wt * (2 * sin_wt**2 - 1),
                 2 * nominal_damping * cos_wt * (cos_wt**2 - 1), 0, 0],
                [-nominal_damping * cos_2wt, -nominal_damping * cos_wt * (2 * cos_wt**2 - 1),
                 nominal_damping * sin_wt * (2 * sin_wt**2 - 1),
                 -0.5 * nominal_damping * (cos_4wt - 4),
                 10 * blade_inertia * rotor_omega - 0.5 * nominal_damping * sin_4wt, 0, 0],
                [-nominal_damping * sin_2wt, 2 * nominal_damping * sin_wt * (sin_wt**2 - 1),
                 2 * nominal_damping * cos_wt * (cos_wt**2 - 1),
                 -0.5 * nominal_damping * sin_4wt - 10 * blade_inertia * rotor_omega,
                 0.5 * nominal_damping * (cos_4wt + 4), 0, 0],
                [0, 0, 0, 0, 0, hub_damping_x, 0],
                [0, 0, 0, 0, 0, 0, hub_damping_y]
            ])

        def stiffness_matrix_K(time: float, rotor_omega: float) -> np.ndarray:
            wt = rotor_omega * time
            cos_wt = np.cos(wt)
            sin_wt = np.sin(wt)
            cos_2wt = np.cos(2 * wt)
            sin_2wt = np.sin(2 * wt)
            cos_4wt = np.cos(4 * wt)
            sin_4wt = np.sin(4 * wt)

            centrifugal_stiffness = blade_static_moment * lag_hinge_offset * rotor_omega**2
            base_stiffness_1 = 2.5 * nominal_stiffness - 2.5 * blade_inertia * rotor_omega**2 + \
                               2.5 * centrifugal_stiffness
            base_stiffness_2 = 2.5 * nominal_stiffness - 10 * blade_inertia * rotor_omega**2 + \
                               2.5 * centrifugal_stiffness

            return np.array([
                [5 * centrifugal_stiffness + 5 * nominal_stiffness,
                 nominal_damping * rotor_omega * sin_wt,
                 -nominal_damping * rotor_omega * cos_wt,
                 2 * nominal_damping * rotor_omega * sin_2wt,
                 -2 * nominal_damping * rotor_omega * cos_2wt, 0, 0],
                [0, base_stiffness_1 + 0.5 * nominal_damping * rotor_omega * sin_2wt,
                 -0.5 * nominal_damping * rotor_omega * (cos_2wt - 4),
                 -4 * nominal_damping * rotor_omega * sin_wt * (sin_wt**2 - 1),
                 -2 * nominal_damping * rotor_omega * cos_wt * (2 * cos_wt**2 - 1), 0, 0],
                [0, -0.5 * nominal_damping * rotor_omega * (cos_2wt + 4),
                 base_stiffness_1 - 0.5 * nominal_damping * rotor_omega * sin_2wt,
                 -4 * nominal_damping * rotor_omega * cos_wt * (cos_wt**2 - 1),
                 2 * nominal_damping * rotor_omega * sin_wt * (2 * sin_wt**2 - 1), 0, 0],
                [0, -nominal_damping * rotor_omega * sin_wt * (2 * sin_wt**2 - 1),
                 -nominal_damping * rotor_omega * cos_wt * (2 * cos_wt**2 - 1),
                 base_stiffness_2 + nominal_damping * rotor_omega * sin_4wt,
                 -nominal_damping * rotor_omega * (cos_4wt - 4), 0, 0],
                [0, -2 * nominal_damping * rotor_omega * cos_wt * (cos_wt**2 - 1),
                 2 * nominal_damping * rotor_omega * sin_wt * (sin_wt**2 - 1),
                 -nominal_damping * rotor_omega * (cos_4wt + 4),
                 base_stiffness_2 - nominal_damping * rotor_omega * sin_4wt, 0, 0],
                [0, 0, 0, 0, 0, hub_stiffness_x, 0],
                [0, 0, 0, 0, 0, 0, hub_stiffness_y]
            ])
    else:
        raise ValueError(f"Unsupported damper_activation for 5-blade: {damper_activation}")

    return damping_matrix_C, stiffness_matrix_K


if __name__ == "__main__":
    # Test matrix construction
    from ..config.problem_definition import create_default_problem
    
    problem = create_default_problem()
    
    # Test mass matrix
    M = build_mass_matrix(problem)
    print(f"Mass matrix shape: {M.shape}")
    
    # Test damping and stiffness matrices
    C_func, K_func = build_damping_matrix(problem)
    C = C_func(0.0, 25.0)
    K = K_func(0.0, 25.0)
    print(f"Damping matrix shape: {C.shape}")
    print(f"Stiffness matrix shape: {K.shape}")
