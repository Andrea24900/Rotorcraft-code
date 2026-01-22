"""Problem definition module for rotor dynamics analysis.

This module defines the data structures for problem configuration and
rotor characteristics, converted from MATLAB problem_definition.m and
rotor_definition_script.m.
"""

from dataclasses import dataclass, field
from typing import Literal
import numpy as np


@dataclass
class RotorCharacteristics:
    """Physical characteristics of the rotor system."""
    
    # Blade properties
    blade_mass_Mb: float = 94.9  # kg
    blade_static_Sb: float = 289.1  # kg*m
    blade_inertia_Ib: float = 1085.0  # kg*m^2
    lag_hinge_offset_e: float = 0.3048  # m
    
    # Hub properties
    hub_mass_Mx: float = 8027.0  # kg
    hub_mass_My: float = 3284.0  # kg
    hub_damping_Cx: float = 51420.0  # Ns/m
    hub_damping_Cy: float = 25710.0  # Ns/m
    hub_stiffness_Kx: float = 1241000.0  # N/m
    hub_stiffness_Ky: float = 1241000.0  # N/m
    
    # Damper properties (will be set based on damper_connection)
    nominal_damping_Cd: float = 4067.0  # Nms/rad or Ns/m
    nominal_stiffness_Kd: float = 0.0  # Nm/rad or N/m
    
    # Geometric properties (for B2B configuration)
    delta_psi: float = 0.0  # rad
    phi: float = 0.0  # rad
    hub_to_A_a: float = 0.46  # m
    hub_to_B_b: float = 0.23  # m
    hinge_to_hinge_c: float = 0.0  # m
    length_undef_l0: float = 0.0  # m
    
    # Equilibrium variables (for B2B configuration)
    zeta_E: float = 0.0  # rad
    length_eq_lE: float = 0.0  # m
    gamma_E: float = 0.0  # rad
    
    # Internal flag
    isotropic_hub: int = 0  # deprecated, do not change


@dataclass
class ProblemDefinition:
    """Complete problem definition for rotor dynamics analysis."""
    
    # Rotor configuration
    number_blades: int = 4
    damper_connection: Literal["H2B", "B2B"] = "H2B"
    
    damper_activation: Literal["ALL", "ODI", "CUSTOM"] = "ALL"  # ODI = One Damper Inoperative
    
    # implement a custom damper activation logic (to be bettered)
    
    damper_activation_vector: np.ndarray = field(default_factory=lambda: np.array([1,1,1,1]))

    # Solution parameters
    lower_rotor_RPM: float = 20.0
    higher_rotor_RPM: float = 400.0
    number_points: int = 200
    
    # Solver configuration
    required_solver: Literal["LTI", "LTP", "HD"] = "LTI"
    number_harmonics: int = 1  # for HD and modal participation
    hd_use_complex: bool = True  # Use complex formulation for HD solver
    
    # Continuation method parameters
    continuation: Literal["YES", "NO"] = "NO"
    continuation_tolerance: float = 1e-3
    continuation_max_iter: int = 100
    step_h: float = 1e-4
    time_samples: int = 10  # must be chosen in accordance with Nyquist criterion
    # based on the number of harmonics that are required for the fft of A and the
    # harmonic decomposition
    solution_direction: Literal["L2R", "R2L"] = "R2L"

    # ODE integration tolerances (for LTP monodromy matrix computation)
    ode_rtol: float = 1e-5  # relative tolerance
    ode_atol: float = 1e-7  # absolute tolerance
    
    # Rotor physical characteristics
    rotor_characteristics: RotorCharacteristics = field(default_factory=RotorCharacteristics)
    
    def __post_init__(self):
        """Initialize derived parameters after instantiation."""
        self._initialize_rotor_characteristics()
        if self.damper_activation == "CUSTOM":
            if len(self.damper_activation_vector) != self.number_blades:
                raise ValueError(f"damper_activation_vector length must match number_blades")

    
    def _initialize_rotor_characteristics(self):
        """Initialize rotor characteristics based on configuration."""
        import numpy as np
        
        # Set damper properties based on connection type
        if self.damper_connection == "H2B":
            # Hub-to-Blade configuration
            self.rotor_characteristics.nominal_damping_Cd = 4067  # Nms/rad #nominal is 4067
            self.rotor_characteristics.nominal_stiffness_Kd = 0.0  # Nm/rad
        
        elif self.damper_connection == "B2B":
            # Blade-to-Blade configuration
            self.rotor_characteristics.nominal_damping_Cd = 38000.0  # Ns/m
            self.rotor_characteristics.nominal_stiffness_Kd = 0.0  # N/m
            
            # Geometric relationships for B2B
            self.rotor_characteristics.delta_psi = 2 * np.pi / self.number_blades
            self.rotor_characteristics.phi = (
                np.pi - self.rotor_characteristics.delta_psi
            ) / 2
            
            self.rotor_characteristics.hinge_to_hinge_c = (
                2 * self.rotor_characteristics.lag_hinge_offset_e *
                np.sin(self.rotor_characteristics.delta_psi / 2)
            )
            
            self.rotor_characteristics.length_undef_l0 = np.sqrt(
                self.rotor_characteristics.hub_to_B_b**2 +
                self.rotor_characteristics.hub_to_A_a**2
            )
            
            # Equilibrium variables (zeta_E is set to 0.0 by default)
            zeta_E = self.rotor_characteristics.zeta_E
            phi = self.rotor_characteristics.phi
            a = self.rotor_characteristics.hub_to_A_a
            b = self.rotor_characteristics.hub_to_B_b
            c = self.rotor_characteristics.hinge_to_hinge_c
            
            self.rotor_characteristics.length_eq_lE = np.sqrt(
                (c + a * np.cos(zeta_E - phi) + b * np.cos(zeta_E + phi))**2 +
                (a * np.sin(zeta_E - phi) + b * np.sin(zeta_E + phi))**2
            )
            
            self.rotor_characteristics.gamma_E = np.arctan2(
                -(a * np.sin(zeta_E - phi) + b * np.sin(zeta_E + phi)),
                c + a * np.cos(zeta_E - phi) + b * np.cos(zeta_E + phi)
            )


def create_default_problem() -> ProblemDefinition:
    """Create a default problem definition.
    
    Returns:
        ProblemDefinition: Default problem configuration
    """
    return ProblemDefinition()


# Example usage
if __name__ == "__main__":
    # Create default problem
    problem = create_default_problem()
    print(f"Number of blades: {problem.number_blades}")
    print(f"Solver: {problem.required_solver}")
    print(f"RPM range: {problem.lower_rotor_RPM} - {problem.higher_rotor_RPM}")
    print(f"Damper connection: {problem.damper_connection}")
    print(f"Nominal damping: {problem.rotor_characteristics.nominal_damping_Cd}")
