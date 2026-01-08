"""Unit tests for problem_definition module.

Run with: pytest tests/test_problem_definition.py
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from config.problem_definition import ProblemDefinition, RotorCharacteristics, create_default_problem


class TestRotorCharacteristics:
    """Test RotorCharacteristics dataclass."""
    
    def test_default_values(self):
        """Test that default values are set correctly."""
        rotor_char = RotorCharacteristics()
        
        assert rotor_char.blade_mass_Mb == 94.9
        assert rotor_char.blade_static_Sb == 289.1
        assert rotor_char.blade_inertia_Ib == 1085.0
        assert rotor_char.hub_mass_Mx == 8027.0
    
    def test_custom_values(self):
        """Test custom initialization."""
        rotor_char = RotorCharacteristics(
            blade_mass_Mb=100.0,
            hub_mass_Mx=8000.0
        )
        
        assert rotor_char.blade_mass_Mb == 100.0
        assert rotor_char.hub_mass_Mx == 8000.0


class TestProblemDefinition:
    """Test ProblemDefinition dataclass."""
    
    def test_default_creation(self):
        """Test default problem creation."""
        problem = create_default_problem()
        
        assert problem.number_blades == 4
        assert problem.damper_connection == "H2B"
        assert problem.required_solver == "HD"
        assert problem.lower_rotor_RPM == 20.0
        assert problem.higher_rotor_RPM == 400.0
    
    def test_h2b_configuration(self):
        """Test Hub-to-Blade damper configuration."""
        problem = ProblemDefinition(damper_connection="H2B")
        
        assert problem.rotor_characteristics.nominal_damping_Cd == 4067.0
        assert problem.rotor_characteristics.nominal_stiffness_Kd == 0.0
    
    def test_b2b_configuration(self):
        """Test Blade-to-Blade damper configuration."""
        problem = ProblemDefinition(
            damper_connection="B2B",
            number_blades=4
        )
        
        # Check damper values
        assert problem.rotor_characteristics.nominal_damping_Cd == 36500.0
        
        # Check geometric calculations
        assert problem.rotor_characteristics.delta_psi == pytest.approx(np.pi / 2)
        assert problem.rotor_characteristics.phi > 0
        assert problem.rotor_characteristics.hinge_to_hinge_c > 0
        assert problem.rotor_characteristics.length_undef_l0 > 0
    
    def test_number_blades_validation(self):
        """Test that valid blade numbers are accepted."""
        for n_blades in [3, 4, 5, 7]:
            problem = ProblemDefinition(number_blades=n_blades)
            assert problem.number_blades == n_blades
    
    def test_solver_types(self):
        """Test different solver configurations."""
        for solver in ["LTI", "LTP", "HD"]:
            problem = ProblemDefinition(required_solver=solver)
            assert problem.required_solver == solver
    
    def test_rpm_range(self):
        """Test RPM range configuration."""
        problem = ProblemDefinition(
            lower_rotor_RPM=50.0,
            higher_rotor_RPM=500.0,
            number_points=100
        )
        
        assert problem.lower_rotor_RPM == 50.0
        assert problem.higher_rotor_RPM == 500.0
        assert problem.number_points == 100
    
    def test_harmonics_setting(self):
        """Test number of harmonics for HD analysis."""
        problem = ProblemDefinition(
            required_solver="HD",
            number_harmonics=3
        )
        
        assert problem.number_harmonics == 3


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
