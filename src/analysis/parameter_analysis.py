"""Parametric stability analysis module.

This module provides tools for sweeping a parameter (e.g., damper damping)
across a range of values while computing stability at each (OMEGA, parameter) point.

Key features:
- Automatic detection of periodic vs. time-invariant systems
- Solver validation to ensure correct method is used for the problem type
- Extraction of least damped eigenvalue for stability mapping
- 2D visualization of stability boundaries
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List, Literal, Optional, Callable, Union
from matplotlib.figure import Figure

from ..config.problem_definition import ProblemDefinition
from ..rotor_build import RotorBuild
from ..stability_analysis.lti_stability import LTIStability
from ..stability_analysis.ltp_stability import LTPStability
from ..stability_analysis.hd_stability import HDStability


@dataclass
class ParametricPoint:
    """Result at a single (OMEGA, parameter) point."""

    OMEGA: float = 0.0  # rad/s
    OMEGA_RPM: float = 0.0  # RPM
    parameter_value: float = 0.0
    least_damped_eigenvalue: complex = 0.0 + 0.0j  # Full eigenvalue
    damping: float = 0.0  # Real part (stability indicator)
    frequency: float = 0.0  # Imaginary part (Hz or rad/s)


class ParametricAnalysis:
    """Parametric sweep analysis over OMEGA and a system parameter to find the least damped mode.

    This class performs stability analysis across a 2D grid of (OMEGA, parameter)
    values. It automatically validates that the chosen solver is appropriate
    for the problem type (LTI for symmetric systems, LTP/HD for periodic systems).

    Attributes:
        base_problem: Reference problem definition
        parameter_name: Name of the parameter being swept
        parameter_values: Array of parameter values to sweep
        omega_values: Array of OMEGA values (rad/s)
        results: 2D list of ParametricPoint results [param_idx][omega_idx]
        solver_type: The solver used ("LTI", "LTP", or "HD")

    Example:
        >>> problem = create_default_problem()
        >>> analysis = ParametricAnalysis(problem)
        >>> analysis.sweep(
        ...     solver_type="LTI",
        ...     omega_values_RPM=np.linspace(50, 400, 50),
        ...     parameter_name="nominal_damping_Cd",
        ...     parameter_values=np.linspace(1000, 8000, 20)
        ... )
        >>> fig = analysis.plot_stability_map()
    """

    def __init__(self, base_problem: ProblemDefinition):
        """Initialize parametric analysis.

        Args:
            base_problem: Base problem definition to use as template.
                         Will be deep-copied for each parameter value.
        """
        self.base_problem = base_problem
        self.parameter_name: str = ""
        self.parameter_values: np.ndarray = np.array([])
        self.omega_values: np.ndarray = np.array([])
        self.omega_values_RPM: np.ndarray = np.array([])
        self.results: List[List[ParametricPoint]] = []
        self.solver_type: str = ""

    @staticmethod
    def check_system_periodicity(problem: ProblemDefinition) -> bool:
        """Check if the system configuration results in a periodic system.

        A system is periodic (time-varying) when damper properties are
        asymmetric across blades. This occurs when:
        - damper_activation is "ODI" (One Damper Inoperative)
        - damper_activation is "CUSTOM" with non-uniform activation vector

        A system is time-invariant (LTI) when:
        - damper_activation is "ALL" (all dampers have same properties)
        - damper_activation is "CUSTOM" with uniform activation vector
          (e.g., [0.5, 0.5, 0.5, 0.5] - all equal, just scaled)

        The damper_activation_vector is a multiplier applied to nominal_damping_Cd:
            effective_damping[i] = nominal_damping_Cd * damper_activation_vector[i]

        Examples:
            - [1, 1, 1, 1]     -> all dampers at 100% -> symmetric (LTI)
            - [0.5, 0.5, 0.5, 0.5] -> all at 50% -> symmetric (LTI)
            - [0.5, 1, 1, 1]  -> first damper at 50%, others 100% -> periodic (LTP)
            - [0, 1, 1, 1]    -> first damper off -> periodic (LTP)

        Args:
            problem: Problem definition to check

        Returns:
            True if system is periodic (requires LTP/HD), False if LTI is valid
        """
        if problem.rotor_characteristics.damper_activation == "ALL":
            return False  # Symmetric system -> LTI valid

        elif problem.rotor_characteristics.damper_activation == "ODI":
            return True  # One damper out -> periodic system

        elif problem.rotor_characteristics.damper_activation == "CUSTOM":
            # Check if activation vector is uniform (all elements equal)
            # The actual value doesn't matter - only uniformity determines periodicity
            activation = problem.rotor_characteristics.damper_activation_vector
            if len(activation) == 0:
                return True  # Default to periodic (safer)

            # System is periodic if not all dampers have the same multiplier
            is_uniform = np.allclose(activation, activation[0])
            return not is_uniform

        return True  # Default to periodic (LTP) - safer assumption

    @staticmethod
    def recommend_solver(problem: ProblemDefinition) -> Literal["LTI", "LTP"]:
        """Recommend the appropriate solver for the given problem configuration.

        Returns the most efficient valid solver:
        - LTI if system is symmetric (faster, no ODE integration)
        - LTP if system is periodic (more general, handles asymmetric dampers)

        Note: HD is also valid for periodic systems but LTP is recommended
        as the default periodic solver.

        Args:
            problem: Problem definition to check

        Returns:
            "LTI" for symmetric systems, "LTP" for periodic systems
        """
        is_periodic = ParametricAnalysis.check_system_periodicity(problem)
        return "LTP" if is_periodic else "LTI"

    @staticmethod
    def validate_solver_for_problem(
        solver_type: Literal["LTI", "LTP", "HD"],
        problem: ProblemDefinition,
        raise_error: bool = True
    ) -> bool:
        """Validate that the chosen solver is appropriate for the problem.

        LTI solver requires a time-invariant (non-periodic) system.
        LTP and HD solvers can handle both periodic and non-periodic systems.

        Args:
            solver_type: The solver to validate ("LTI", "LTP", or "HD")
            problem: Problem definition to check
            raise_error: If True, raise ValueError for invalid combinations.
                        If False, just return False.

        Returns:
            True if solver is valid for the problem

        Raises:
            ValueError: If solver_type is "LTI" but system is periodic
                       (only if raise_error=True)
        """
        is_periodic = ParametricAnalysis.check_system_periodicity(problem)

        if solver_type == "LTI" and is_periodic:
            msg = (
                f"Cannot use LTI solver for periodic system.\n"
                f"  - damper_activation: {problem.rotor_characteristics.damper_activation}\n"
            )
            if problem.rotor_characteristics.damper_activation == "CUSTOM":
                msg += f"  - damper_activation_vector: {problem.rotor_characteristics.damper_activation_vector}\n"
            msg += (
                f"The system has asymmetric damper configuration, making it time-periodic.\n"
                f"Use LTP or HD solver instead."
            )
            if raise_error:
                raise ValueError(msg)
            return False

        return True

    def sweep(
        self,
        solver_type: Literal["LTI", "LTP", "HD"],
        omega_values_RPM: np.ndarray,
        parameter_name: str,
        parameter_values: np.ndarray,
        parameter_setter: Optional[Callable[[ProblemDefinition, float], None]] = None,
        verbose: bool = True
    ) -> "ParametricAnalysis":
        """Run parametric sweep over OMEGA and parameter values.

        For each parameter value:
        1. Create a copy of the problem and modify the parameter
        2. Validate that the solver is appropriate for the resulting system
        3. Build the rotor and run single-point analysis at each OMEGA
        4. Extract and store the least damped eigenvalue

        Args:
            solver_type: Solver to use ("LTI", "LTP", or "HD")
            omega_values_RPM: Array of rotor speeds in RPM
            parameter_name: Name of parameter to sweep. Can be:
                           - Attribute of rotor_characteristics (e.g., "nominal_damping_Cd")
                           - "damper_activation_vector[i]" for individual damper control
            parameter_values: Array of parameter values to sweep
            parameter_setter: Optional custom function to set parameter.
                            Signature: (problem, value) -> None.
                            If None, uses setattr on rotor_characteristics.
            verbose: Print progress information

        Returns:
            Self with populated results

        Raises:
            ValueError: If solver_type is invalid for any parameter configuration
        """
        self.solver_type = solver_type
        self.parameter_name = parameter_name
        self.parameter_values = parameter_values
        self.omega_values_RPM = omega_values_RPM
        self.omega_values = omega_values_RPM * (2 * np.pi / 60)  # Convert to rad/s
        self.results = []

        n_params = len(parameter_values)
        n_omegas = len(omega_values_RPM)

        if verbose:
            print(f"Parametric Analysis: {parameter_name}")
            print(f"  Solver: {solver_type}")
            print(f"  OMEGA range: {omega_values_RPM[0]:.1f} - {omega_values_RPM[-1]:.1f} RPM ({n_omegas} points)")
            print(f"  Parameter range: {parameter_values[0]:.4g} - {parameter_values[-1]:.4g} ({n_params} points)")
            print(f"  Total evaluations: {n_params * n_omegas}")
            print()

        for i, param_val in enumerate(parameter_values):
            if verbose:
                print(f"  [{i+1}/{n_params}] {parameter_name} = {param_val:.4g}", end="")

            # Create deep copy of problem and modify parameter
            problem_copy = copy.deepcopy(self.base_problem)

            if parameter_setter is not None:
                parameter_setter(problem_copy, param_val)
            else:
                self._set_parameter(problem_copy, parameter_name, param_val)

            # Validate solver is appropriate for this configuration
            is_periodic = self.check_system_periodicity(problem_copy)
            self.validate_solver_for_problem(solver_type, problem_copy, raise_error=True)

            if verbose and is_periodic:
                print(" [PERIODIC]", end="")

            # Build rotor with modified problem
            rotor = RotorBuild(problem_copy, auto_build=True)

            # Create solver
            solver = self._create_solver(solver_type, rotor)

            # Compute at each OMEGA point
            param_results = []
            for j, omega in enumerate(self.omega_values):
                result = self._compute_single_point(solver, omega, param_val, solver_type)
                param_results.append(result)

            self.results.append(param_results)

            if verbose:
                # Report max damping for this parameter value
                max_damping = max(r.damping for r in param_results)
                status = "UNSTABLE" if max_damping > 0 else "stable"
                print(f" -> max damping: {max_damping:.4f} ({status})")

        if verbose:
            print("\nSweep complete.")

        return self

    def _set_parameter(
        self,
        problem: ProblemDefinition,
        parameter_name: str,
        value: float
    ) -> None:
        """Set a parameter value on the problem.

        Handles both simple attributes and indexed access.

        Args:
            problem: Problem to modify
            parameter_name: Parameter name (e.g., "nominal_damping_Cd" or
                          "damper_activation_vector[0]")
            value: Value to set
        """
        # Check for indexed access pattern
        if "[" in parameter_name and "]" in parameter_name:
            # Parse "array_name[index]"
            base_name = parameter_name.split("[")[0]
            index = int(parameter_name.split("[")[1].split("]")[0])

            if hasattr(problem.rotor_characteristics, base_name):
                arr = getattr(problem.rotor_characteristics, base_name)
                arr[index] = value
            elif hasattr(problem, base_name):
                arr = getattr(problem, base_name)
                arr[index] = value
            else:
                raise AttributeError(
                    f"Neither problem nor rotor_characteristics has attribute '{base_name}'"
                )
        else:
            # Simple attribute access
            if hasattr(problem.rotor_characteristics, parameter_name):
                setattr(problem.rotor_characteristics, parameter_name, value)
            elif hasattr(problem, parameter_name):
                setattr(problem, parameter_name, value)
            else:
                raise AttributeError(
                    f"Neither problem nor rotor_characteristics has attribute '{parameter_name}'"
                )

    def _create_solver(
        self,
        solver_type: Literal["LTI", "LTP", "HD"],
        rotor: RotorBuild
    ):
        """Create appropriate solver instance.

        Args:
            solver_type: Type of solver to create
            rotor: Rotor build object

        Returns:
            Solver instance (LTIStability, LTPStability, or HDStability)
        """
        if solver_type == "LTI":
            return LTIStability(rotor)
        elif solver_type == "LTP":
            return LTPStability(rotor)
        elif solver_type == "HD":
            return HDStability(rotor)
        else:
            raise ValueError(f"Unknown solver type: {solver_type}")

    def _compute_single_point(
        self,
        solver,
        omega: float,
        param_val: float,
        solver_type: str
    ) -> ParametricPoint:
        """Compute stability at a single (OMEGA, parameter) point.

        Args:
            solver: Solver instance
            omega: Rotor speed in rad/s
            param_val: Current parameter value
            solver_type: Type of solver being used

        Returns:
            ParametricPoint with least damped eigenvalue extracted
        """
        # Call appropriate single-point method
        if solver_type == "LTI":
            solution = solver.eigen_single_point(omega)
            eigenvalues = solution['eigenvalues']
        elif solver_type == "LTP":
            solution = solver.ltp_single_point(omega)
            # For LTP, construct complex eigenvalues from characteristic exponents
            eigenvalues = solution['real_char_exp'] + 1j * solution['imag_char_exp']
        elif solver_type == "HD":
            solution = solver.hd_single_point(omega)
            eigenvalues = solution['eigenvalues']

        # Find least damped eigenvalue (most positive real part)
        damping_values = np.real(eigenvalues)
        least_damped_idx = np.argmax(damping_values)
        least_damped = eigenvalues[least_damped_idx]

        return ParametricPoint(
            OMEGA=omega,
            OMEGA_RPM=omega * 60 / (2 * np.pi),
            parameter_value=param_val,
            least_damped_eigenvalue=least_damped,
            damping=np.real(least_damped),
            frequency=np.imag(least_damped)
        )

    def get_damping_matrix(self) -> np.ndarray:
        """Get 2D array of damping values.

        Returns:
            Array of shape (n_params, n_omegas) with damping values
        """
        n_params = len(self.parameter_values)
        n_omegas = len(self.omega_values)

        damping = np.zeros((n_params, n_omegas))
        for i in range(n_params):
            for j in range(n_omegas):
                damping[i, j] = self.results[i][j].damping

        return damping

    def get_frequency_matrix(self) -> np.ndarray:
        """Get 2D array of frequency values.

        Returns:
            Array of shape (n_params, n_omegas) with frequency values
        """
        n_params = len(self.parameter_values)
        n_omegas = len(self.omega_values)

        frequency = np.zeros((n_params, n_omegas))
        for i in range(n_params):
            for j in range(n_omegas):
                frequency[i, j] = self.results[i][j].frequency

        return frequency

    def plot_stability_map(
        self,
        figsize: tuple = (10, 8),
        cmap: str = "RdYlGn_r",
        levels: int = 50,
        show_stability_boundary: bool = True,
        title: Optional[str] = None
    ) -> Figure:
        """Create contour plot of damping vs (OMEGA, parameter).

        Args:
            figsize: Figure size (width, height)
            cmap: Colormap name (default: RdYlGn_r, red=unstable, green=stable)
            levels: Number of contour levels
            show_stability_boundary: If True, draw black line at damping=0
            title: Optional custom title

        Returns:
            Matplotlib Figure object
        """
        damping = self.get_damping_matrix()

        # Create meshgrid for plotting
        OMEGA_grid, PARAM_grid = np.meshgrid(self.omega_values_RPM, self.parameter_values)

        fig, ax = plt.subplots(figsize=figsize)

        # Contour plot
        cf = ax.contourf(OMEGA_grid, PARAM_grid, damping, levels=levels, cmap=cmap)
        cbar = plt.colorbar(cf, ax=ax)
        cbar.set_label("Damping (Real part of eigenvalue) [1/s]")

        # Stability boundary
        if show_stability_boundary:
            cs = ax.contour(OMEGA_grid, PARAM_grid, damping, levels=[0],
                           colors='black', linewidths=2, linestyles='--')
            ax.clabel(cs, fmt='σ=0', fontsize=10)

        ax.set_xlabel("Rotor Speed [RPM]")
        ax.set_ylabel(self._format_parameter_label())

        if title is None:
            title = f"Stability Map ({self.solver_type} solver)"
        ax.set_title(title)

        ax.grid(True, alpha=0.3)

        return fig

    def plot_frequency_map(
        self,
        figsize: tuple = (10, 8),
        cmap: str = "viridis",
        levels: int = 50,
        title: Optional[str] = None
    ) -> Figure:
        """Create contour plot of frequency vs (OMEGA, parameter).

        Args:
            figsize: Figure size (width, height)
            cmap: Colormap name
            levels: Number of contour levels
            title: Optional custom title

        Returns:
            Matplotlib Figure object
        """
        frequency = self.get_frequency_matrix()

        # Create meshgrid for plotting
        OMEGA_grid, PARAM_grid = np.meshgrid(self.omega_values_RPM, self.parameter_values)

        fig, ax = plt.subplots(figsize=figsize)

        # Contour plot
        cf = ax.contourf(OMEGA_grid, PARAM_grid, frequency, levels=levels, cmap=cmap)
        cbar = plt.colorbar(cf, ax=ax)
        cbar.set_label("Frequency (Imag part of eigenvalue) [rad/s]")

        ax.set_xlabel("Rotor Speed [RPM]")
        ax.set_ylabel(self._format_parameter_label())

        if title is None:
            title = f"Frequency Map ({self.solver_type} solver)"
        ax.set_title(title)

        ax.grid(True, alpha=0.3)

        return fig

    def plot_damping_slices(
        self,
        param_indices: Optional[List[int]] = None,
        figsize: tuple = (10, 6),
        title: Optional[str] = None
    ) -> Figure:
        """Plot damping vs OMEGA for selected parameter values.

        Args:
            param_indices: List of parameter indices to plot.
                          If None, plots 5 evenly spaced values.
            figsize: Figure size
            title: Optional custom title

        Returns:
            Matplotlib Figure object
        """
        if param_indices is None:
            n = len(self.parameter_values)
            param_indices = list(np.linspace(0, n-1, min(5, n), dtype=int))

        fig, ax = plt.subplots(figsize=figsize)

        for idx in param_indices:
            damping = [self.results[idx][j].damping for j in range(len(self.omega_values))]
            label = f"{self.parameter_name} = {self.parameter_values[idx]:.4g}"
            ax.plot(self.omega_values_RPM, damping, label=label, linewidth=1.5)

        ax.axhline(y=0, color='k', linestyle='--', linewidth=1, label='Stability boundary')
        ax.set_xlabel("Rotor Speed [RPM]")
        ax.set_ylabel("Damping [1/s]")

        if title is None:
            title = f"Damping vs Rotor Speed ({self.solver_type} solver)"
        ax.set_title(title)

        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)

        return fig

    def _format_parameter_label(self) -> str:
        """Format parameter name for axis label."""
        # Common parameter name mappings
        label_map = {
            "nominal_damping_Cd": "Damper Coefficient [Nms/rad]",
            "nominal_stiffness_Kd": "Damper Stiffness [Nm/rad]",
            "hub_damping_Cx": "Hub Damping X [Ns/m]",
            "hub_damping_Cy": "Hub Damping Y [Ns/m]",
            "hub_stiffness_Kx": "Hub Stiffness X [N/m]",
            "hub_stiffness_Ky": "Hub Stiffness Y [N/m]",
        }

        return label_map.get(self.parameter_name, self.parameter_name)

    def find_stability_boundary(
        self,
        omega_RPM: float,
        tol: float = 1e-3
    ) -> Optional[float]:
        """Find parameter value at stability boundary for given OMEGA.

        Uses linear interpolation between computed points.

        Args:
            omega_RPM: Rotor speed in RPM
            tol: Tolerance for considering damping as zero

        Returns:
            Parameter value at stability boundary, or None if not found
        """
        # Find closest OMEGA index
        omega_idx = np.argmin(np.abs(self.omega_values_RPM - omega_RPM))

        # Get damping values at this OMEGA for all parameters
        dampings = [self.results[i][omega_idx].damping for i in range(len(self.parameter_values))]

        # Find sign changes (stability boundary crossings)
        for i in range(len(dampings) - 1):
            if dampings[i] * dampings[i+1] < 0:  # Sign change
                # Linear interpolation
                p1, p2 = self.parameter_values[i], self.parameter_values[i+1]
                d1, d2 = dampings[i], dampings[i+1]

                # p at d=0: p = p1 + (p2-p1) * (-d1)/(d2-d1)
                p_boundary = p1 + (p2 - p1) * (-d1) / (d2 - d1)
                return p_boundary

        return None


if __name__ == "__main__":
    # Example usage
    from ..config.problem_definition import create_default_problem

    print("=" * 60)
    print("Parametric Analysis Example")
    print("=" * 60)

    # Create base problem
    problem = create_default_problem()
    problem.rotor_characteristics.damper_activation = "ALL"  # Symmetric system -> LTI valid

    # Create analysis object
    analysis = ParametricAnalysis(problem)

    # Define sweep parameters
    omega_range = np.linspace(50, 400, 30)  # RPM
    damping_range = np.linspace(500, 8000, 20)  # Nms/rad

    # Run sweep
    analysis.sweep(
        solver_type="LTI",
        omega_values_RPM=omega_range,
        parameter_name="nominal_damping_Cd",
        parameter_values=damping_range,
        verbose=True
    )

    # Plot results
    fig1 = analysis.plot_stability_map()
    fig2 = analysis.plot_damping_slices()

    plt.show()
