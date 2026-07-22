"""Plotting utilities and configuration.

Converted from MATLAB plot_properties.m and my_plot.m.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from dataclasses import dataclass
from typing import Optional, List, Tuple

# Configure matplotlib to use LaTeX-style fonts (Computer Modern)
# and publication-quality settings
plt.rcParams.update({
    'text.usetex': False,  # Set True if LaTeX is installed for full LaTeX rendering
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',  # Computer Modern for math text
    'axes.unicode_minus': False,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
})


@dataclass
class PlotProperties:
    """Plot styling properties for publication-quality figures."""

    marker_size: int = 6
    line_width: float = 3.0
    fontsize_legend: int = 22
    dash_width: float = 2.0
    fontsize_label: int = 26


@dataclass
class ColorScheme:
    """Color scheme for plotting, matching MATLAB defaults."""
    
    # Original MATLAB default colors
    blue_mat: Tuple[float, float, float] = (0.0, 0.4470, 0.7410)
    orange_mat: Tuple[float, float, float] = (0.8500, 0.3250, 0.0980)
    yellow_mat: Tuple[float, float, float] = (0.9290, 0.6940, 0.1250)
    purple_mat: Tuple[float, float, float] = (0.4940, 0.1840, 0.5560)
    green_mat: Tuple[float, float, float] = (0.4660, 0.6740, 0.1880)
    light_blue_mat: Tuple[float, float, float] = (0.3010, 0.7450, 0.9330)
    red_mat: Tuple[float, float, float] = (0.6350, 0.0780, 0.1840)
    
    # Additional distinct colors
    teal_mat: Tuple[float, float, float] = (0.0, 0.5, 0.5)
    pink_mat: Tuple[float, float, float] = (0.9, 0.4, 0.7)
    brown_mat: Tuple[float, float, float] = (0.6, 0.3, 0.0)
    gray_mat: Tuple[float, float, float] = (0.5, 0.5, 0.5)
    cyan_mat: Tuple[float, float, float] = (0.0, 0.8, 0.8)
    magenta_mat: Tuple[float, float, float] = (0.8, 0.0, 0.8)
    lime_mat: Tuple[float, float, float] = (0.7, 0.9, 0.1)
    gold_mat: Tuple[float, float, float] = (1.0, 0.85, 0.0)
    navy_mat: Tuple[float, float, float] = (0.0, 0.0, 0.5)
    maroon_mat: Tuple[float, float, float] = (0.5, 0.0, 0.0)
    turquoise_mat: Tuple[float, float, float] = (0.25, 0.88, 0.815)
    violet_mat: Tuple[float, float, float] = (0.58, 0.0, 0.83)
    
    def get_color_vector(self) -> np.ndarray:
        """Get array of all colors."""
        return np.array([
            self.blue_mat,
            self.orange_mat,
            self.yellow_mat,
            self.green_mat,
            self.red_mat,
            (0, 0, 0),  # Black
            self.purple_mat,
            self.light_blue_mat,
            self.teal_mat,
            self.pink_mat,
            self.brown_mat,
            self.gray_mat,
            self.cyan_mat,
            self.magenta_mat,
            self.lime_mat,
            self.gold_mat,
            self.navy_mat,
            self.maroon_mat,
            self.turquoise_mat,
            self.violet_mat,
        ])
    
    def get_color_vector_double(self) -> np.ndarray:
        """Get doubled color vector for paired plotting."""
        colors = self.get_color_vector()
        return np.repeat(colors, 2, axis=0)


class MyPlot:
    """Plotting utilities for rotor dynamics analysis.

    This class will contain methods converted from my_plot.m
    """

    @staticmethod
    def plot_damping_generic(
        modal_solution: List,
        xlimits: Optional[Tuple[float, float]] = (0, 400),
        ylimits: Optional[Tuple[float, float]] = (-5, 1),
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot generic damping vs RPM.

        Args:
            modal_solution: List of modal solution structures
            xlimits: Optional x-axis limits (min_rpm, max_rpm)
            ylimits: Optional y-axis limits (min_damping, max_damping)
            figure_handle: Optional existing figure to plot on

        Returns:
            matplotlib Figure object
        """
        # Create or reuse figure
        if figure_handle is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = figure_handle
            ax = fig.gca()

        ax.grid(True)

        # Plot all eigenvalues for all operating points (continuous lines)
        n_modes = len(modal_solution[0].damping)
        for j in range(n_modes):
            rpm_values = [modal_sol.OMEGA_RPM for modal_sol in modal_solution]
            damping_values = [modal_sol.damping[j] for modal_sol in modal_solution]
            ax.plot(rpm_values, damping_values,
                   color='blue', marker='.', linestyle='none', markersize=plot_property.marker_size, linewidth=plot_property.line_width)

        ax.set_xlabel(r'$\Omega$ $\mathrm{[rpm]}$', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'$\lambda$ $\mathrm{[1/s]}$', fontsize=plot_property.fontsize_label)

        if xlimits:
            ax.set_xlim(xlimits)
        if ylimits:
            ax.set_ylim(ylimits)

        return fig

    @staticmethod
    def plot_damping_order(
        modal_solution: List,
        modes: Optional[List[int]] = None,
        xlimits: Optional[Tuple[float, float]] = (0, 400),
        ylimits: Optional[Tuple[float, float]] = (-5, 1),
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot damping for specific mode orders.

        Args:
            modal_solution: List of modal solution structures
            modes: List of mode numbers to plot (1-based indexing)
            xlimits: Optional x-axis limits (min_rpm, max_rpm)
            ylimits: Optional y-axis limits (min_damping, max_damping)
            figure_handle: Optional existing figure to plot on

        Returns:
            matplotlib Figure object
        """
        if figure_handle is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = figure_handle
            ax = fig.gca()

        ax.grid(True)

        # Get color vector
        colors_double = color.get_color_vector_double()

        # Generate ordered mode names
        if modes is not None and len(modes) > 0:
            min_mode_number = min(modes)
            max_mode_number = max(modes)
            mode_names = np.arange(min_mode_number, max_mode_number + 1)
        else:
            modes = None

        # Plot data depending on mode selection
        if modes is None:
            # Plot all modes
            n_modes = len(modal_solution[0].damping)
            for j in range(n_modes):
                rpm_values = [modal_sol.OMEGA_RPM for modal_sol in modal_solution]
                damping_values = [modal_sol.damping[j] for modal_sol in modal_solution]
                ax.plot(rpm_values, damping_values,
                       color=colors_double[j], linewidth=plot_property.line_width)
        else:
            # Plot selected modes (convert 1-based to 0-based for array indexing)
            # Eigenvalues come in pairs: (1,2), (3,4), (5,6), (7,8), (9,10), (11,12)
            # Each pair represents one mode: mode 1, 2, 3, 4, 5, 6
            for j, mode_number in enumerate(modes):
                mode_idx = mode_number - 1  # Convert to 0-based index for array access
                rpm_values = [modal_sol.OMEGA_RPM for modal_sol in modal_solution]
                damping_values = [modal_sol.damping[mode_idx] for modal_sol in modal_solution]

                # Color based on actual mode_idx (not position j)
                # This ensures consistent coloring even with filtered lists
                # Show label for even mode_numbers (2, 4, 6, 8, 10, 12)
                if mode_number % 2 == 0:
                    ax.plot(rpm_values, damping_values,
                           color=colors_double[j], linewidth=plot_property.line_width,
                           label=f'{j // 2+1}')
                else:
                    ax.plot(rpm_values, damping_values,
                           color=colors_double[j], linewidth=plot_property.line_width)

        ax.set_xlabel(r'$\Omega$ $\mathrm{[rpm]}$', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'$\lambda$ $\mathrm{[1/s]}$', fontsize=plot_property.fontsize_label)

        if modes is not None:
            ax.legend(fontsize=plot_property.fontsize_legend, loc='center left',
                     bbox_to_anchor=(1, 0.5))

        if xlimits:
            ax.set_xlim(xlimits)
        if ylimits:
            ax.set_ylim(ylimits)

        return fig
    @staticmethod
    def plot_damping_hd(
        modal_solution: List,
        modes: Optional[List[int]] = None,
        xlimits: Optional[Tuple[float, float]] = (0, 400),
        ylimits: Optional[Tuple[float, float]] = (-5, 1),
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot damping for specific mode orders.

        Args:
            modal_solution: List of modal solution structures
            modes: List of mode numbers to highlight (1-based indexing)
            xlimits: Optional x-axis limits (min_rpm, max_rpm)
            ylimits: Optional y-axis limits (min_damping, max_damping)
            figure_handle: Optional existing figure to plot on

        Returns:
            matplotlib Figure object
        """
        if figure_handle is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = figure_handle
            ax = fig.gca()

        ax.grid(True)

        # Get color vector
        colors_double = color.get_color_vector_double()

        # Generate ordered mode names for the modes to highlight
        if modes is not None and len(modes) > 0:
            min_mode_number = min(modes)
            max_mode_number = max(modes)
            mode_names = np.arange(min_mode_number, max_mode_number + 1)
        else:
            modes = None

        # Plot all data in purple
        n_modes = len(modal_solution[0].damping)
        rpm_values = [modal_sol.OMEGA_RPM for modal_sol in modal_solution]
        damping_values = [modal_sol.damping[0] for modal_sol in modal_solution]
        ax.plot(rpm_values, damping_values,
                    color=colors_double[12], linewidth=plot_property.line_width,label="All HD")
        for j in range(n_modes):
            damping_values = [modal_sol.damping[j] for modal_sol in modal_solution]
            ax.plot(rpm_values, damping_values,
                       color=colors_double[12], linewidth=plot_property.line_width)
        
        # superimpose the highlighted modes from the list
        # take the part to plot the modes from above
        for j in range(0,10): #plot the extra modes, used 12 for the special case 
            mode_number= modes[j]-1 
            damping_values = [modal_sol.damping[mode_number] for modal_sol in modal_solution]
            if j % 2 == 0:
                ax.plot(rpm_values, damping_values,
                       color=colors_double[j], linewidth=plot_property.line_width,label=f"HD {j//2+1}") # remember to change for more than 4 blades
            else:
                ax.plot(rpm_values, damping_values,
                       color=colors_double[j], linewidth=plot_property.line_width)
        
        mode_number = modes[11]
        damping_values = [modal_sol.damping[mode_number] for modal_sol in modal_solution]
        ax.plot(rpm_values, damping_values,
                       color=colors_double[11], linewidth=plot_property.line_width,label=f"HD 6") # remember to change for more than 4 blades

        ax.set_xlabel(r'$\Omega$ $\mathrm{[rpm]}$', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'$\lambda$ $\mathrm{[1/s]}$', fontsize=plot_property.fontsize_label)

        # use the predefined legend template

        ax.legend(fontsize=plot_property.fontsize_legend, loc='center left',
                     bbox_to_anchor=(1, 0.5))

        if xlimits:
            ax.set_xlim(xlimits)
        if ylimits:
            ax.set_ylim(ylimits)

        return fig
    


    @staticmethod
    def plot_frequency_generic(
        modal_solution: List,
        xlimits: Optional[Tuple[float, float]] = (0, 400),
        ylimits: Optional[Tuple[float, float]] = (-5, 1),
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot generic frequency vs RPM.

        Args:
            modal_solution: List of modal solution structures
            xlimits: Optional x-axis limits (min_rpm, max_rpm)
            ylimits: Optional y-axis limits (min_freq, max_freq)
            figure_handle: Optional existing figure to plot on

        Returns:
            matplotlib Figure object
        """
        if figure_handle is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = figure_handle
            ax = fig.gca()

        ax.grid(True)

        # Plot all eigenvalues for all operating points (continuous lines)
        n_modes = len(modal_solution[0].frequency)
        for j in range(n_modes):
            rpm_values = [modal_sol.OMEGA_RPM for modal_sol in modal_solution]
            freq_values = [modal_sol.frequency[j] for modal_sol in modal_solution]
            ax.plot(rpm_values, freq_values,
                   color='blue', marker='.', linestyle='none', markersize=plot_property.marker_size)

        ax.set_xlabel(r'$\Omega$ $\mathrm{[rpm]}$', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'$\omega$ $\mathrm{[rad/s]}$', fontsize=plot_property.fontsize_label)

        if xlimits:
            ax.set_xlim(xlimits)
        if ylimits:
            ax.set_ylim(ylimits)

        return fig

    @staticmethod
    def plot_frequency_order(
        modal_solution: List,
        modes: Optional[List[int]] = None,
        xlimits: Optional[Tuple[float, float]] = (0, 400),
        ylimits: Optional[Tuple[float, float]] = (0, 30),
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot frequency for specific mode orders.

        Args:
            modal_solution: List of modal solution structures
            modes: List of mode numbers to plot (1-based indexing)
            xlimits: Optional x-axis limits (min_rpm, max_rpm)
            ylimits: Optional y-axis limits (min_freq, max_freq)
            figure_handle: Optional existing figure to plot on

        Returns:
            matplotlib Figure object
        """
        if figure_handle is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = figure_handle
            ax = fig.gca()

        ax.grid(True)

        # Get color vector
        colors_double = color.get_color_vector_double()

        # Generate ordered mode names
        if modes is not None and len(modes) > 0:
            min_mode_number = min(modes)
            max_mode_number = max(modes)
            mode_names = np.arange(min_mode_number, max_mode_number + 1)
        else:
            modes = None

        # Plot data depending on mode selection
        if modes is None:
            # Plot all modes
            n_modes = len(modal_solution[0].frequency)
            for j in range(n_modes):
                rpm_values = [modal_sol.OMEGA_RPM for modal_sol in modal_solution]
                freq_values = [modal_sol.frequency[j] for modal_sol in modal_solution]
                ax.plot(rpm_values, freq_values,
                       color=colors_double[j], linewidth=plot_property.line_width)
        else:
            # Plot selected modes (convert 1-based to 0-based for array indexing)
            # Eigenvalues come in pairs: (1,2), (3,4), (5,6), (7,8), (9,10), (11,12)
            # Each pair represents one mode: mode 1, 2, 3, 4, 5, 6
            for j, mode_number in enumerate(modes):
                mode_idx = mode_number - 1  # Convert to 0-based index for array access
                rpm_values = [modal_sol.OMEGA_RPM for modal_sol in modal_solution]
                freq_values = [modal_sol.frequency[mode_idx] for modal_sol in modal_solution]

                # Color based on actual mode_idx (not position j)
                # This ensures consistent coloring even with filtered lists
                # Show label for even mode_numbers (2, 4, 6, 8, 10, 12)
                if mode_number % 2 == 0:
                    ax.plot(rpm_values, freq_values,
                           color=colors_double[mode_idx], linewidth=plot_property.line_width,
                           label=f'{mode_number // 2}')
                else:
                    ax.plot(rpm_values, freq_values,
                           color=colors_double[mode_idx], linewidth=plot_property.line_width)

        ax.set_xlabel(r'$\Omega$ $\mathrm{[rpm]}$', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'$\omega$ $\mathrm{[rad/s]}$', fontsize=plot_property.fontsize_label)

        if modes is not None:
            ax.legend(fontsize=plot_property.fontsize_legend, loc='center left',
                     bbox_to_anchor=(1, 0.5))

        if xlimits:
            ax.set_xlim(xlimits)
        if ylimits:
            ax.set_ylim(ylimits)

        return fig

    @staticmethod
    def plot_mod_part(
        modal_participation: List,
        state_index: int,
        mode_index: int,
        labels: Optional[List[str]] = None,
        xlimits: Optional[Tuple[float, float]] = None,
        plot_sum: bool = True
    ) -> plt.Figure:
        """Plot modal participation factors.

        Args:
            modal_participation: Modal participation data
            state_index: Index of state to plot
            mode_index: Index of mode to plot
            labels: Optional labels for plot
            xlimits: Optional x-axis limits
            plot_sum: Whether to plot sum of participations

        Returns:
            matplotlib Figure object
        """
        fig, ax = plt.subplots(figsize=(10, 6))

        # TODO: Implement plotting logic from MATLAB version
        # This is a placeholder

        ax.set_xlabel(r'$\Omega$ $\mathrm{[rpm]}$', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'Modal Participation $\mathrm{[-]}$', fontsize=plot_property.fontsize_label)
        ax.grid(True)

        if xlimits:
            ax.set_xlim(xlimits)

        if labels:
            ax.legend(labels)

        return fig

    @staticmethod
    def plot_parametric_stability_map(
        parametric_analysis,
        figsize: Tuple[float, float] = (10, 8),
        cmap: str = "RdYlGn_r",
        levels: int = 50,
        show_stability_boundary: bool = True,
        reference_analysis=None,
        reference_label: Optional[str] = None,
        reference_color: str = "blue",
        reference_linestyle: str = "-",
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot stability map from parametric analysis.

        Creates a filled contour plot showing damping (stability) as a function
        of rotor speed and the swept parameter. Red regions indicate instability
        (positive damping), green regions indicate stability (negative damping).

        Args:
            parametric_analysis: ParametricAnalysis object with completed sweep
            figsize: Figure size (width, height) in inches
            cmap: Colormap name. Default "RdYlGn_r" shows red=unstable, green=stable
            levels: Number of contour levels for smooth gradients
            show_stability_boundary: If True, draw dashed black line at damping=0
            reference_analysis: Optional second ParametricAnalysis to overlay its
                               stability boundary (e.g., LTP boundary on HD plot)
            reference_label: Label for the reference boundary in legend
            reference_color: Color for reference boundary line
            reference_linestyle: Line style for reference boundary
            title: Custom title. If None, auto-generated from solver type
            xlabel: Custom x-axis label. If None, uses "Rotor Speed [rpm]"
            ylabel: Custom y-axis label. If None, auto-generated from parameter name
            figure_handle: Optional existing figure to plot on

        Returns:
            matplotlib Figure object

        Example:
            >>> analysis_hd = ParametricAnalysis(problem)
            >>> analysis_hd.sweep(solver_type="HD", ...)
            >>> analysis_ltp = ParametricAnalysis(problem)
            >>> analysis_ltp.sweep(solver_type="LTP", ...)
            >>> fig = MyPlot.plot_parametric_stability_map(
            ...     analysis_hd,
            ...     reference_analysis=analysis_ltp,
            ...     reference_label="LTP boundary"
            ... )
        """
        if figure_handle is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = figure_handle
            ax = fig.gca()

        # Get data matrices
        damping = parametric_analysis.get_damping_matrix()
        omega_rpm = parametric_analysis.omega_values_RPM
        param_vals = parametric_analysis.parameter_values

        # Create meshgrid for contour plotting
        OMEGA_grid, PARAM_grid = np.meshgrid(omega_rpm, param_vals)

        # Filled contour plot
        cf = ax.contourf(OMEGA_grid, PARAM_grid, damping, levels=levels, cmap=cmap)
        cbar = fig.colorbar(cf, ax=ax)
        cbar.set_label(r'$\lambda$ $\mathrm{[1/s]}$ of least damped mode', fontsize=plot_property.fontsize_label)

        # Stability boundary (damping = 0) for main analysis
        if show_stability_boundary:
            cs = ax.contour(OMEGA_grid, PARAM_grid, damping, levels=[0],
                           colors='black', linewidths=2, linestyles='--')
            # Add to legend manually
            boundary_label = f"{parametric_analysis.solver_type} stability boundary"
            ax.plot([], [], 'k--', linewidth=plot_property.line_width, label=boundary_label)

        # Overlay reference analysis stability boundary if provided
        if reference_analysis is not None:
            ref_damping = reference_analysis.get_damping_matrix()
            ref_omega_rpm = reference_analysis.omega_values_RPM
            ref_param_vals = reference_analysis.parameter_values

            REF_OMEGA_grid, REF_PARAM_grid = np.meshgrid(ref_omega_rpm, ref_param_vals)

            cs_ref = ax.contour(REF_OMEGA_grid, REF_PARAM_grid, ref_damping, levels=[0],
                               colors=reference_color, linewidths=2.5,
                               linestyles=reference_linestyle)

            # Add to legend
            if reference_label is None:
                reference_label = f"{reference_analysis.solver_type} boundary"
            ax.plot([], [], color=reference_color, linestyle=reference_linestyle,
                   linewidth=plot_property.line_width, label=reference_label)

        # Labels
        if xlabel is None:
            xlabel = r'$\Omega$ $\mathrm{[rpm]}$'
        ax.set_xlabel(xlabel, fontsize=plot_property.fontsize_label)

        if ylabel is None:
            ylabel = MyPlot._format_parameter_label(parametric_analysis.parameter_name)
        ax.set_ylabel(ylabel, fontsize=plot_property.fontsize_label)

        if title is None:
            title = f"Stability Map ({parametric_analysis.solver_type} solver)"
        ax.set_title(title, fontsize=plot_property.fontsize_label + 2)

        ax.grid(True, alpha=0.3)

        # Show legend if we have boundaries
        if show_stability_boundary or reference_analysis is not None:
            ax.legend(loc='best', fontsize=plot_property.fontsize_legend)

        return fig

    @staticmethod
    def plot_parametric_frequency_map(
        parametric_analysis,
        figsize: Tuple[float, float] = (10, 8),
        cmap: str = "viridis",
        levels: int = 50,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot frequency map from parametric analysis.

        Creates a filled contour plot showing the frequency of the least damped
        mode as a function of rotor speed and the swept parameter.

        Args:
            parametric_analysis: ParametricAnalysis object with completed sweep
            figsize: Figure size (width, height) in inches
            cmap: Colormap name
            levels: Number of contour levels
            title: Custom title. If None, auto-generated from solver type
            xlabel: Custom x-axis label
            ylabel: Custom y-axis label
            figure_handle: Optional existing figure to plot on

        Returns:
            matplotlib Figure object
        """
        if figure_handle is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = figure_handle
            ax = fig.gca()

        # Get data matrices
        frequency = parametric_analysis.get_frequency_matrix()
        omega_rpm = parametric_analysis.omega_values_RPM
        param_vals = parametric_analysis.parameter_values

        # Create meshgrid
        OMEGA_grid, PARAM_grid = np.meshgrid(omega_rpm, param_vals)

        # Filled contour plot
        cf = ax.contourf(OMEGA_grid, PARAM_grid, frequency, levels=levels, cmap=cmap)
        cbar = fig.colorbar(cf, ax=ax)
        cbar.set_label(r'Frequency $\omega$ $\mathrm{[rad/s]}$', fontsize=plot_property.fontsize_label)

        # Labels
        if xlabel is None:
            xlabel = r'$\Omega$ $\mathrm{[rpm]}$'
        ax.set_xlabel(xlabel, fontsize=plot_property.fontsize_label)

        if ylabel is None:
            ylabel = MyPlot._format_parameter_label(parametric_analysis.parameter_name)
        ax.set_ylabel(ylabel, fontsize=plot_property.fontsize_label)

        if title is None:
            title = f"Frequency Map ({parametric_analysis.solver_type} solver)"
        ax.set_title(title, fontsize=plot_property.fontsize_label + 2)

        ax.grid(True, alpha=0.3)

        return fig

    @staticmethod
    def plot_parametric_damping_slices(
        parametric_analysis,
        param_indices: Optional[List[int]] = None,
        figsize: Tuple[float, float] = (10, 6),
        title: Optional[str] = None,
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot damping vs rotor speed for selected parameter values.

        Creates line plots showing how damping evolves with rotor speed
        for specific values of the swept parameter.

        Args:
            parametric_analysis: ParametricAnalysis object with completed sweep
            param_indices: List of parameter indices to plot. If None, plots
                          5 evenly spaced values across the parameter range.
            figsize: Figure size (width, height) in inches
            title: Custom title
            figure_handle: Optional existing figure to plot on

        Returns:
            matplotlib Figure object
        """
        if figure_handle is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = figure_handle
            ax = fig.gca()

        # Select parameter indices if not provided
        if param_indices is None:
            n = len(parametric_analysis.parameter_values)
            param_indices = list(np.linspace(0, n-1, min(5, n), dtype=int))

        colors_vec = color.get_color_vector()

        for i, idx in enumerate(param_indices):
            damping = [parametric_analysis.results[idx][j].damping
                      for j in range(len(parametric_analysis.omega_values))]
            param_val = parametric_analysis.parameter_values[idx]
            label = f"{parametric_analysis.parameter_name} = {param_val:.4g}"

            ax.plot(parametric_analysis.omega_values_RPM, damping,
                   color=colors_vec[i % len(colors_vec)],
                   linewidth=plot_property.line_width,
                   label=label)

        # Stability boundary
        ax.axhline(y=0, color='k', linestyle='--', linewidth=plot_property.dash_width,
                  label='Stability boundary')

        ax.set_xlabel(r'$\Omega$ $\mathrm{[rpm]}$', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'Damping $\lambda$ $\mathrm{[1/s]}$', fontsize=plot_property.fontsize_label)

        if title is None:
            title = f"Damping vs Rotor Speed ({parametric_analysis.solver_type} solver)"
        ax.set_title(title, fontsize=plot_property.fontsize_label + 2)

        ax.legend(fontsize=plot_property.fontsize_legend, loc='best')
        ax.grid(True, alpha=0.3)

        return fig

    @staticmethod
    def plot_parametric_stability_boundary(
        parametric_analysis,
        figsize: Tuple[float, float] = (10, 6),
        fill_unstable: bool = True,
        title: Optional[str] = None,
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot stability boundary curve.

        Extracts and plots the stability boundary (damping=0 contour) showing
        the critical parameter value as a function of rotor speed.

        Args:
            parametric_analysis: ParametricAnalysis object with completed sweep
            figsize: Figure size (width, height) in inches
            fill_unstable: If True, shade the unstable region
            title: Custom title
            figure_handle: Optional existing figure to plot on

        Returns:
            matplotlib Figure object
        """
        if figure_handle is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = figure_handle
            ax = fig.gca()

        damping = parametric_analysis.get_damping_matrix()
        omega_rpm = parametric_analysis.omega_values_RPM
        param_vals = parametric_analysis.parameter_values

        OMEGA_grid, PARAM_grid = np.meshgrid(omega_rpm, param_vals)

        # Extract stability boundary contour
        cs = ax.contour(OMEGA_grid, PARAM_grid, damping, levels=[0],
                       colors='black', linewidths=2)

        # Fill unstable region if requested
        if fill_unstable:
            ax.contourf(OMEGA_grid, PARAM_grid, damping, levels=[0, damping.max()],
                       colors=['red'], alpha=0.3)
            ax.contourf(OMEGA_grid, PARAM_grid, damping, levels=[damping.min(), 0],
                       colors=['green'], alpha=0.3)

        ax.set_xlabel(r'$\Omega$ $\mathrm{[rpm]}$', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(MyPlot._format_parameter_label(parametric_analysis.parameter_name),
                     fontsize=plot_property.fontsize_label)

        if title is None:
            title = f"Stability Boundary ({parametric_analysis.solver_type} solver)"
        ax.set_title(title, fontsize=plot_property.fontsize_label + 2)

        ax.grid(True, alpha=0.3)

        return fig

    @staticmethod
    def root_locus(
        modal_solution: List,
        figsize: Tuple[float, float] = (10, 6),
        title: Optional[str] = None,
        figure_handle: Optional[plt.Figure] = None
            ) -> plt.Figure:
        """Plot the locus of characteristic multipliers and the unit circle,"""

        # Create or reuse figure and axis
        if figure_handle is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = figure_handle
            ax = fig.gca()

        # Normalize input: accept
        # - a list of ModalSolution-like objects (with .char_solution or .eigensolution),
        # - a list/array of complex eigenvalue arrays,
        # - a single dict returned by LTPStability.ltp_single_point_given_M
        mult_list = []

        # Single dict case (char_solutions)
        if isinstance(modal_solution, dict):
            if 'char_multipliers' in modal_solution:
                vals = modal_solution.get('char_multipliers')
            elif 'eigenvalues' in modal_solution:
                vals = modal_solution.get('eigenvalues')
            else:
                vals = None
            if vals is None:
                raise ValueError('Dict input must contain "char_multipliers" or "eigenvalues"')
            mult_list.append(np.asarray(vals))

        # Numpy array of complex eigenvalues (single point)
        elif isinstance(modal_solution, np.ndarray) and modal_solution.ndim == 1:
            mult_list.append(np.asarray(modal_solution))

        else:
            # Assume iterable of solutions or arrays
            for sol in modal_solution:
                vals = None
                if isinstance(sol, dict):
                    if 'char_multipliers' in sol:
                        vals = sol.get('char_multipliers')
                    elif 'eigenvalues' in sol:
                        vals = sol.get('eigenvalues')
                    else:
                        vals = None
                else:
                    if hasattr(sol, 'char_solution') and sol.char_solution:
                        vals = sol.char_solution.get('char_multipliers')
                    if vals is None and hasattr(sol, 'eigensolution') and sol.eigensolution:
                        vals = sol.eigensolution.get('eigenvalues') or sol.eigensolution.get('eigvals')
                    # If sol itself is array-like, accept it
                    if vals is None:
                        try:
                            arr = np.asarray(sol)
                            if arr.ndim == 1 and np.iscomplexobj(arr):
                                vals = arr
                        except Exception:
                            vals = None

                if vals is None:
                    raise ValueError('Each modal_solution entry must contain char_multipliers/eigenvalues or be an array of complex eigenvalues')
                mult_list.append(np.asarray(vals))

        # Ensure consistent shape
        n_points = len(mult_list)
        n_modes = mult_list[0].size
        for arr in mult_list:
            if arr.size != n_modes:
                raise ValueError('Inconsistent number of modes across modal_solution entries')

        # Stack into (n_points, n_modes)
        mults = np.vstack(mult_list)

        # Plot shaded unstable region outside unit circle
        max_radius = max(1.0, np.max(np.abs(mults)) * 1.2)
        outer = Circle((0, 0), max_radius, color='lightcoral', alpha=0.12, zorder=0)
        inner = Circle((0, 0), 1.0, color=ax.get_facecolor(), zorder=1)
        ax.add_patch(outer)
        ax.add_patch(inner)

        # Unit circle boundary
        theta = np.linspace(0, 2 * np.pi, 400)
        ax.plot(np.cos(theta), np.sin(theta), color='red', linestyle='--', linewidth=2.0, zorder=2)

        # Choose colors for mode branches
        cmap = plt.get_cmap('tab20')

        # Plot each mode branch (connect points across rotor speeds)
        for j in range(n_modes):
            real_parts = np.real(mults[:, j])
            imag_parts = np.imag(mults[:, j])
            color_idx = j % cmap.N
            ax.plot(real_parts, imag_parts, '-', color=cmap(color_idx), linewidth=plot_property.line_width * 0.6, alpha=0.9, zorder=3)
            ax.plot(real_parts, imag_parts, '.', color=cmap(color_idx), markersize=plot_property.marker_size, zorder=4)

        ax.set_xlabel('Real', fontsize=plot_property.fontsize_label)
        ax.set_ylabel('Imag', fontsize=plot_property.fontsize_label)
        ax.set_title(title or 'Root Locus (Characteristic Multipliers)', fontsize=plot_property.fontsize_label + 2)
        ax.grid(True, which='both', linestyle=':', alpha=0.6)
        ax.set_aspect('equal', 'box')

        # Set limits with some margin
        lim = max_radius
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)

        return fig


    @staticmethod
    def _format_parameter_label(parameter_name: str) -> str:
        """Format parameter name for axis label.

        Args:
            parameter_name: Raw parameter name from analysis

        Returns:
            Formatted label with units if known
        """
        label_map = {
            "nominal_damping_Cd": r"Nominal Damping $C_d$ $\mathrm{[Nms/rad]}$",
            "nominal_stiffness_Kd": r"Damper Stiffness $K_d$ $\mathrm{[Nm/rad]}$",
            "hub_damping_Cx": r"Hub Damping $C_x$ $\mathrm{[Ns/m]}$",
            "hub_damping_Cy": r"Hub Damping $C_y$ $\mathrm{[Ns/m]}$",
            "hub_stiffness_Kx": r"Hub Stiffness $K_x$ $\mathrm{[N/m]}$",
            "hub_stiffness_Ky": r"Hub Stiffness $K_y$ $\mathrm{[N/m]}$",
            "blade_mass_Mb": r"Blade Mass $M_b$ $\mathrm{[kg]}$",
            "blade_inertia_Ib": r"Blade Inertia $I_b$ $\mathrm{[kg \cdot m^2]}$",
            "lag_hinge_offset_e": r"Lag Hinge Offset $e$ $\mathrm{[m]}$",
            "damper_1_ratio": r"First Damper Effectiveness $\varepsilon_1 = C_{d,1}/C_{d,nom}$ $\mathrm{[-]}$",
            "damper_ratio": r"Damper Effectiveness $\varepsilon = C_d/C_{d,nom}$ $\mathrm{[-]}$",
        }
        return label_map.get(parameter_name, parameter_name)


# Initialize global instances for convenience
plot_property = PlotProperties()
color = ColorScheme()


if __name__ == "__main__":
    # Test color scheme
    colors = color.get_color_vector()
    print(f"Number of colors: {len(colors)}")
    print(f"First color (blue): {colors[0]}")
    
    # Test plotting utilities
    print(f"Line width: {plot_property.line_width}")
