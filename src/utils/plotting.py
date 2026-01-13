"""Plotting utilities and configuration.

Converted from MATLAB plot_properties.m and my_plot.m.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Optional, List, Tuple


@dataclass
class PlotProperties:
    """Plot styling properties."""
    
    marker_size: int = 4
    line_width: float = 2.0
    fontsize_legend: int = 11
    dash_width: float = 1.5
    fontsize_label: int = 15


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

    def __init__(self):
        """Initialize plotting utilities."""
        self.properties = PlotProperties()
        self.colors = ColorScheme()

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

        # Plot all eigenvalues for all operating points (dots only, no lines)
        n_modes = len(modal_solution[0].damping)
        for j in range(n_modes):
            rpm_values = [modal_sol.OMEGA_RPM for modal_sol in modal_solution]
            damping_values = [modal_sol.damping[j] for modal_sol in modal_solution]
            ax.plot(rpm_values, damping_values,
                   color='blue', marker='*', markersize=plot_property.marker_size,
                   linestyle='None')

        ax.set_xlabel(r'$\Omega$ [rpm]', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'$\lambda$ [-]', fontsize=plot_property.fontsize_label)

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
                       color=colors_double[j], marker='o',
                       markersize=plot_property.marker_size,
                       linewidth=plot_property.line_width)
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
                           color=colors_double[mode_idx % len(colors_double)], marker='o',
                           markersize=plot_property.marker_size,
                           linewidth=plot_property.line_width,
                           label=f'{mode_number // 2}')
                else:
                    ax.plot(rpm_values, damping_values,
                           color=colors_double[mode_idx % len(colors_double)], marker='o',
                           markersize=plot_property.marker_size,
                           linewidth=plot_property.line_width)

        ax.set_xlabel(r'$\Omega$ [rpm]', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'$\lambda$ [-]', fontsize=plot_property.fontsize_label)

        if modes is not None:
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

        # Plot all eigenvalues for all operating points (dots only, no lines)
        n_modes = len(modal_solution[0].frequency)
        for j in range(n_modes):
            rpm_values = [modal_sol.OMEGA_RPM for modal_sol in modal_solution]
            freq_values = [modal_sol.frequency[j] for modal_sol in modal_solution]
            ax.plot(rpm_values, freq_values,
                   color='blue', marker='*', markersize=plot_property.marker_size,
                   linestyle='None')

        ax.set_xlabel(r'$\Omega$ [rpm]', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'$\omega$ [rad/s]', fontsize=plot_property.fontsize_label)

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
                       color=colors_double[j], marker='o',
                       markersize=plot_property.marker_size,
                       linewidth=plot_property.line_width)
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
                           color=colors_double[mode_idx % len(colors_double)], marker='o',
                           markersize=plot_property.marker_size,
                           linewidth=plot_property.line_width,
                           label=f'{mode_number // 2}')
                else:
                    ax.plot(rpm_values, freq_values,
                           color=colors_double[mode_idx % len(colors_double)], marker='o',
                           markersize=plot_property.marker_size,
                           linewidth=plot_property.line_width)

        ax.set_xlabel(r'$\Omega$ [rpm]', fontsize=plot_property.fontsize_label)
        ax.set_ylabel(r'$\omega$ [rad/s]', fontsize=plot_property.fontsize_label)

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
        xlimits: Optional[Tuple[float, float]] = (0, 400),
        ylimits: Optional[Tuple[float, float]] = (0, 1),
        show_legend: bool = True,
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot modal participation factors for individual harmonics.

        Plots participation factors vs RPM for each harmonic component separately.
        Harmonics are styled as: solid lines (negative), dash-dot (zero), dashed (positive).

        Args:
            modal_participation: List of ModalParticipation objects
            state_index: Index of state to plot (0-based)
            mode_index: Index of mode to plot (0-based)
            xlimits: Optional x-axis limits (min_rpm, max_rpm)
            ylimits: Optional y-axis limits
            show_legend: Whether to show legend
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

        # Extract data for the specified state and mode
        rpm_values = np.array([mp.OMEGA_RPM for mp in modal_participation])
        n_harmonics = modal_participation[0].phi_jkn.shape[2]

        # Extract participation factors: phi_1mode1state[n, rpm]
        phi_1mode1state = np.zeros((n_harmonics, len(rpm_values)))
        for i, mp in enumerate(modal_participation):
            phi_1mode1state[:, i] = mp.phi_jkn[state_index, mode_index, :]

        # Get colors (single vector without duplication)
        colors_single = color.get_color_vector()

        # Determine middle index (zero harmonic)
        mid_idx = n_harmonics // 2

        # Create harmonic labels and plot
        n_harm_side = n_harmonics // 2

        # Plot each harmonic with consistent color and appropriate line style
        for i in range(n_harmonics):
            # Calculate harmonic number
            n_harmonic = i - mid_idx  # e.g., -4, -3, -2, -1, 0, 1, 2, 3, 4

            # Create label
            if n_harmonic == 0:
                label = '$n=0$'
                linestyle = '-.'  # Dash-dot for zero harmonic
            elif n_harmonic < 0:
                label = f'$n={n_harmonic:d}$'
                linestyle = '-'   # Solid for negative harmonics
            else:
                label = f'$n=+{n_harmonic:d}$'
                linestyle = '--'  # Dashed for positive harmonics

            # Use consistent color based on harmonic index
            color_idx = i % len(colors_single)

            ax.plot(rpm_values, phi_1mode1state[i, :],
                   linestyle=linestyle,
                   linewidth=plot_property.line_width,
                   color=colors_single[color_idx],
                   label=label)

        ax.set_xlabel(r'$\Omega$ [rpm]', fontsize=plot_property.fontsize_label)
        # Create clear y-label with state, mode indices
        ax.set_ylabel(rf'$\phi_{{{state_index},{mode_index},n}}$', fontsize=18)

        if show_legend:
            ax.legend(fontsize=plot_property.fontsize_legend, loc='center left',
                     bbox_to_anchor=(1, 0.5))

        if xlimits:
            ax.set_xlim(xlimits)
        if ylimits:
            ax.set_ylim(ylimits)

        fig.tight_layout()
        return fig

    @staticmethod
    def plot_mod_part_sum(
        modal_participation: List,
        state_index: int,
        mode_index: int,
        xlimits: Optional[Tuple[float, float]] = (0, 400),
        ylimits: Optional[Tuple[float, float]] = (0, 1),
        show_legend: bool = True,
        figure_handle: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """Plot summed modal participation factors (±n harmonics combined).

        Plots participation factors vs RPM with positive and negative harmonics summed.
        For example, n=±2 combines the contributions of n=-2 and n=+2.

        Args:
            modal_participation: List of ModalParticipation objects
            state_index: Index of state to plot (0-based)
            mode_index: Index of mode to plot (0-based)
            xlimits: Optional x-axis limits (min_rpm, max_rpm)
            ylimits: Optional y-axis limits
            show_legend: Whether to show legend
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

        # Extract data for the specified state and mode
        rpm_values = np.array([mp.OMEGA_RPM for mp in modal_participation])
        n_harmonics = modal_participation[0].phi_jkn.shape[2]

        # Extract participation factors
        phi_1mode1state = np.zeros((n_harmonics, len(rpm_values)))
        for i, mp in enumerate(modal_participation):
            phi_1mode1state[:, i] = mp.phi_jkn[state_index, mode_index, :]

        # Sum positive and negative harmonics
        mid_idx = n_harmonics // 2
        n_summed = mid_idx + 1  # Includes zero harmonic
        phi_1mode1state_sum = np.zeros((n_summed, len(rpm_values)))

        k = 0
        for j in range(mid_idx + 1):
            # Sum symmetric harmonics: phi[j] + phi[n_harmonics-1-j]
            if j == mid_idx:
                # Zero harmonic (no pair to sum)
                phi_1mode1state_sum[k, :] = phi_1mode1state[j, :]
            else:
                phi_1mode1state_sum[k, :] = (phi_1mode1state[j, :] +
                                             phi_1mode1state[n_harmonics - 1 - j, :])
            k += 1

        # Get colors (single vector without duplication)
        colors_single = color.get_color_vector()

        # Create labels and plot summed harmonics
        n_harm_side = n_harmonics // 2

        # Plot each summed harmonic with consistent color and style
        for i in range(n_summed):
            # Calculate harmonic number (counting down from max to 0)
            n_harmonic = n_harm_side - i

            # Determine label and linestyle
            if n_harmonic == 0:
                label = '$n=0$'
                linestyle = '-.'  # Dash-dot for zero harmonic
            else:
                label = f'$n=\\pm {n_harmonic}$'
                linestyle = '-'   # Solid for ±n harmonics

            # Use consistent color
            color_idx = i % len(colors_single)

            ax.plot(rpm_values, phi_1mode1state_sum[i, :],
                   linestyle=linestyle,
                   linewidth=plot_property.line_width,
                   color=colors_single[color_idx],
                   label=label)

        ax.set_xlabel(r'$\Omega$ [rpm]', fontsize=plot_property.fontsize_label)
        # Create clear y-label showing it's the sum of ±n harmonics
        ax.set_ylabel(rf'$\phi_{{{state_index},{mode_index},\pm n}}$', fontsize=18)

        if show_legend:
            ax.legend(fontsize=plot_property.fontsize_legend, loc='center left',
                     bbox_to_anchor=(1, 0.5))

        if xlimits:
            ax.set_xlim(xlimits)
        if ylimits:
            ax.set_ylim(ylimits)

        fig.tight_layout()
        return fig


# Initialize global instances for convenience
plot_property = PlotProperties()
color = ColorScheme()


if __name__ == "__main__":
    # Test color scheme
    colors = color.get_color_vector()
    print(f"Number of colors: {len(colors)}")
    print(f"First color (blue): {colors[0]}")
    
    # Test plotting utilities
    plotter = MyPlot()
    print(f"Line width: {plotter.properties.line_width}")
