classdef my_plot

    properties
        rotor_build
        data_to_plot
        plot_properties
        colors
    end


    methods (Static)

        function figure_handle = plot_damping_generic(data_to_plot, varargin)
            % Load plot properties
            run("plot_properties.m")

            % Create input parser
            p = inputParser;

            % Define optional parameters with defaults
            addParameter(p, 'XLimits', [0 400], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'YLimits', [-5 1], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'FigureHandle', [], @(h) isempty(h) || (ishandle(h) && strcmp(get(h, 'Type'), 'figure')));

            % Parse inputs
            parse(p, varargin{:});

            % Extract parameters
            x_limits = p.Results.XLimits;
            y_limits = p.Results.YLimits;
            figure_handle = p.Results.FigureHandle;

            % Handle figure creation or reuse
            if isempty(figure_handle)
                figure_handle = figure;
            else
                figure(figure_handle); % bring existing figure to focus
            end

            hold on
            grid on

            % Plotting loop
            for i = 1:width(data_to_plot)
                plot(data_to_plot(i).OMEGA_RPM, data_to_plot(i).damping, ...
                    'Color', 'blue', 'Marker', '*', 'MarkerSize', plot_property.marker_size, ...
                    'LineWidth', plot_property.line_width,'HandleVisibility','off');
            end

            % Stability line
            %yline(0, '--r', 'Upper stability bound', ...
             %   'LineWidth', plot_property.line_width, ...
             %  'FontSize', plot_property.fontsize_label, ...
             %  'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex','HandleVisibility','off');

            xlabel('$\Omega$ [rpm]', 'Interpreter', 'latex', 'FontSize', plot_property.fontsize_label);
            ylabel('$\lambda$ [-]', 'Interpreter', 'latex', 'FontSize', plot_property.fontsize_label);

            % Apply axis limits
            xlim(x_limits);
            ylim(y_limits);
        end


        function figure_handle = plot_frequency_generic(data_to_plot,varargin)
            % Load plot properties
            run("plot_properties.m")
            % Create input parser
            p = inputParser;

            % Define optional parameters with defaults
            addParameter(p, 'XLimits', [0 400], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'YLimits', [-5 1], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'FigureHandle', [], @(h) isempty(h) || (ishandle(h) && strcmp(get(h, 'Type'), 'figure')));

            % Parse inputs
            parse(p, varargin{:});

            % Extract parameters
            x_limits = p.Results.XLimits;
            y_limits = p.Results.YLimits;
            figure_handle = p.Results.FigureHandle;

            % Handle figure creation or reuse
            if isempty(figure_handle)
                figure_handle = figure;
            else
                figure(figure_handle); % bring existing figure to focus
            end
            hold on
            grid on
            for i = 1:width(data_to_plot)
                plot(data_to_plot(i).OMEGA_RPM,data_to_plot(i).frequency,...
                    "Color",'blue','Marker','*','MarkerSize',plot_property.marker_size,...
                    'LineWidth',plot_property.line_width,'HandleVisibility','off');
            end
            xlabel('$\Omega$ [rpm]','Interpreter','latex','FontSize',plot_property.fontsize_label)
            ylabel('$\omega$ [rad/s]','Interpreter','latex','FontSize',plot_property.fontsize_label)
            % Apply axis limits
            xlim(x_limits);
            ylim(y_limits);
        end

        function figure_handle = plot_damping_order(data_to_plot, varargin)
            run("plot_properties.m")

            % Set up input parser
            p = inputParser;
            addParameter(p, 'XLimits', [0 400], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'YLimits', [-5 1], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'Modes', [], @(x) isempty(x) || isnumeric(x));
            addParameter(p, 'FigureHandle', [], @(h) isempty(h) || (ishandle(h) && strcmp(get(h, 'Type'), 'figure')));

            parse(p, varargin{:});

            x_limits = p.Results.XLimits;
            y_limits = p.Results.YLimits;
            figure_handle = p.Results.FigureHandle;
            modes_to_plot = p.Results.Modes;

            % Create or reuse figure
            if isempty(figure_handle)
                figure_handle = figure;
            else
                figure(figure_handle);
            end

            hold on
            grid on

            %generation of ordered mode names
            min_mode_number = min(modes_to_plot);
            max_mode_number = max(modes_to_plot);
            mode_names = (min_mode_number:max_mode_number);
            % Plot data depending on mode selection
            if isempty(modes_to_plot)
                for j = 1:height(data_to_plot(1).damping)
                    for i = 1:width(data_to_plot)
                        plot(data_to_plot(i).OMEGA_RPM, data_to_plot(i).damping(j), ...
                            "Color", color.color_vect_double(j,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, 'HandleVisibility', 'off');
                    end
                end
            else
                for j = 1:length(modes_to_plot)
                    for i = 1:width(data_to_plot) - 1
                        plot(data_to_plot(i).OMEGA_RPM, data_to_plot(i).damping(modes_to_plot(j)), ...
                            "Color", color.color_vect_double(j,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, 'HandleVisibility', 'off');
                    end
                    % Decide legend visibility based on mode index
                    if mod(j,2) == 0

                        plot(data_to_plot(end).OMEGA_RPM, data_to_plot(end).damping(modes_to_plot(j)), ...
                            "Color", color.color_vect_double(j,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, ...
                            'DisplayName', sprintf('%d', mode_names(j)/2));
                    else
                        plot(data_to_plot(end).OMEGA_RPM, data_to_plot(end).damping(modes_to_plot(j)), ...
                            "Color", color.color_vect_double(j,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, 'HandleVisibility', 'off');
                    end
                end
            end

            % Stability line
            %yline(0, '--r', 'Upper stability bound', ...
            %    'LineWidth', plot_property.line_width, ...
            %    'FontSize', plot_property.fontsize_label, ...
             %   'LabelHorizontalAlignment', 'left', ...
              %  'Interpreter', 'latex', 'HandleVisibility', 'off');

            % Axes labels
            xlabel('$\Omega$ [rpm]', 'Interpreter', 'latex', 'FontSize', plot_property.fontsize_label);
            ylabel('$\lambda$ [-]', 'Interpreter', 'latex', 'FontSize', plot_property.fontsize_label);

            % Legend
            legend('Interpreter', 'latex', 'FontSize', plot_property.fontsize_legend, 'Location', 'eastoutside');

            % Axis limits
            xlim(x_limits);
            ylim(y_limits);
        end


        function figure_handle = plot_frequency_order(data_to_plot, varargin)
            run("plot_properties.m")

            % Input parser
            p = inputParser;
            addParameter(p, 'XLimits', [0 400], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'YLimits', [0 30], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'FigureHandle', [], @(h) isempty(h) || (ishandle(h) && strcmp(get(h, 'Type'), 'figure')));
            addParameter(p, 'ModesToPlot', [], @(x) isempty(x) || isnumeric(x));
            parse(p, varargin{:});

            x_limits = p.Results.XLimits;
            y_limits = p.Results.YLimits;
            figure_handle = p.Results.FigureHandle;
            modes_to_plot = p.Results.ModesToPlot;

            % Create or reuse figure
            if isempty(figure_handle)
                figure_handle = figure;
            else
                figure(figure_handle);
            end

            hold on
            grid on
            %generation of ordered mode names
            min_mode_number = min(modes_to_plot);
            max_mode_number = max(modes_to_plot);
            mode_names = (min_mode_number:max_mode_number);


            if isempty(modes_to_plot)
                for j = 1:height(data_to_plot(1).frequency)
                    for i = 1:width(data_to_plot) - 1
                        plot(data_to_plot(i).OMEGA_RPM, data_to_plot(i).frequency(j), ...
                            "Color", color.color_vect_double(j,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, 'HandleVisibility', 'off');
                    end
                end
            else
                for j = 1:length(modes_to_plot)
                    for i = 1:width(data_to_plot)
                        plot(data_to_plot(i).OMEGA_RPM, data_to_plot(i).frequency(modes_to_plot(j)), ...
                            "Color", color.color_vect_double(j,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, 'HandleVisibility', 'off');
                    end

                    % Legend only on even-indexed modes
                    if mod(j, 2) == 0
                        plot(data_to_plot(end).OMEGA_RPM, data_to_plot(end).frequency(modes_to_plot(j)), ...
                            "Color", color.color_vect_double(j,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, ...
                            'DisplayName', sprintf('%d', mode_names(j)/2));
                    else
                        plot(data_to_plot(end).OMEGA_RPM, data_to_plot(end).frequency(modes_to_plot(j)), ...
                            "Color", color.color_vect_double(j,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, 'HandleVisibility', 'off');
                    end
                end
            end

            % Axes labels and legend
            xlabel('$\Omega$ [rpm]', 'Interpreter', 'latex', 'FontSize', plot_property.fontsize_label);
            ylabel('$\omega$ [rad/s]', 'Interpreter', 'latex', 'FontSize', plot_property.fontsize_label);
            legend('Interpreter', 'latex', 'FontSize', plot_property.fontsize_legend, 'Location', 'eastoutside');

            % Apply axis limits
            xlim(x_limits);
            ylim(y_limits);
        end

                function figure_handle = plot_mod_part(data_to_plot,state,mode,mode_name, varargin)
            run("plot_properties.m")

            % Set up input parser
            p = inputParser;
            addParameter(p, 'XLimits', [0 400], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'YLimits', [0 1], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'PlotSum', false, @(x) islogical(x) || isnumeric(x));
            addParameter(p, 'FigureHandle', [], @(h) isempty(h) || (ishandle(h) && strcmp(get(h, 'Type'), 'figure')));

            parse(p, varargin{:});

            x_limits = p.Results.XLimits;
            y_limits = p.Results.YLimits;
            figure_handle = p.Results.FigureHandle;


            % Create or reuse figure
            if isempty(figure_handle)
                figure_handle = figure;
            else
                figure(figure_handle);
            end

            hold on
            grid on



            if ~p.Results.PlotSum
            % Plot data depending on mode selection
            for n=1:size(data_to_plot(1).phi_jkn,3)
                    n_legend = n-ceil(size(data_to_plot(1).phi_jkn,3)/2);
                    for i = 1:max(size(data_to_plot))-1
                        plot(data_to_plot(i).OMEGA_RPM, data_to_plot(i).phi_jkn(state,mode,n), ...
                            "Color", color.color_vect(n,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, ...
                            'HandleVisibility','off');
                    end
                    plot(data_to_plot(end).OMEGA_RPM, data_to_plot(end).phi_jkn(state,mode,n), ...
                            "Color", color.color_vect(n,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, ...
                            'DisplayName', sprintf('$n$ = %d',n_legend));
            end
            elseif p.Results.PlotSum

                for n=1:ceil(size(data_to_plot(1).phi_jkn,3)/2)
                    n_legend = -n+ceil(size(data_to_plot(1).phi_jkn,3)/2);
                    for i = 1:max(size(data_to_plot))-1
                        sum_to_plot =  data_to_plot(i).phi_jkn(state,mode,n) + ...
                            data_to_plot(i).phi_jkn(state,mode,end+1-n);
                        plot(data_to_plot(i).OMEGA_RPM,sum_to_plot, ...
                            "Color", color.color_vect(n,:), ...
                            'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                            'LineWidth', plot_property.line_width, ...
                            'HandleVisibility','off');
                    end

                    sum_to_plot =  data_to_plot(end).phi_jkn(state,mode,n) + ...
                            data_to_plot(end).phi_jkn(state,mode,end+1-n);
                    plot(data_to_plot(end).OMEGA_RPM, sum_to_plot, ...
                        "Color", color.color_vect(n,:), ...
                        'Marker', 'o', 'MarkerSize', plot_property.marker_size, ...
                        'LineWidth', plot_property.line_width, ...
                        'DisplayName', sprintf('$n$ = +/- %d',n_legend));
                end
            end


            % Axes labels
            xlabel('$\Omega$ [rpm] ', 'Interpreter', 'latex', 'FontSize', plot_property.fontsize_label);
            ylabel(mode_name, 'Interpreter', 'latex', 'FontSize', 1.5*plot_property.fontsize_label);

            % Legend
            legend('Interpreter', 'latex', 'FontSize', plot_property.fontsize_legend, 'Location', 'eastoutside');

            % Axis limits
            xlim(x_limits);
            ylim(y_limits);
        end


    end
end
