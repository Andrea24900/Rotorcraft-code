%% plot settings
plot_property.marker_size = 4;
plot_property.line_width = 2;
plot_property.fontsize_legend = 11;
plot_property.dash_width = 1.5;
plot_property.fontsize_label=15;


% Original MATLAB default colors
color.blue_mat       = [0 0.4470 0.7410];
color.orange_mat     = [0.8500 0.3250 0.0980];
color.yellow_mat     = [0.9290 0.6940 0.1250];
color.purple_mat     = [0.4940 0.1840 0.5560];
color.green_mat      = [0.4660 0.6740 0.1880];
color.light_blue_mat = [0.3010 0.7450 0.9330];
color.red_mat        = [0.6350 0.0780 0.1840];

% ✅ Additional distinct colors (well-spaced)
color.teal_mat       = [0.0000 0.5000 0.5000]; % Teal
color.pink_mat       = [0.9000 0.4000 0.7000]; % Pink
color.brown_mat      = [0.6000 0.3000 0.0000]; % Brown
color.gray_mat       = [0.5000 0.5000 0.5000]; % Gray
color.cyan_mat       = [0.0000 0.8000 0.8000]; % Cyan
color.magenta_mat    = [0.8000 0.0000 0.8000]; % Magenta
color.lime_mat       = [0.7000 0.9000 0.1000]; % Lime
color.gold_mat       = [1.0000 0.8500 0.0000]; % Gold
color.navy_mat       = [0.0000 0.0000 0.5000]; % Navy
color.maroon_mat     = [0.5000 0.0000 0.0000]; % Maroon
color.turquoise_mat  = [0.2500 0.8800 0.8150]; % Turquoise
color.violet_mat     = [0.5800 0.0000 0.8300]; % Violet

% ✅ Combine into color vectors
color.color_vect = [
    color.blue_mat;
    color.orange_mat;
    color.yellow_mat;
    color.green_mat;
    color.red_mat;
    [0 0 0]; % Black
    color.purple_mat;
    color.light_blue_mat;
    color.teal_mat;
    color.pink_mat;
    color.brown_mat;
    color.gray_mat;
    color.cyan_mat;
    color.magenta_mat;
    color.lime_mat;
    color.gold_mat;
    color.navy_mat;
    color.maroon_mat;
    color.turquoise_mat;
    color.violet_mat
];

% ✅ Double vector (pairs for plotting)
color.color_vect_double = reshape(repelem(color.color_vect(:), 2), [], 3);
