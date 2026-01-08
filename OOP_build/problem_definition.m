%% OVERALL PROBLEM DEFINITION (must be compatible with the existing system matrices)

% number of rotor blades (3, 4, 5, 7)
problem.number_blades = 4;
% type of blade damper connection ("H2B" or "B2B")
problem.damper_connection = "H2B";
% activation of damper ("ALL", "ODI", "2DI ADJ", "2DI OPP", "3DI")
problem.damper_activation = "ODI";
% solution interval in RPM
problem.lower_rotor_RPM = 20;
problem.higher_rotor_RPM = 400;
problem.number_points = 50;
% type of solver ("LTI", "LTP", "HD")
problem.required_solver = "HD";
% definition of the number of harmonics in the HD, equal to the number of
% harmonics in the modal participation
problem.number_harmonics = 1;

%% CONTINUATION METHOD
problem.continuation = "NO";
problem.continuation_tolerance = 1e-2;
problem.continuation_max_iter = 100;
problem.step_h=1e-4;
% the following parameter is for HD only
problem.time_samples = 100;
% how the continuation solver moves along the interval ("L2R" or "R2L")
problem.solution_direction = "R2L";
%% acquisition of rotor characteristics
run("rotor_definition_script.m")





