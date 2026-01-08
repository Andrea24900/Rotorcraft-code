%% DO NOT CHANGE HERE
problem.rotor_characteristics.isotropic_hub=0; % deprecated, do not change


%% CHANGEABLE CHARACTERISTICS (CHANGE NUMBERS ONLY)
problem.rotor_characteristics.blade_mass_Mb=94.9; %[kg]
problem.rotor_characteristics.blade_static_Sb = 289.1; %[kg m]
problem.rotor_characteristics.blade_inertia_Ib = 1085; %[kg m^2]
problem.rotor_characteristics.lag_hinge_offset_e = 0.3048; %[m]

problem.rotor_characteristics.hub_mass_Mx = 8027; %[kg]

if problem.rotor_characteristics.isotropic_hub==0
    problem.rotor_characteristics.hub_mass_My = 3284; %[kg]
else
    problem.rotor_characteristics.hub_mass_My = problem.rotor_characteristics.hub_mass_Mx;
end

problem.rotor_characteristics.hub_damping_Cx = 51420; %[Ns/m]
if problem.rotor_characteristics.isotropic_hub==0
    problem.rotor_characteristics.hub_damping_Cy = 25710; %[kg]
else
    problem.rotor_characteristics.hub_damping_Cy = problem.rotor_characteristics.hub_damping_Cx;
end

problem.rotor_characteristics.hub_stiffness_Kx = 1241000; %[N/m]
problem.rotor_characteristics.hub_stiffness_Ky = problem.rotor_characteristics.hub_stiffness_Kx;

%% this script contains other damper values for the two configurations


if problem.damper_connection == "H2B"
    %% HUB2BLADE
    problem.rotor_characteristics.nominal_damping_Cd = 4067; %[Nms/rad] % nominal 4067
    problem.rotor_characteristics.nominal_stiffness_Kd = 0; %[Nm/rad]
elseif problem.damper_connection == "B2B"
    %% BLADE2BLADE
    problem.rotor_characteristics.nominal_damping_Cd =   36500;    %[Ns/m]
    problem.rotor_characteristics.nominal_stiffness_Kd = 0; %[N/m]

    %this script contains geometric relationships for the B2B rotor
    problem.rotor_characteristics.delta_psi=2*pi/problem.number_blades;
    problem.rotor_characteristics.phi = (pi-problem.rotor_characteristics.delta_psi)/2;

    % GEOMETRIC CHARACTERISATION (may be changed)
    problem.rotor_characteristics.hub_to_A_a = 0.46; %[m]
    problem.rotor_characteristics.hub_to_B_b = 0.23; %[m]
    problem.rotor_characteristics.hinge_to_hinge_c = 2*problem.rotor_characteristics.lag_hinge_offset_e...
        *sin(problem.rotor_characteristics.delta_psi/2);
    problem.rotor_characteristics.length_undef_l0 = sqrt(problem.rotor_characteristics.hub_to_B_b^2 ...
        +problem.rotor_characteristics.hub_to_A_a^2);

    % equilibrium variables (required only in this formulation)
    % SET
    problem.rotor_characteristics.zeta_E = 0; %[rad]
    % computed
    problem.rotor_characteristics.length_eq_lE = ...
        sqrt((problem.rotor_characteristics.hinge_to_hinge_c+...
        problem.rotor_characteristics.hub_to_A_a*cos(problem.rotor_characteristics.zeta_E-...
        problem.rotor_characteristics.phi)+problem.rotor_characteristics.hub_to_B_b*...
        cos(problem.rotor_characteristics.zeta_E+problem.rotor_characteristics.phi))^2+...
    (problem.rotor_characteristics.hub_to_A_a*sin(problem.rotor_characteristics.zeta_E-problem.rotor_characteristics.phi)...
    +problem.rotor_characteristics.hub_to_B_b*sin(problem.rotor_characteristics.zeta_E+problem.rotor_characteristics.phi))^2);
    problem.rotor_characteristics.gamma_E = atan2 ( - ( problem.rotor_characteristics.hub_to_A_a*sin(problem.rotor_characteristics.zeta_E...
        -problem.rotor_characteristics.phi)+problem.rotor_characteristics.hub_to_B_b*sin(problem.rotor_characteristics.zeta_E+problem.rotor_characteristics.phi)),...
    (problem.rotor_characteristics.hinge_to_hinge_c+problem.rotor_characteristics.hub_to_A_a*cos(problem.rotor_characteristics.zeta_E-problem.rotor_characteristics.phi)...
    +problem.rotor_characteristics.hub_to_B_b*cos(problem.rotor_characteristics.zeta_E+problem.rotor_characteristics.phi)));

end







