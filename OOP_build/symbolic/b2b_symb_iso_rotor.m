clear all 
close all
clc 



syms  i

rotor_OMEGA = sym('rotor_OMEGA','real');

t = sym('t','real');

azimuth_psi (t,i) = rotor_OMEGA*t+2*sym(pi)/4*(i-1);

T(t) = [1 cos(azimuth_psi(t,1)) sin(azimuth_psi(t,1)) -1;
    1 cos(azimuth_psi(t,2)) sin(azimuth_psi(t,2)) 1;
    1 cos(azimuth_psi(t,3)) sin(azimuth_psi(t,3)) -1;
    1 cos(azimuth_psi(t,4)) sin(azimuth_psi(t,4)) 1];


T_T = simplify(T');



T_dot = diff(T,t);

T_dotdot = diff(T_dot,t);

syms blade_mass_Mb hub_mass_Mx hub_mass_My blade_static_Sb blade_inertia_Ib number_blades_Nb

M_R = [blade_inertia_Ib 0 0 0 blade_static_Sb*sin(azimuth_psi(t,1)) -blade_static_Sb*cos(azimuth_psi(t,1));
    0 blade_inertia_Ib 0 0 blade_static_Sb*sin(azimuth_psi(t,2)) -blade_static_Sb*cos(azimuth_psi(t,2));
    0 0 blade_inertia_Ib 0 blade_static_Sb*sin(azimuth_psi(t,3)) -blade_static_Sb*cos(azimuth_psi(t,3));
    0 0 0 blade_inertia_Ib blade_static_Sb*sin(azimuth_psi(t,4)) -blade_static_Sb*cos(azimuth_psi(t,4));
    blade_static_Sb*sin(azimuth_psi(t,1)) blade_static_Sb*sin(azimuth_psi(t,2)) blade_static_Sb*sin(azimuth_psi(t,3)) blade_static_Sb*sin(azimuth_psi(t,4)) hub_mass_Mx+number_blades_Nb*blade_mass_Mb 0;
    -blade_static_Sb*cos(azimuth_psi(t,1)) -blade_static_Sb*cos(azimuth_psi(t,2)) -blade_static_Sb*cos(azimuth_psi(t,3)) -blade_static_Sb*cos(azimuth_psi(t,4)) 0 hub_mass_My+number_blades_Nb*blade_mass_Mb];

% damping matrix without blade-to-blade dampers
syms hub_damping_Cx hub_damping_Cy

C_R = [0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    blade_static_Sb*2*rotor_OMEGA*cos(azimuth_psi(t,1)) blade_static_Sb*2*rotor_OMEGA*cos(azimuth_psi(t,2)) blade_static_Sb*2*rotor_OMEGA*cos(azimuth_psi(t,3)) blade_static_Sb*2*rotor_OMEGA*cos(azimuth_psi(t,4)) hub_damping_Cx 0;
    blade_static_Sb*2*rotor_OMEGA*sin(azimuth_psi(t,1)) blade_static_Sb*2*rotor_OMEGA*sin(azimuth_psi(t,2)) blade_static_Sb*2*rotor_OMEGA*sin(azimuth_psi(t,3)) blade_static_Sb*2*rotor_OMEGA*sin(azimuth_psi(t,4)) 0 hub_damping_Cy];

% additional damping 
syms zeta_E gamma_E length_eq_lE phi length_undef_l0 damping_cd stiffness_kd hub_to_A_a hub_to_B_b 

C_zetad = hub_to_A_a^2*damping_cd*sin(zeta_E-phi+gamma_E)^2+hub_to_B_b^2*damping_cd*sin(zeta_E+phi+gamma_E)^2;

C_zetaed = hub_to_A_a*hub_to_B_b*damping_cd*sin(zeta_E-phi+gamma_E)*sin(zeta_E+phi+gamma_E);

syms C_zetad C_zetaed
C_R_ADD = [C_zetad C_zetaed 0 C_zetaed 0 0;
           C_zetaed C_zetad C_zetaed 0 0 0;
           0 C_zetaed C_zetad C_zetaed 0 0;
           C_zetaed 0 C_zetaed C_zetad 0 0;
           0 0 0 0 0 0;
           0 0 0 0 0 0];

C_R = C_R_ADD +C_R;

% stiffness matrix without blade-to-blade springs
syms hub_stiffness_Kx hub_stiffness_Ky lag_hinge_offset_e

K_R = [lag_hinge_offset_e*rotor_OMEGA^2*blade_static_Sb 0 0 0 0 0;
    0 lag_hinge_offset_e*rotor_OMEGA^2*blade_static_Sb 0 0 0 0;
    0 0 lag_hinge_offset_e*rotor_OMEGA^2*blade_static_Sb 0 0 0;
    0 0 0 lag_hinge_offset_e*rotor_OMEGA^2*blade_static_Sb 0 0;
    -blade_static_Sb*rotor_OMEGA^2*sin(azimuth_psi(t,1)) -blade_static_Sb*rotor_OMEGA^2*sin(azimuth_psi(t,2)) -blade_static_Sb*rotor_OMEGA^2*sin(azimuth_psi(t,3)) -blade_static_Sb*rotor_OMEGA^2*sin(azimuth_psi(t,4)) hub_stiffness_Kx 0;
    blade_static_Sb*rotor_OMEGA^2*cos(azimuth_psi(t,1)) blade_static_Sb*rotor_OMEGA^2*cos(azimuth_psi(t,2)) blade_static_Sb*rotor_OMEGA^2*cos(azimuth_psi(t,3)) blade_static_Sb*rotor_OMEGA^2*cos(azimuth_psi(t,4)) 0 hub_stiffness_Ky];

% additional stiffness

K_zetad=hub_to_A_a^2*stiffness_kd*sin(zeta_E-phi+gamma_E)+hub_to_B_b^2*stiffness_kd*sin(zeta_E+phi+gamma_E)-...
        -hub_to_A_a*stiffness_kd*(length_eq_lE-length_undef_l0)/length_eq_lE*cos(zeta_E-phi+gamma_E)*...
        (length_eq_lE-hub_to_A_a*cos(zeta_E-phi+gamma_E))-hub_to_B_b*stiffness_kd*...
        (length_eq_lE-length_undef_l0)/length_eq_lE*cos(zeta_E+phi+gamma_E)*...
        (length_eq_lE-hub_to_B_b*cos(zeta_E+phi+gamma_E));

K_zetaed = hub_to_A_a*hub_to_B_b*stiffness_kd*sin(zeta_E-phi+gamma_E)*sin(zeta_E+phi+gamma_E)+...
    hub_to_A_a*hub_to_B_b*stiffness_kd*(length_eq_lE-length_undef_l0)/length_undef_l0*...
    cos(zeta_E-phi+gamma_E)*cos(zeta_E+phi+gamma_E);

syms K_zetad K_zetaed
K_R_ADD = [K_zetad K_zetaed 0 K_zetaed 0 0;
           K_zetaed K_zetad K_zetaed 0 0 0;
           0 K_zetaed K_zetad K_zetaed 0 0;
           K_zetaed 0 K_zetaed K_zetad 0 0;
           0 0 0 0 0 0;
           0 0 0 0 0 0];

K_R = K_R_ADD + K_R;









T_F = [T zeros(4,2);
    zeros(2,4) [1 0; 0 1]];

T_F_dot = simplify(diff(T_F,t));

T_F_dotdot = simplify(diff(T_F_dot,t));

T_F_T = T_F';

M_NR = simplify(T_F_T* M_R * T_F);

C_NR = simplify(T_F_T*(2*M_R*T_F_dot+C_R*T_F));

K_NR = simplify(T_F_T*(M_R*T_F_dotdot+C_R*T_F_dot+K_R*T_F));


