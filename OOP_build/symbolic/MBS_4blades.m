clear all 
close all
clc 



syms  i

rotor_OMEGA = sym('rotor_OMEGA','real');

t = sym('t','real');

psi (t,i) = rotor_OMEGA*t+2*sym(pi)/4*(i-1);
azimuth_psi (t,i) = rotor_OMEGA*t+2*sym(pi)/4*(i-1);

T(t) = [1 cos(psi(t,1)) sin(psi(t,1)) -1;
    1 cos(psi(t,2)) sin(psi(t,2)) 1;
    1 cos(psi(t,3)) sin(psi(t,3)) -1;
    1 cos(psi(t,4)) sin(psi(t,4)) 1];


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

syms blade_damping_C hub_damping_Cx hub_damping_Cy blade_damping_C1

C_R = [0 0 0 0 0 0;
    0 blade_damping_C 0 0 0 0;
    0 0 blade_damping_C 0 0 0;
    0 0 0 blade_damping_C 0 0;
    blade_static_Sb*2*rotor_OMEGA*cos(azimuth_psi(t,1)) blade_static_Sb*2*rotor_OMEGA*cos(azimuth_psi(t,2)) blade_static_Sb*2*rotor_OMEGA*cos(azimuth_psi(t,3)) blade_static_Sb*2*rotor_OMEGA*cos(azimuth_psi(t,4)) hub_damping_Cx 0;
    blade_static_Sb*2*rotor_OMEGA*sin(azimuth_psi(t,1)) blade_static_Sb*2*rotor_OMEGA*sin(azimuth_psi(t,2)) blade_static_Sb*2*rotor_OMEGA*sin(azimuth_psi(t,3)) blade_static_Sb*2*rotor_OMEGA*sin(azimuth_psi(t,4)) 0 hub_damping_Cy];

syms hub_stiffness_Kx hub_stiffness_Ky lag_hinge_offset_e blade_stiffness_K

K_R = [(lag_hinge_offset_e*rotor_OMEGA^2*blade_static_Sb+blade_stiffness_K) 0 0 0 0 0;
    0 (lag_hinge_offset_e*rotor_OMEGA^2*blade_static_Sb +blade_stiffness_K) 0 0 0 0;
    0 0 (lag_hinge_offset_e*rotor_OMEGA^2*blade_static_Sb+blade_stiffness_K) 0 0 0;
    0 0 0 (lag_hinge_offset_e*rotor_OMEGA^2*blade_static_Sb+blade_stiffness_K) 0 0;
    -blade_static_Sb*rotor_OMEGA^2*sin(azimuth_psi(t,1)) -blade_static_Sb*rotor_OMEGA^2*sin(azimuth_psi(t,2)) -blade_static_Sb*rotor_OMEGA^2*sin(azimuth_psi(t,3)) -blade_static_Sb*rotor_OMEGA^2*sin(azimuth_psi(t,4)) hub_stiffness_Kx 0;
    blade_static_Sb*rotor_OMEGA^2*cos(azimuth_psi(t,1)) blade_static_Sb*rotor_OMEGA^2*cos(azimuth_psi(t,2)) blade_static_Sb*rotor_OMEGA^2*cos(azimuth_psi(t,3)) blade_static_Sb*rotor_OMEGA^2*cos(azimuth_psi(t,4)) 0 hub_stiffness_Ky];
T_F = [T zeros(4,2);
    zeros(2,4) [1 0; 0 1]];

T_F_dot = simplify(diff(T_F,t));

T_F_dotdot = simplify(diff(T_F_dot,t));

T_F_T = T_F';

M_NR = simplify(T_F_T* M_R * T_F);

C_NR = simplify(T_F_T*(2*M_R*T_F_dot+C_R*T_F));

K_NR = simplify(T_F_T*(M_R*T_F_dotdot+C_R*T_F_dot+K_R*T_F));




M_NR_ODI=M_NR

C_NR_ODI = simplify(subs(C_NR,[blade_damping_C1 ],...
    [0]))

K_NR_ODI = simplify(subs(K_NR,[blade_damping_C1 ],...
    [0]))


A_NR_ODI = [-M_NR_ODI\eye(6)*C_NR_ODI  -M_NR_ODI\eye(6)*K_NR_ODI;
   eye(6) zeros(6)];

A_NR_ODI = simplify(A_NR_ODI)
