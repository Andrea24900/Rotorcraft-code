function state_matrix_A=function_hammond_matrices (case_name)

run("hammond_parameters_US.m")
    run("extra_parameters_US.m")

switch case_name
    case 'interblade_NR_matrices_iso'

        mass_matrix_M = [4*blade_inertia_Ib,        0,         0, 0,             0,            0;
            0,        2*blade_inertia_Ib,         0, 0,             0, -2*blade_static_Sb;
            0,        0,         2*blade_inertia_Ib, 0, 2*blade_static_Sb,            0;
            0,        0,         0, 4*blade_inertia_Ib,             0,            0;
            0,        0, 2*blade_static_Sb, 0,             hub_mass_Mx + blade_mass_Mb*number_blades_Nb,            0;
            0, -2*blade_static_Sb,         0, 0,             0,            hub_mass_My + blade_mass_Mb*number_blades_Nb];


        %%
        C_zetad = hub_to_A_a^2*damping_cd*sin(zeta_E-phi+gamma_E)^2+hub_to_B_b^2*damping_cd*sin(zeta_E+phi+gamma_E)^2;

        C_zetaed = hub_to_A_a*hub_to_B_b*damping_cd*sin(zeta_E-phi+gamma_E)*sin(zeta_E+phi+gamma_E);

        damping_matrix_C=@(rotor_OMEGA) [4*C_zetad + 8*C_zetaed,                               0,                              0,                      0,              0,              0;
            0,                       2*C_zetad, 4*blade_inertia_Ib*rotor_OMEGA,                      0,              0,              0;
            0, -4*blade_inertia_Ib*rotor_OMEGA,                      2*C_zetad,                      0,              0,              0;
            0,                               0,                              0, 4*C_zetad - 8*C_zetaed,              0,              0;
            0,                               0,                              0,                      0, hub_damping_Cx,              0;
            0,                               0,                              0,                      0,              0, hub_damping_Cy];


        % stiffness

        K_zetad=hub_to_A_a^2*stiffness_kd*sin(zeta_E-phi+gamma_E)+hub_to_B_b^2*stiffness_kd*sin(zeta_E+phi+gamma_E)-...
            -hub_to_A_a*stiffness_kd*(length_eq_lE-length_undef_l0)/length_eq_lE*cos(zeta_E-phi+gamma_E)*...
            (length_eq_lE-hub_to_A_a*cos(zeta_E-phi+gamma_E))-hub_to_B_b*stiffness_kd*...
            (length_eq_lE-length_undef_l0)/length_eq_lE*cos(zeta_E+phi+gamma_E)*...
            (length_eq_lE-hub_to_B_b*cos(zeta_E+phi+gamma_E));

        K_zetaed = hub_to_A_a*hub_to_B_b*stiffness_kd*sin(zeta_E-phi+gamma_E)*sin(zeta_E+phi+gamma_E)+...
            hub_to_A_a*hub_to_B_b*stiffness_kd*(length_eq_lE-length_undef_l0)/length_undef_l0*...
            cos(zeta_E-phi+gamma_E)*cos(zeta_E+phi+gamma_E);

        stiffness_matrix_K =@(rotor_OMEGA) ...
            [4*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + 4*K_zetad + 8*K_zetaed,                                                                                                 0,                                                                                                 0,                                                                           0,                0,                0;
            0, 2*K_zetad - 2*blade_inertia_Ib*rotor_OMEGA^2 + 2*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2,                                                                             2*C_zetad*rotor_OMEGA,                                                                           0,                0,                0;
            0,                                                                            -2*C_zetad*rotor_OMEGA, 2*K_zetad - 2*blade_inertia_Ib*rotor_OMEGA^2 + 2*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2,                                                                           0,                0,                0;
            0,                                                                                                 0,                                                                                                 0, 4*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + 4*K_zetad - 8*K_zetaed,                0,                0;
            0,                                                                                                 0,                                                                                                 0,                                                                           0, hub_stiffness_Kx,                0;
            0,                                                                                                 0,                                                                                                 0,                                                                           0,                0, hub_stiffness_Ky];

        state_matrix_A =@(rotor_OMEGA) [-mass_matrix_M\eye(width(mass_matrix_M))*damping_matrix_C(rotor_OMEGA) -mass_matrix_M\eye(width(mass_matrix_M))*stiffness_matrix_K(rotor_OMEGA);
            eye(width(mass_matrix_M)) zeros(width(mass_matrix_M))];


    case 'interblade_NR_matrices_ODI'

        % mass matrix
        mass_matrix_M = [4*blade_inertia_Ib,        0,         0, 0,             0,            0;
            0,        2*blade_inertia_Ib,         0, 0,             0, -2*blade_static_Sb;
            0,        0,         2*blade_inertia_Ib, 0, 2*blade_static_Sb,            0;
            0,        0,         0, 4*blade_inertia_Ib,             0,            0;
            0,        0, 2*blade_static_Sb, 0,             hub_mass_Mx + blade_mass_Mb*number_blades_Nb,            0;
            0, -2*blade_static_Sb,         0, 0,             0,            hub_mass_My + blade_mass_Mb*number_blades_Nb];



        %damping matrix
        C_zetad = hub_to_A_a^2*damping_cd*sin(zeta_E-phi+gamma_E)^2+hub_to_B_b^2*damping_cd*sin(zeta_E+phi+gamma_E)^2;

        C_zetaed = hub_to_A_a*hub_to_B_b*damping_cd*sin(zeta_E-phi+gamma_E)*sin(zeta_E+phi+gamma_E);

        C_11 = hub_to_B_b^2*damping_cd*sin(zeta_E+phi+gamma_E)^2;
        C_22 = hub_to_A_a^2*damping_cd*sin(zeta_E-phi+gamma_E)^2;


        damping_matrix_C=@(t,rotor_OMEGA)...
            [                                                                                                                                   C_11 + C_22 + 2*C_zetad + 6*C_zetaed, C_11*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) - C_zetaed*cos(rotor_OMEGA*t) - C_22*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) + C_zetaed*sin(rotor_OMEGA*t), C_22*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) - C_zetaed*cos(rotor_OMEGA*t) + C_11*sin(rotor_OMEGA*t) - C_zetad*sin(rotor_OMEGA*t) - C_zetaed*sin(rotor_OMEGA*t),                                                                                                                                                             C_22 - C_11,              0,              0;
            C_11*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) - C_zetaed*cos(rotor_OMEGA*t) - C_22*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) + C_zetaed*sin(rotor_OMEGA*t),                                               C_11/2 + C_22/2 + C_zetad + (C_11*cos(2*rotor_OMEGA*t))/2 - (C_22*cos(2*rotor_OMEGA*t))/2 + C_zetaed*sin(2*rotor_OMEGA*t),                                          4*blade_inertia_Ib*rotor_OMEGA - C_zetaed*cos(2*rotor_OMEGA*t) + (C_11*sin(2*rotor_OMEGA*t))/2 - (C_22*sin(2*rotor_OMEGA*t))/2, C_zetad*cos(rotor_OMEGA*t) - C_11*cos(rotor_OMEGA*t) - C_zetaed*cos(rotor_OMEGA*t) - C_22*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) - C_zetaed*sin(rotor_OMEGA*t),              0,              0;
            C_22*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) - C_zetaed*cos(rotor_OMEGA*t) + C_11*sin(rotor_OMEGA*t) - C_zetad*sin(rotor_OMEGA*t) - C_zetaed*sin(rotor_OMEGA*t),                                          (C_11*sin(2*rotor_OMEGA*t))/2 - C_zetaed*cos(2*rotor_OMEGA*t) - 4*blade_inertia_Ib*rotor_OMEGA - (C_22*sin(2*rotor_OMEGA*t))/2,                                               C_11/2 + C_22/2 + C_zetad - (C_11*cos(2*rotor_OMEGA*t))/2 + (C_22*cos(2*rotor_OMEGA*t))/2 - C_zetaed*sin(2*rotor_OMEGA*t), C_22*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) + C_zetaed*cos(rotor_OMEGA*t) - C_11*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) - C_zetaed*sin(rotor_OMEGA*t),              0,              0;
            C_22 - C_11, C_zetad*cos(rotor_OMEGA*t) - C_11*cos(rotor_OMEGA*t) - C_zetaed*cos(rotor_OMEGA*t) - C_22*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) - C_zetaed*sin(rotor_OMEGA*t), C_22*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) + C_zetaed*cos(rotor_OMEGA*t) - C_11*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) - C_zetaed*sin(rotor_OMEGA*t),                                                                                                                                    C_11 + C_22 + 2*C_zetad - 6*C_zetaed,              0,              0;
            0,                                                                                                                                                                       0,                                                                                                                                                                       0,                                                                                                                                                                       0, hub_damping_Cx,              0;
            0,                                                                                                                                                                       0,                                                                                                                                                                       0,                                                                                                                                                                       0,              0, hub_damping_Cy];

        %stiffness matrix
        K_zetad=hub_to_A_a^2*stiffness_kd*sin(zeta_E-phi+gamma_E)+hub_to_B_b^2*stiffness_kd*sin(zeta_E+phi+gamma_E)-...
            -hub_to_A_a*stiffness_kd*(length_eq_lE-length_undef_l0)/length_eq_lE*cos(zeta_E-phi+gamma_E)*...
            (length_eq_lE-hub_to_A_a*cos(zeta_E-phi+gamma_E))-hub_to_B_b*stiffness_kd*...
            (length_eq_lE-length_undef_l0)/length_eq_lE*cos(zeta_E+phi+gamma_E)*...
            (length_eq_lE-hub_to_B_b*cos(zeta_E+phi+gamma_E));

        K_zetaed = hub_to_A_a*hub_to_B_b*stiffness_kd*sin(zeta_E-phi+gamma_E)*sin(zeta_E+phi+gamma_E)+...
            hub_to_A_a*hub_to_B_b*stiffness_kd*(length_eq_lE-length_undef_l0)/length_undef_l0*...
            cos(zeta_E-phi+gamma_E)*cos(zeta_E+phi+gamma_E);
        stiffness_matrix_K =@(t,rotor_OMEGA) ...
            [4*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + 4*K_zetad + 8*K_zetaed,                                                 rotor_OMEGA*(C_zetad*cos(rotor_OMEGA*t) - C_22*cos(rotor_OMEGA*t) + C_zetaed*cos(rotor_OMEGA*t) - C_11*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) + C_zetaed*sin(rotor_OMEGA*t)),                                                 rotor_OMEGA*(C_11*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) - C_zetaed*cos(rotor_OMEGA*t) - C_22*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) + C_zetaed*sin(rotor_OMEGA*t)),                                                                           0,                0,                0;
            0, 2*K_zetad - 2*blade_inertia_Ib*rotor_OMEGA^2 + 2*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + C_zetaed*rotor_OMEGA*cos(2*rotor_OMEGA*t) - (C_11*rotor_OMEGA*sin(2*rotor_OMEGA*t))/2 + (C_22*rotor_OMEGA*sin(2*rotor_OMEGA*t))/2,                                                                                                   (rotor_OMEGA*(C_11 + C_22 + 2*C_zetad + C_11*cos(2*rotor_OMEGA*t) - C_22*cos(2*rotor_OMEGA*t) + 2*C_zetaed*sin(2*rotor_OMEGA*t)))/2,                                                                           0,                0,                0;
            0,                                                                                                  -(rotor_OMEGA*(C_11 + C_22 + 2*C_zetad - C_11*cos(2*rotor_OMEGA*t) + C_22*cos(2*rotor_OMEGA*t) - 2*C_zetaed*sin(2*rotor_OMEGA*t)))/2, 2*K_zetad - 2*blade_inertia_Ib*rotor_OMEGA^2 + 2*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 - C_zetaed*rotor_OMEGA*cos(2*rotor_OMEGA*t) + (C_11*rotor_OMEGA*sin(2*rotor_OMEGA*t))/2 - (C_22*rotor_OMEGA*sin(2*rotor_OMEGA*t))/2,                                                                           0,                0,                0;
            0,                                                -rotor_OMEGA*(C_22*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) + C_zetaed*cos(rotor_OMEGA*t) - C_11*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) - C_zetaed*sin(rotor_OMEGA*t)),                                                -rotor_OMEGA*(C_11*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) + C_zetaed*cos(rotor_OMEGA*t) + C_22*sin(rotor_OMEGA*t) - C_zetad*sin(rotor_OMEGA*t) + C_zetaed*sin(rotor_OMEGA*t)), 4*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + 4*K_zetad - 8*K_zetaed,                0,                0;
            0,                                                                                                                                                                                                                                     0,                                                                                                                                                                                                                                     0,                                                                           0, hub_stiffness_Kx,                0;
            0,                                                                                                                                                                                                                                     0,                                                                                                                                                                                                                                     0,                                                                           0,                0, hub_stiffness_Ky];


        %state matrix
        state_matrix_A =@(t,rotor_OMEGA) [-mass_matrix_M\eye(6)*damping_matrix_C(t,rotor_OMEGA) -mass_matrix_M\eye(6)*stiffness_matrix_K(t,rotor_OMEGA);
            eye(6) zeros(6)];

    case 'interblade_NR_matrices_2DI_adjacent'

        % mass matrix
        mass_matrix_M = [4*blade_inertia_Ib,        0,         0, 0,             0,            0;
            0,        2*blade_inertia_Ib,         0, 0,             0, -2*blade_static_Sb;
            0,        0,         2*blade_inertia_Ib, 0, 2*blade_static_Sb,            0;
            0,        0,         0, 4*blade_inertia_Ib,             0,            0;
            0,        0, 2*blade_static_Sb, 0,             hub_mass_Mx + blade_mass_Mb*number_blades_Nb,            0;
            0, -2*blade_static_Sb,         0, 0,             0,            hub_mass_My + blade_mass_Mb*number_blades_Nb];



        %damping matrix
        C_zetad = hub_to_A_a^2*damping_cd*sin(zeta_E-phi+gamma_E)^2+hub_to_B_b^2*damping_cd*sin(zeta_E+phi+gamma_E)^2;

        C_zetaed = hub_to_A_a*hub_to_B_b*damping_cd*sin(zeta_E-phi+gamma_E)*sin(zeta_E+phi+gamma_E);

        C_11 = hub_to_B_b^2*damping_cd*sin(zeta_E+phi+gamma_E)^2;
        C_22 = hub_to_A_a^2*damping_cd*sin(zeta_E-phi+gamma_E)^2;


        damping_matrix_C=@(t,rotor_OMEGA)...
[                                                                                    C_11 + 2*C_zetad + 4*C_zetaed, C_11*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) + 2*C_zetaed*sin(rotor_OMEGA*t), C_11*sin(rotor_OMEGA*t) - 2*C_zetaed*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) - C_zetad*sin(rotor_OMEGA*t),                                                                                                             -C_11,              0,              0;
C_11*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) + 2*C_zetaed*sin(rotor_OMEGA*t),                                                                      - C_11*sin(rotor_OMEGA*t)^2 + C_11 + C_zetad,                                                    4*blade_inertia_Ib*rotor_OMEGA + (C_11*sin(2*rotor_OMEGA*t))/2, C_zetad*cos(rotor_OMEGA*t) - C_11*cos(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) - 2*C_zetaed*sin(rotor_OMEGA*t),              0,              0;
C_11*sin(rotor_OMEGA*t) - 2*C_zetaed*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) - C_zetad*sin(rotor_OMEGA*t),                                                    (C_11*sin(2*rotor_OMEGA*t))/2 - 4*blade_inertia_Ib*rotor_OMEGA,                                                                               C_11*sin(rotor_OMEGA*t)^2 + C_zetad, 2*C_zetaed*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) - C_11*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t),              0,              0;
                                                                                                            -C_11, C_zetad*cos(rotor_OMEGA*t) - C_11*cos(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) - 2*C_zetaed*sin(rotor_OMEGA*t), 2*C_zetaed*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) - C_11*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t),                                                                                     C_11 + 2*C_zetad - 4*C_zetaed,              0,              0;
                                                                                                                0,                                                                                                                 0,                                                                                                                 0,                                                                                                                 0, hub_damping_Cx,              0;
                                                                                                                0,                                                                                                                 0,                                                                                                                 0,                                                                                                                 0,              0, hub_damping_Cy];

        %stiffness matrix
        K_zetad=hub_to_A_a^2*stiffness_kd*sin(zeta_E-phi+gamma_E)+hub_to_B_b^2*stiffness_kd*sin(zeta_E+phi+gamma_E)-...
            -hub_to_A_a*stiffness_kd*(length_eq_lE-length_undef_l0)/length_eq_lE*cos(zeta_E-phi+gamma_E)*...
            (length_eq_lE-hub_to_A_a*cos(zeta_E-phi+gamma_E))-hub_to_B_b*stiffness_kd*...
            (length_eq_lE-length_undef_l0)/length_eq_lE*cos(zeta_E+phi+gamma_E)*...
            (length_eq_lE-hub_to_B_b*cos(zeta_E+phi+gamma_E));

        K_zetaed = hub_to_A_a*hub_to_B_b*stiffness_kd*sin(zeta_E-phi+gamma_E)*sin(zeta_E+phi+gamma_E)+...
            hub_to_A_a*hub_to_B_b*stiffness_kd*(length_eq_lE-length_undef_l0)/length_undef_l0*...
            cos(zeta_E-phi+gamma_E)*cos(zeta_E+phi+gamma_E);
        stiffness_matrix_K =@(t,rotor_OMEGA) ...
[4*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + 4*K_zetad + 8*K_zetaed,               rotor_OMEGA*(C_zetad*cos(rotor_OMEGA*t) + 2*C_zetaed*cos(rotor_OMEGA*t) - C_11*sin(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t)),               rotor_OMEGA*(C_11*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) + C_zetad*sin(rotor_OMEGA*t) + 2*C_zetaed*sin(rotor_OMEGA*t)),                                                                           0,                0,                0;
                                                                          0, 2*K_zetad - 2*blade_inertia_Ib*rotor_OMEGA^2 + 2*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 - (C_11*rotor_OMEGA*sin(2*rotor_OMEGA*t))/2,                                                                                (rotor_OMEGA*(C_11 + 2*C_zetad + C_11*cos(2*rotor_OMEGA*t)))/2,                                                                           0,                0,                0;
                                                                          0,                                                                                            -rotor_OMEGA*(C_11*sin(rotor_OMEGA*t)^2 + C_zetad), 2*K_zetad - 2*blade_inertia_Ib*rotor_OMEGA^2 + 2*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + (C_11*rotor_OMEGA*sin(2*rotor_OMEGA*t))/2,                                                                           0,                0,                0;
                                                                          0,               rotor_OMEGA*(C_zetad*cos(rotor_OMEGA*t) - 2*C_zetaed*cos(rotor_OMEGA*t) + C_11*sin(rotor_OMEGA*t) - C_zetad*sin(rotor_OMEGA*t)),              -rotor_OMEGA*(C_11*cos(rotor_OMEGA*t) - C_zetad*cos(rotor_OMEGA*t) - C_zetad*sin(rotor_OMEGA*t) + 2*C_zetaed*sin(rotor_OMEGA*t)), 4*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + 4*K_zetad - 8*K_zetaed,                0,                0;
                                                                          0,                                                                                                                                             0,                                                                                                                                             0,                                                                           0, hub_stiffness_Kx,                0;
                                                                          0,                                                                                                                                             0,                                                                                                                                             0,                                                                           0,                0, hub_stiffness_Ky];


        %state matrix
        state_matrix_A =@(t,rotor_OMEGA) [-mass_matrix_M\eye(6)*damping_matrix_C(t,rotor_OMEGA) -mass_matrix_M\eye(6)*stiffness_matrix_K(t,rotor_OMEGA);
            eye(6) zeros(6)];
    case 'interblade_NR_matrices_2DI_opposite'

        % mass matrix
        mass_matrix_M = [4*blade_inertia_Ib,        0,         0, 0,             0,            0;
            0,        2*blade_inertia_Ib,         0, 0,             0, -2*blade_static_Sb;
            0,        0,         2*blade_inertia_Ib, 0, 2*blade_static_Sb,            0;
            0,        0,         0, 4*blade_inertia_Ib,             0,            0;
            0,        0, 2*blade_static_Sb, 0,             hub_mass_Mx + blade_mass_Mb*number_blades_Nb,            0;
            0, -2*blade_static_Sb,         0, 0,             0,            hub_mass_My + blade_mass_Mb*number_blades_Nb];



        %damping matrix
        C_zetad = hub_to_A_a^2*damping_cd*sin(zeta_E-phi+gamma_E)^2+hub_to_B_b^2*damping_cd*sin(zeta_E+phi+gamma_E)^2;

        C_zetaed = hub_to_A_a*hub_to_B_b*damping_cd*sin(zeta_E-phi+gamma_E)*sin(zeta_E+phi+gamma_E);

        C_11 = hub_to_B_b^2*damping_cd*sin(zeta_E+phi+gamma_E)^2;
        C_22 = hub_to_A_a^2*damping_cd*sin(zeta_E-phi+gamma_E)^2;


        damping_matrix_C=@(t,rotor_OMEGA)...
[2*C_11 + 2*C_22 + 4*C_zetaed,                                                                                                                        0,                                                                                                                        0,              2*C_22 - 2*C_11,              0,              0;
                           0,             2*C_11*cos(rotor_OMEGA*t)^2 + 2*C_22*sin(rotor_OMEGA*t)^2 + 4*C_zetaed*cos(rotor_OMEGA*t)*sin(rotor_OMEGA*t), 4*blade_inertia_Ib*rotor_OMEGA - 2*C_zetaed*cos(2*rotor_OMEGA*t) + C_11*sin(2*rotor_OMEGA*t) - C_22*sin(2*rotor_OMEGA*t),                            0,              0,              0;
                           0, C_11*sin(2*rotor_OMEGA*t) - 2*C_zetaed*cos(2*rotor_OMEGA*t) - 4*blade_inertia_Ib*rotor_OMEGA - C_22*sin(2*rotor_OMEGA*t),             2*C_22*cos(rotor_OMEGA*t)^2 + 2*C_11*sin(rotor_OMEGA*t)^2 - 4*C_zetaed*cos(rotor_OMEGA*t)*sin(rotor_OMEGA*t),                            0,              0,              0;
             2*C_22 - 2*C_11,                                                                                                                        0,                                                                                                                        0, 2*C_11 + 2*C_22 - 4*C_zetaed,              0,              0;
                           0,                                                                                                                        0,                                                                                                                        0,                            0, hub_damping_Cx,              0;
                           0,                                                                                                                        0,                                                                                                                        0,                            0,              0, hub_damping_Cy];
 
        %stiffness matrix
        K_zetad=hub_to_A_a^2*stiffness_kd*sin(zeta_E-phi+gamma_E)+hub_to_B_b^2*stiffness_kd*sin(zeta_E+phi+gamma_E)-...
            -hub_to_A_a*stiffness_kd*(length_eq_lE-length_undef_l0)/length_eq_lE*cos(zeta_E-phi+gamma_E)*...
            (length_eq_lE-hub_to_A_a*cos(zeta_E-phi+gamma_E))-hub_to_B_b*stiffness_kd*...
            (length_eq_lE-length_undef_l0)/length_eq_lE*cos(zeta_E+phi+gamma_E)*...
            (length_eq_lE-hub_to_B_b*cos(zeta_E+phi+gamma_E));

        K_zetaed = hub_to_A_a*hub_to_B_b*stiffness_kd*sin(zeta_E-phi+gamma_E)*sin(zeta_E+phi+gamma_E)+...
            hub_to_A_a*hub_to_B_b*stiffness_kd*(length_eq_lE-length_undef_l0)/length_undef_l0*...
            cos(zeta_E-phi+gamma_E)*cos(zeta_E+phi+gamma_E);
        stiffness_matrix_K =@(t,rotor_OMEGA) ...
[4*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + 4*K_zetad + 8*K_zetaed,                                                                                                                                                                                                                               0,                                                                                                                                                                                                                               0,                                                                           0,                0,                0;
                                                                          0, 2*K_zetad - 2*blade_inertia_Ib*rotor_OMEGA^2 + 2*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + 2*C_zetaed*rotor_OMEGA*cos(2*rotor_OMEGA*t) - C_11*rotor_OMEGA*sin(2*rotor_OMEGA*t) + C_22*rotor_OMEGA*sin(2*rotor_OMEGA*t),                                                                                2*C_11*rotor_OMEGA*cos(rotor_OMEGA*t)^2 + 2*C_22*rotor_OMEGA*sin(rotor_OMEGA*t)^2 + 4*C_zetaed*rotor_OMEGA*cos(rotor_OMEGA*t)*sin(rotor_OMEGA*t),                                                                           0,                0,                0;
                                                                          0,                                                                                4*C_zetaed*rotor_OMEGA*cos(rotor_OMEGA*t)*sin(rotor_OMEGA*t) - 2*C_11*rotor_OMEGA*sin(rotor_OMEGA*t)^2 - 2*C_22*rotor_OMEGA*cos(rotor_OMEGA*t)^2, 2*K_zetad - 2*blade_inertia_Ib*rotor_OMEGA^2 + 2*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 - 2*C_zetaed*rotor_OMEGA*cos(2*rotor_OMEGA*t) + C_11*rotor_OMEGA*sin(2*rotor_OMEGA*t) - C_22*rotor_OMEGA*sin(2*rotor_OMEGA*t),                                                                           0,                0,                0;
                                                                         0,                                                                                                                                                                                                                               0,                                                                                                                                                                                                                               0, 4*blade_static_Sb*lag_hinge_offset_e*rotor_OMEGA^2 + 4*K_zetad - 8*K_zetaed,                0,                0;
                                                                          0,                                                                                                                                                                                                                               0,                                                                                                                                                                                                                               0,                                                                           0, hub_stiffness_Kx,                0;
                                                                          0,                                                                                                                                                                                                                               0,                                                                                                                                                                                                                               0,                                                                           0,                0, hub_stiffness_Ky];


        %state matrix
        state_matrix_A =@(t,rotor_OMEGA) [-mass_matrix_M\eye(6)*damping_matrix_C(t,rotor_OMEGA) -mass_matrix_M\eye(6)*stiffness_matrix_K(t,rotor_OMEGA);
            eye(6) zeros(6)];




end

end