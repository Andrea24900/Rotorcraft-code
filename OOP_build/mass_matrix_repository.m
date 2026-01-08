function mass_matrix_M = mass_matrix_repository(obj)
switch obj.problem.number_blades
    case 3
        
        mass_matrix_M = [3* obj.problem.rotor_characteristics.blade_inertia_Ib,        0,         0,             0,            0;
            0,        3/2* obj.problem.rotor_characteristics.blade_inertia_Ib,         0, 0,      -3/2* obj.problem.rotor_characteristics.blade_static_Sb;
            0,        0,         3/2* obj.problem.rotor_characteristics.blade_inertia_Ib,  3/2* obj.problem.rotor_characteristics.blade_static_Sb,            0;
            0,        0, 3/2* obj.problem.rotor_characteristics.blade_static_Sb,               obj.problem.rotor_characteristics.hub_mass_Mx +  obj.problem.rotor_characteristics.blade_mass_Mb*obj.problem.number_blades,            0;
            0, -3/2* obj.problem.rotor_characteristics.blade_static_Sb,         0,              0,             obj.problem.rotor_characteristics.hub_mass_My +  obj.problem.rotor_characteristics.blade_mass_Mb*obj.problem.number_blades];
    case 4
        
        mass_matrix_M = [4* obj.problem.rotor_characteristics.blade_inertia_Ib,        0,         0, 0,             0,            0;
            0,        2* obj.problem.rotor_characteristics.blade_inertia_Ib,         0, 0,             0, -2* obj.problem.rotor_characteristics.blade_static_Sb;
            0,        0,         2* obj.problem.rotor_characteristics.blade_inertia_Ib, 0, 2* obj.problem.rotor_characteristics.blade_static_Sb,            0;
            0,        0,         0, 4* obj.problem.rotor_characteristics.blade_inertia_Ib,             0,            0;
            0,        0, 2* obj.problem.rotor_characteristics.blade_static_Sb, 0,              obj.problem.rotor_characteristics.hub_mass_Mx +  obj.problem.rotor_characteristics.blade_mass_Mb*obj.problem.number_blades,            0;
            0, -2* obj.problem.rotor_characteristics.blade_static_Sb,         0, 0,             0,             obj.problem.rotor_characteristics.hub_mass_My +  obj.problem.rotor_characteristics.blade_mass_Mb*obj.problem.number_blades];


    case 5
        mass_matrix_M = [5* obj.problem.rotor_characteristics.blade_inertia_Ib,        0,        0,        0,        0,         0,         0;
            0, (5* obj.problem.rotor_characteristics.blade_inertia_Ib)/2, 0,        0,        0,         0, -(5* obj.problem.rotor_characteristics.blade_static_Sb)/2;
            0,        0, (5* obj.problem.rotor_characteristics.blade_inertia_Ib)/2, 0,        0, (5* obj.problem.rotor_characteristics.blade_static_Sb)/2, 0;
            0,        0,        0, (5* obj.problem.rotor_characteristics.blade_inertia_Ib)/2, 0,         0,         0;
            0,        0,        0,        0, (5* obj.problem.rotor_characteristics.blade_inertia_Ib)/2, 0,         0;
            0,        0, (5* obj.problem.rotor_characteristics.blade_static_Sb)/2, 0,        0, obj.problem.rotor_characteristics.hub_mass_Mx + obj.problem.rotor_characteristics.blade_mass_Mb*obj.problem.number_blades, 0;
            0, -(5* obj.problem.rotor_characteristics.blade_static_Sb)/2, 0,        0,        0,         0, obj.problem.rotor_characteristics.hub_mass_My + obj.problem.rotor_characteristics.blade_mass_Mb*obj.problem.number_blades];

    case 7
        mass_matrix_M = [
            7*obj.problem.rotor_characteristics.blade_inertia_Ib, 0, 0, 0, 0, 0, 0, 0, 0;
            0, (7*obj.problem.rotor_characteristics.blade_inertia_Ib)/2, (obj.problem.rotor_characteristics.blade_inertia_Ib*sin(2*rotor_OMEGA*t)*(4*cos(pi/7) + 4*cos(pi/7)^2 - 8*cos(pi/7)^3 - 1))/2, 0, 0, 0, 0, (obj.problem.rotor_characteristics.blade_static_Sb*sin(2*rotor_OMEGA*t)*(4*cos(pi/7) + 4*cos(pi/7)^2 - 8*cos(pi/7)^3 - 1))/2, -(7*obj.problem.rotor_characteristics.blade_static_Sb)/2;
            0, (obj.problem.rotor_characteristics.blade_inertia_Ib*sin(2*rotor_OMEGA*t)*(4*cos(pi/7) + 4*cos(pi/7)^2 - 8*cos(pi/7)^3 - 1))/2, (7*obj.problem.rotor_characteristics.blade_inertia_Ib)/2, 0, 0, 0, 0, (7*obj.problem.rotor_characteristics.blade_static_Sb)/2, -(obj.problem.rotor_characteristics.blade_static_Sb*sin(2*rotor_OMEGA*t)*(4*cos(pi/7) + 4*cos(pi/7)^2 - 8*cos(pi/7)^3 - 1))/2;
            0, 0, 0, (7*obj.problem.rotor_characteristics.blade_inertia_Ib)/2, (obj.problem.rotor_characteristics.blade_inertia_Ib*sin(4*rotor_OMEGA*t)*(4*cos(pi/7) + 4*cos(pi/7)^2 - 8*cos(pi/7)^3 - 1))/2, 0, 0, 0, 0;
            0, 0, 0, (obj.problem.rotor_characteristics.blade_inertia_Ib*sin(4*rotor_OMEGA*t)*(4*cos(pi/7) + 4*cos(pi/7)^2 - 8*cos(pi/7)^3 - 1))/2, (7*obj.problem.rotor_characteristics.blade_inertia_Ib)/2, 0, 0, 0, 0;
            0, 0, 0, 0, 0, (7*obj.problem.rotor_characteristics.blade_inertia_Ib)/2, (obj.problem.rotor_characteristics.blade_inertia_Ib*sin(6*rotor_OMEGA*t)*(4*cos(pi/7) + 4*cos(pi/7)^2 - 8*cos(pi/7)^3 - 1))/2, 0, 0;
            0, 0, 0, 0, 0, (obj.problem.rotor_characteristics.blade_inertia_Ib*sin(6*rotor_OMEGA*t)*(4*cos(pi/7) + 4*cos(pi/7)^2 - 8*cos(pi/7)^3 - 1))/2, (7*obj.problem.rotor_characteristics.blade_inertia_Ib)/2, 0, 0;
            0, (obj.problem.rotor_characteristics.blade_static_Sb*sin(2*rotor_OMEGA*t)*(4*cos(pi/7) + 4*cos(pi/7)^2 - 8*cos(pi/7)^3 - 1))/2, (7*obj.problem.rotor_characteristics.blade_static_Sb)/2, 0, 0, 0, 0, 7*obj.problem.rotor_characteristics.blade_mass_Mb + obj.problem.rotor_characteristics.hub_mass_Mx, 0;
            0, -(7*obj.problem.rotor_characteristics.blade_static_Sb)/2, -(obj.problem.rotor_characteristics.blade_static_Sb*sin(2*rotor_OMEGA*t)*(4*cos(pi/7) + 4*cos(pi/7)^2 - 8*cos(pi/7)^3 - 1))/2, 0, 0, 0, 0, 0, 7*obj.problem.rotor_characteristics.blade_mass_Mb + obj.problem.rotor_characteristics.hub_mass_My
            ];

end
end



