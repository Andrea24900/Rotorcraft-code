clear all 
close all
clc


syms  i

OMEGA = sym('OMEGA','real');

t = sym('t','real');

psi (t,i) = OMEGA*t+2*sym(pi)/7*(i-1);


T(t) = [1 cos(psi(t,1)) sin(psi(t,1)) cos(2*psi(t,1)) sin(2*psi(t,1)) cos(3*psi(t,1)) sin(3*psi(t,1));
    1 cos(psi(t,2)) sin(psi(t,2)) cos(2*psi(t,2)) sin(2*psi(t,2)) cos(3*psi(t,2)) sin(3*psi(t,2));
    1 cos(psi(t,3)) sin(psi(t,3)) cos(2*psi(t,3)) sin(2*psi(t,3)) cos(3*psi(t,3)) sin(3*psi(t,3));
    1 cos(psi(t,4)) sin(psi(t,4))  cos(2*psi(t,4)) sin(2*psi(t,4)) cos(3*psi(t,4)) sin(3*psi(t,4));
    1 cos(psi(t,5)) sin(psi(t,5))  cos(2*psi(t,5)) sin(2*psi(t,5)) cos(3*psi(t,5)) sin(3*psi(t,5));
    1 cos(psi(t,6)) sin(psi(t,6))  cos(2*psi(t,6)) sin(2*psi(t,6)) cos(3*psi(t,6)) sin(3*psi(t,6));
    1 cos(psi(t,7)) sin(psi(t,7))  cos(2*psi(t,7)) sin(2*psi(t,7)) cos(3*psi(t,7)) sin(3*psi(t,7))];


T_T = simplify(T');



T_dot = diff(T,t);

T_dotdot = diff(T_dot,t);

syms Mb Mx My Sb Ib Nb

M_R = [
    Ib, 0, 0, 0, 0, 0, 0, Sb*sin(psi(t,1)), -Sb*cos(psi(t,1));
    0, Ib, 0, 0, 0, 0, 0, Sb*sin(psi(t,2)), -Sb*cos(psi(t,2));
    0, 0, Ib, 0, 0, 0, 0, Sb*sin(psi(t,3)), -Sb*cos(psi(t,3));
    0, 0, 0, Ib, 0, 0, 0, Sb*sin(psi(t,4)), -Sb*cos(psi(t,4));
    0, 0, 0, 0, Ib, 0, 0, Sb*sin(psi(t,5)), -Sb*cos(psi(t,5));
    0, 0, 0, 0, 0, Ib, 0, Sb*sin(psi(t,6)), -Sb*cos(psi(t,6));
    0, 0, 0, 0, 0, 0, Ib, Sb*sin(psi(t,7)), -Sb*cos(psi(t,7));
    Sb*sin(psi(t,1)), Sb*sin(psi(t,2)), Sb*sin(psi(t,3)), Sb*sin(psi(t,4)), Sb*sin(psi(t,5)), Sb*sin(psi(t,6)), Sb*sin(psi(t,7)), Mx+7*Mb, 0;
    -Sb*cos(psi(t,1)), -Sb*cos(psi(t,2)), -Sb*cos(psi(t,3)), -Sb*cos(psi(t,4)), -Sb*cos(psi(t,5)), -Sb*cos(psi(t,6)), -Sb*cos(psi(t,7)), 0, My+7*Mb
];

syms C Cx Cy C1

C_R = [
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, C, 0, 0, 0, 0, 0, 0, 0;
    0, 0, C, 0, 0, 0, 0, 0, 0;
    0, 0, 0, C, 0, 0, 0, 0, 0;
    0, 0, 0, 0, C, 0, 0, 0, 0;
    0, 0, 0, 0, 0, C, 0, 0, 0;
    0, 0, 0, 0, 0, 0, C, 0, 0;
    Sb*2*OMEGA*cos(psi(t,1)), Sb*2*OMEGA*cos(psi(t,2)), Sb*2*OMEGA*cos(psi(t,3)), Sb*2*OMEGA*cos(psi(t,4)), Sb*2*OMEGA*cos(psi(t,5)), Sb*2*OMEGA*cos(psi(t,6)), Sb*2*OMEGA*cos(psi(t,7)), Cx, 0;
    Sb*2*OMEGA*sin(psi(t,1)), Sb*2*OMEGA*sin(psi(t,2)), Sb*2*OMEGA*sin(psi(t,3)), Sb*2*OMEGA*sin(psi(t,4)), Sb*2*OMEGA*sin(psi(t,5)), Sb*2*OMEGA*sin(psi(t,6)), Sb*2*OMEGA*sin(psi(t,7)), 0, Cy
];


syms Kx Ky e K
K_R = [
    e*OMEGA^2*Sb + K, 0, 0, 0, 0, 0, 0, 0, 0;
    0, e*OMEGA^2*Sb + K, 0, 0, 0, 0, 0, 0, 0;
    0, 0, e*OMEGA^2*Sb + K, 0, 0, 0, 0, 0, 0;
    0, 0, 0, e*OMEGA^2*Sb + K, 0, 0, 0, 0, 0;
    0, 0, 0, 0, e*OMEGA^2*Sb + K, 0, 0, 0, 0;
    0, 0, 0, 0, 0, e*OMEGA^2*Sb + K, 0, 0, 0;
    0, 0, 0, 0, 0, 0, e*OMEGA^2*Sb + K, 0, 0;
    -Sb*OMEGA^2*sin(psi(t,1)), -Sb*OMEGA^2*sin(psi(t,2)), -Sb*OMEGA^2*sin(psi(t,3)), -Sb*OMEGA^2*sin(psi(t,4)), -Sb*OMEGA^2*sin(psi(t,5)), -Sb*OMEGA^2*sin(psi(t,6)), -Sb*OMEGA^2*sin(psi(t,7)), Kx, 0;
    Sb*OMEGA^2*cos(psi(t,1)), Sb*OMEGA^2*cos(psi(t,2)), Sb*OMEGA^2*cos(psi(t,3)), Sb*OMEGA^2*cos(psi(t,4)), Sb*OMEGA^2*cos(psi(t,5)), Sb*OMEGA^2*cos(psi(t,6)), Sb*OMEGA^2*cos(psi(t,7)), 0, Ky
];


T_F = [T zeros(7,2);
    zeros(2,7) [1 0; 0 1]];

T_F_dot = simplify(diff(T_F,t));

T_F_dotdot = simplify(diff(T_F_dot,t));

T_F_T = T_F';

M_NR = simplify(T_F_T* M_R * T_F)

C_NR = simplify(T_F_T*(2*M_R*T_F_dot+C_R*T_F))

K_NR = simplify(T_F_T*(M_R*T_F_dotdot+C_R*T_F_dot+K_R*T_F))




