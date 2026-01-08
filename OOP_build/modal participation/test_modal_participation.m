clear all
close all
clc

%test for modal coefficients

OMEGA= 5;
T = 2*pi/OMEGA;
t= (0:0.0001:T);
% wave_1 = sin(OMEGA.*t) + 0.5*sin(2.*OMEGA.*t) + 0.25* sin(3*OMEGA.*t);

wave_1 = 5 + sin(OMEGA.*t) + 0.5 * cos(3*OMEGA.*t)- 4 * exp(1i*5*OMEGA.*t);
wave_1 = cos(OMEGA.*t) + wave_1;
wave_1 = wave_1 + 4.*wave_1 + 5*sin (OMEGA.*t);
wave_1 = wave_1 + wave_1/2;

plot(t,wave_1,'x-')
hold on
grid on

NH = 6;
complex_coeffs = complex_coeff (T,wave_1,NH,t);
nh=(-NH:NH);

for i=1:length(complex_coeffs)
    app_vect(i,:) = complex_coeffs(i)*exp(1i*nh(i)*OMEGA*t);
end

wave_1_app=sum(app_vect,1);

plot(t,wave_1_app,'-','color','r')


