clear all, close all, clc


%% INTRODUCTION OF NON-LINEAR FUNCTION TO SOLVE


OM=27;

rotor_build = rotor_build.build_all;

%% Generation of data (snapshots) and SVD
delta_T = 1e-2; %NO effect
time = (0:delta_T:7); %good fit (with higher values we get a worse result)
x_0 = 0.1*ones(12,1); % IF outside limit cycle: good prediction
                       % IF inside limit cycle: bad prediction caused by
                       % snapshots on divergent part of simulation
[t,x] = ode45(@(t,x) odefun(t,x,rotor_build,OM),time,x_0);
x=x';

X = x(:,1:end-1);
X2 = x(:,2:end);
% X = x(:,ceil(length(time)/2):end-1); %phase change with respect to time dynamics, but overall good fit
% X2 = x(:,ceil(length(time)/2)+1:end);
[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)
r = 12;  % right now best value
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,Lambda] = eig(Atilde);
Phi = X2*V*inv(S)*W;

%% Rebuild state (discrete)
alpha1=S*V(1,:)';
b =(W*Lambda)\alpha1;
x_rebuilt = Phi*Lambda*b;

for i=2:length(time)

    x_rebuilt(:,i) = Phi*Lambda^(i-1)*b;

end
%% Rebuild state (continuous)

%Omega = diag(1./delta_T.*log(diag(Lambda)));
Omega = 1/delta_T*log(Lambda);

x_reb_cont = @(t) Phi*exp(Omega.*t)*b;

for i=1:length(t)
    x_rebuilt_con(:,i)=x_reb_cont(t(i));
    err(:,i) = abs(x(:,i)-x_rebuilt_con(:,i));
end


plot_properties

%% PLOT

figure
plot(t,x(7,:),'linewidth',2,'Color',color.blue_mat,'DisplayName','LTI simulation')
hold on
grid on
plot(t,x(8:12,:),'linewidth',2,'Color',color.blue_mat,'HandleVisibility','off')
plot(t,x_rebuilt(7,:),'square','MarkerSize',12,'linewidth',2,'color',color.orange_mat,'DisplayName','DMD discrete reconstruction')
plot(t,x_rebuilt(8:12,:),'square','MarkerSize',12,'linewidth',2,'color',color.orange_mat,'HandleVisibility','off')
plot(t,x_rebuilt_con(7,:),'o','LineWidth',2,'color',color.yellow_mat,'DisplayName','DMD continuous reconstruction')
plot(t,x_rebuilt_con(8:12,:),'o','LineWidth',2,'color',color.yellow_mat,'HandleVisibility','off')
legend
xlabel('t')
ylabel('States')



%% STATE PREDICTION

time = (0:delta_T:20);

[t_2,x_2] = ode45(@(t,x) odefun(t,x,rotor_build,OM),time,x_0);
x_2=x_2';

for i=1:length(t_2)
    x_rebuilt_con_2(:,i)=x_reb_cont(t_2(i));
    err_2(:,i) = abs(x_2(:,i)-x_rebuilt_con_2(:,i));
end

x_rebuilt_con_2(:,1) = x_0;
diff_2 = x_2-x_rebuilt_con_2;
FIT_2 = 1 - norm(diff_2,'fro')/norm(x_2,'fro');

figure
plot(t_2,x_2(11,:),'o')
hold on
grid on
plot(t_2,x_rebuilt_con_2(11,:),'LineWidth',2)
xline(5,'LineWidth',2)
legend('LTI Simulation','DMD prediction','Start of DMD prediction')
text(1.5,-3.5, sprintf('FIT (all states) = %d percent', FIT_2),'FontSize',15)
xlabel('t')
ylabel('x_H')

figure
plot(x_2(11,:),x_2(5,:),'o')
hold on
grid on
plot(x_rebuilt_con_2(11,:),x_rebuilt_con_2(5,:),'LineWidth',2)

legend('LTI Simulation','DMD prediction')
xlabel('x_H')
ylabel('x_H dot')



%% BUILD zeta_1_dot

for i=1:length(t_2)

    zeta_1_dot_rebuilt(i) = build_z1d(t_2(i),x_rebuilt_con_2(:,i),OM);
    zeta_1_dot(i) = build_z1d(t_2(i),x_2(:,i),OM);
end

figure
plot(t_2,zeta_1_dot,'o')
hold on
grid on
plot(t_2,zeta_1_dot_rebuilt,'LineWidth',2)
legend('LTI Simulation','DMD prediction','Start of DMD prediction')
xlabel('t')
ylabel('zeta_1 dot')


%% EIGENVALUES and MODES

figure
scatter(real(Omega(1,1)),imag(Omega(1,1)),'LineWidth',2,'DisplayName','1')
hold on 
grid on
scatter(real(Omega(3,3)),imag(Omega(3,3)),'LineWidth',2,'DisplayName','2')
scatter(real(Omega(5,5)),imag(Omega(5,5)),'LineWidth',2,'DisplayName','3')
scatter(real(Omega(7,7)),imag(Omega(7,7)),'LineWidth',2,'DisplayName','4')
scatter(real(Omega(9,9)),imag(Omega(9,9)),'LineWidth',2,'DisplayName','5')
scatter(real(Omega(11,11)),imag(Lambda(11,11)),'LineWidth',2,'DisplayName','6')
legend

%% RANGE OF OMEGA

OM_vect=linspace(1,40,80);
for j=1:length(OM_vect)

% A = rotor_build.state_matrix_A_handles(0,OM_vect(j));
% [eigvect,eigval(j,:)] = eig(A,'vector');
% Generation of data



[t,x] = ode45(@(t,x) odefun(t,x,rotor_build,OM_vect(j)),time,x_0);
x=x';

X = x(:,1:end-1);
X2 = x(:,2:end);
% X = x(:,ceil(length(time)/2):end-1); %phase change with respect to time dynamics, but overall good fit
% X2 = x(:,ceil(length(time)/2)+1:end);

[U,S,V] = svd(X,'econ');

%  Compute DMD (Phi are eigenvectors)

U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,Lambda] = eig(Atilde);
Phi = X2*V*inv(S)*W;


%Rebuild state (continuous)
%Omega = diag(1./delta_T.*log(diag(Lambda)));
Omega = 1/delta_T*log(Lambda);


% EIGENVALUES OF OMEGA vs A

eigval_Omega(j,:) = diag(Omega);

end


figure
plot(real(eigval_Omega(:,1)),imag(eigval_Omega(:,1)),'o','MarkerSize',10,'linewidth',3,'color',color.blue_mat,'DisplayName','DMD eigenvalues')
hold on 
grid on
plot(real(eigval_Omega(:,2:end)),imag(eigval_Omega(:,2:end)),'o','MarkerSize',10,'linewidth',3,...
    'color',color.blue_mat,'HandleVisibility','off')
legend

%
figure
plot(OM_vect,real(eigval_Omega(:,1)),'o','MarkerSize',5,'linewidth',3,'color',color.blue_mat,'DisplayName','DMD eigenvalues')
hold on 
grid on
plot(OM_vect,real(eigval_Omega(:,2:end)),'o','MarkerSize',5,'linewidth',3,...
    'color',color.blue_mat,'HandleVisibility','off')
legend
ylim([-5 1])
%
figure
plot(OM_vect,real(eigval_Omega(:,1)),'o','MarkerSize',5,'linewidth',3,'color',color.blue_mat,'DisplayName','DMD eigenvalues')
hold on 
grid on
plot(OM_vect,real(eigval_Omega(:,2:end)),'o','MarkerSize',5,'linewidth',3,...
    'color',color.blue_mat,'HandleVisibility','off')
legend
ylim([-0.3 0.2])

%
figure
plot(OM_vect,imag(eigval_Omega(:,1)),'o','MarkerSize',5,'linewidth',3,'color',color.blue_mat,'DisplayName','DMD eigenvalues')
hold on 
grid on
plot(OM_vect,imag(eigval_Omega(:,2:end)),'o','MarkerSize',5,'linewidth',3,...
    'color',color.blue_mat,'HandleVisibility','off')
legend
ylim([0 50])



%% FUNCS
function dx_dt = odefun (t,x,rotor_build,OM)

Z_0_dot = x(1);
Z_1c_dot = x(2);
Z_1s_dot = x(3);
Z_d_dot = x(4);
x_H_dot = x(5);
y_H_dot = x(6);

Z_0= x(7);
Z_1c = x(8);
Z_1s = x(9);
Z_d = x(10);
x_H = x(11);
y_H= x(12);


A =@(t) rotor_build.state_matrix_A_handles(t,OM);

M_NR = [4* rotor_build.problem.rotor_characteristics.blade_inertia_Ib,        0,         0, 0,             0,            0;
            0,        2* rotor_build.problem.rotor_characteristics.blade_inertia_Ib,         0, 0,             0, -2* rotor_build.problem.rotor_characteristics.blade_static_Sb;
            0,        0,         2* rotor_build.problem.rotor_characteristics.blade_inertia_Ib, 0, 2* rotor_build.problem.rotor_characteristics.blade_static_Sb,            0;
            0,        0,         0, 4* rotor_build.problem.rotor_characteristics.blade_inertia_Ib,             0,            0;
            0,        0, 2* rotor_build.problem.rotor_characteristics.blade_static_Sb, 0,              rotor_build.problem.rotor_characteristics.hub_mass_Mx +  rotor_build.problem.rotor_characteristics.blade_mass_Mb*rotor_build.problem.number_blades,            0;
            0, -2* rotor_build.problem.rotor_characteristics.blade_static_Sb,         0, 0,             0,             rotor_build.problem.rotor_characteristics.hub_mass_My +  rotor_build.problem.rotor_characteristics.blade_mass_Mb*rotor_build.problem.number_blades];


M_inv = M_NR\eye(6);

dx_dt = A(t)*x+[M_inv;zeros(6,6)]*nl_damper(t,x,OM);

end

function F1_vect = nl_damper(t,x,OM)

Z_0_dot = x(1);
Z_1c_dot = x(2);
Z_1s_dot = x(3);
Z_d_dot = x(4);

Z_1c = x(8);
Z_1s = x(9);

% use MBC to obtain zeta_1_dot

zeta_1_dot = Z_0_dot + Z_1c_dot*cos(OM*t) + Z_1s_dot*sin(OM*t) - Z_d_dot - Z_1c*OM*sin(OM*t) + Z_1s*OM*cos(OM*t);
% use zeta_1_dot to compute the NL output force
F1 = NL_freeplay(zeta_1_dot,0.3,4067);
% use the MBC to split the force in the four components
F1_vect = [1/4;1/2*cos(OM*t);1/2*sin(OM*t);-1/4;zeros(2,1)].*F1;
% produce a 12*1 output vector to link to the equations


end

function F1 = NL_freeplay (zeta_1_dot,delta,C1)

    if zeta_1_dot<=-delta

        F1 = -C1.*(zeta_1_dot+delta);

    elseif zeta_1_dot<delta&&zeta_1_dot>-delta

        F1 = 0;

    elseif zeta_1_dot>=delta

        F1 = -C1.*(zeta_1_dot-delta);

    end

end

function zeta_1_dot = build_z1d(t,x,OM)

Z_0_dot = x(1);
Z_1c_dot = x(2);
Z_1s_dot = x(3);
Z_d_dot = x(4);

Z_1c = x(8);
Z_1s = x(9);

zeta_1_dot = Z_0_dot + Z_1c_dot*cos(OM*t) + Z_1s_dot*sin(OM*t) - Z_d_dot - Z_1c*OM*sin(OM*t) + Z_1s*OM*cos(OM*t);

end