clear all, close all, clc


%% INTRODUCTION OF NON-LINEAR FUNCTION TO SOLVE


OM=25;

rotor_build = rotor_build.build_all;

%% Generation of data (snapshots) and SVD
delta_T = 1e-2; %NO effect
time = (0:delta_T:5); %best fit ?? (with higher values we get a worse result)
x_0 = ones(12,1); %probable effect of initial condition
[t,x] = ode45(@(t,x) odefun(t,x,rotor_build,OM),time,x_0);
x=x';

X = x(:,1:end-1);
X2 = x(:,2:end);
[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)
r = 10;  % right now best value
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

figure
plot(t,err)
grid on



%% STATE PREDICTION

time = (0:delta_T:10);
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
xline(1,'LineWidth',2)
legend('LTI Simulation','DMD prediction','Start of DMD prediction')
text(1.5,-3.5, sprintf('FIT (all states) = %d percent', FIT_2),'FontSize',15)
xlabel('t')
ylabel('x_H')


figure
plot(t_2,err_2)
grid on




%% RANGE OF OMEGA

OM_vect=linspace(1,40,39);
for j=1:length(OM_vect)

% A = rotor_build.state_matrix_A_handles(0,OM_vect(j));
% [eigvect,eigval(j,:)] = eig(A,'vector');
%% Generation of data



[t,x] = ode45(@(t,x) odefun(t,x,rotor_build,OM_vect(j)),time,x_0);
x=x';

X = x(:,1:end-1);
X2 = x(:,2:end);
[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)

U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,Lambda] = eig(Atilde);
Phi = X2*V*inv(S)*W;


%% Rebuild state (continuous)

%Omega = diag(1./delta_T.*log(diag(Lambda)));
Omega = 1/delta_T*log(Lambda);


%% EIGENVALUES OF OMEGA vs A

eigval_Omega(j,:) = diag(Omega);

end


figure
plot(real(eigval_Omega(:,1)),imag(eigval_Omega(:,1)),'o','MarkerSize',10,'linewidth',3,'color',color.blue_mat,'DisplayName','DMD eigenvalues')
hold on 
grid on
plot(real(eigval_Omega(:,2:end)),imag(eigval_Omega(:,2:end)),'o','MarkerSize',10,'linewidth',3,...
    'color',color.blue_mat,'HandleVisibility','off')
% plot(real(eigval(:,1)),imag(eigval(:,1)),'*','MarkerSize',6,'linewidth',4,'DisplayName','LTI eigenvalues','color',color.orange_mat)
% plot(real(eigval(:,2:end)),imag(eigval(:,2:end)),'*','MarkerSize',6,'linewidth',4,'handleVisibility','off' ,'color',color.orange_mat)
legend

%%
figure
plot(OM_vect,real(eigval_Omega(:,1)),'o','MarkerSize',10,'linewidth',3,'color',color.blue_mat,'DisplayName','DMD eigenvalues')
hold on 
grid on
plot(OM_vect,real(eigval_Omega(:,2:end)),'o','MarkerSize',10,'linewidth',3,...
    'color',color.blue_mat,'HandleVisibility','off')
% plot(OM_vect,real(eigval(:,1)),'*','MarkerSize',8,'linewidth',4,...
%     'color',color.orange_mat,'DisplayName','LTI eigenvalues')
% plot(OM_vect,real(eigval(:,2:end)),'*','MarkerSize',8,'linewidth',4,...
%     'color',color.orange_mat,'HandleVisibility','off')
legend
ylim([-5 1])
%%
figure
plot(OM_vect,imag(eigval_Omega(:,1)),'o','MarkerSize',10,'linewidth',3,'color',color.blue_mat,'DisplayName','DMD eigenvalues')
hold on 
grid on
plot(OM_vect,imag(eigval_Omega(:,2:end)),'o','MarkerSize',10,'linewidth',3,...
    'color',color.blue_mat,'HandleVisibility','off')
% plot(OM_vect,imag(eigval(:,1)),'*','MarkerSize',8,'linewidth',4,...
%     'color',color.orange_mat,'DisplayName','LTI eigenvalues')
% plot(OM_vect,imag(eigval(:,2:end)),'*','MarkerSize',8,'linewidth',4,...
%     'color',color.orange_mat,'HandleVisibility','off')
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

dx_dt = A(t)*x;

end