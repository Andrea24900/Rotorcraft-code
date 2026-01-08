function [phi_jkn,plot_handles] = modal_participation(A_handle_all,char_exponents,eigenvectors,...
    Nh,state_names,state_numbers,mode_names,modes_to_plot,colors_to_plot,plot_type,OMEGA,par,varargin)

% this function produces values for the modal participation and relative plots
% inputs:
% - A_handle_all, function handle of A in all necessary variables
% (t,par,OMEGA)
% - char_exponents, these are the eta_k, computed from the monodromy matrix
% - eigenvectors, these are the V_k eigenvectors from the monodromy matrix,
% in 3d notation
% - Nh, number of required harmonics
% - state_names, used in plots
% - state_numbers, the states to plot
% - mode_names, used in plots
% - modes_to_plot, used to determine which modes are to be plotted,
% numbered as in the original problem
% - colors_to_plot, used to assign colors
% - OMEGA, the frequency of the oscillations of the problem
% - par (optional), the value for the parameter in the problem


% determination of the period T
period_T = 2*pi/OMEGA;
% computation of the LTP state matrix A(t)
if nargin == 12
    state_matrix_A = @(t) A_handle_all(t,par,OMEGA);
elseif nargin == 11
    state_matrix_A = @(t) A_handle_all(t,OMEGA);
end

% determination of the size of A
size_A = width(state_matrix_A(0));

% extraction of eta_M eigenvalues from eigenvalues

eta_M = char_exponents;

% creation of exponential matrix function handle

exp_mat = @(t) expm(-diag(eta_M).*t);


% compute phi(t) for a set of time instants, as in monodromy_computer

phi_0_matrix = eye(size_A);

if period_T > 1
    ode_options=odeset('RelTol',1e-7,'AbsTol',1e-9);
else
    ode_options = odeset('RelTol',1e-5,'AbsTol',1e-7);
end

time_phi = linspace(0,period_T,500);
%for t = 2:length(time_phi)
for j=1:size_A

    [time,phi] = ode78 (@(t,phi) odefun(t,phi,state_matrix_A),[0 period_T],phi_0_matrix(:,j),ode_options);
    
    phi = phi'; % phi is composed by number_states rows and number_time columns
    
    % we want to create a 3d matrix where each row is a state, each columns
    % is one of the integrations, each layer is one time instant

    %phi_collection(:,j,t) = phi(:,end);
    for i=1:size(phi,1)
        phi_interpolated=interp1(time,phi(i,:),time_phi);
        for t = 1:length(time_phi)
        phi_collection(i,j,t) = phi_interpolated(t);
        end
    end

end
%end

%phi_collection(:,:,1) = phi_0_matrix;
% computation of all time varying eigenvectors V(t)

for i=1:length(time_phi)

    V_collection(:,:,i) = phi_collection(:,:,i) * eigenvectors * exp_mat(time_phi(i));

end

% check that V(t) at 0 and T is equal to the eigenvector of the monodromy
% matrix

% eigenvectors
% V_collection(:,:,1)
% V_collection(:,:,end)



%% computation of the fourier series coefficients of V(t)


for j=1:size_A %j-th state, follows the rows of V(t)
    for k=1:size_A % k-th mode, follows the columns of V(t)

        c_jk = complex_coeff (V_collection(j,k,:),Nh);

        c_jk_matrix (j,k,:) = c_jk;


    end
end


%% normalisation and computation of modal participation factors

c_jk_matrix = abs(c_jk_matrix);

N_harmonics = 2*Nh+1;

for j=1:size_A
    for k=1:size_A
        for n=1:N_harmonics
            phi_jkn (j,k,n) = c_jk_matrix (j,k,n) / sum(c_jk_matrix(j,k,:));
        end
    end
end


plot_handles = plot_call (Nh,state_names,state_numbers,mode_names,modes_to_plot,colors_to_plot,plot_type,phi_jkn);

end









function dphi_dt=odefun (t,phi,state_matrix_A)
% used to simulate phi
dphi_dt =  state_matrix_A(t)*phi;

end



