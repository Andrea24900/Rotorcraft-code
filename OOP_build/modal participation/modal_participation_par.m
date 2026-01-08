function phi_1mode1state = modal_participation_par(A_handle_all,char_exponents,eigenvectors,Nh,mode,state,OMEGA,par,varargin)

% this function is equal to modal_participation, but does not produce plots
% inputs:
% - par (optional), the value for the parameter in the problem
% - A_handle_all, function handle of A in all necessary variables
% (t,par,OMEGA) 
% - char_exponents, these are the eta_k, computed from the monodromy matrix
% - eigenvectors, these are the V_k eigenvectors from the monodromy matrix,
% in 3d notation
% - Nh, number of required harmonics
% - mode, required mode
% - state, required state
% - OMEGA, the frequency of the oscillations of the problem

% determination of the period T 
period_T = 2*pi/OMEGA;
% computation of the LTP state matrix A(t)

if nargin == 8
state_matrix_A = @(t) A_handle_all(t,par,OMEGA);
elseif nargin == 7
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

time_step = period_T/1000;

time_phi = (0:time_step:period_T);

for j=1:size_A

    [time,phi] = ode78 (@(t,phi) odefun(t,phi,state_matrix_A),[0 period_T],phi_0_matrix(:,j),ode_options);

    phi = phi'; % phi is composed by number_states rows and number_time columns

    % we want to create a 3d matrix where each row is a state, each columns
    % is one of the integrations, each layer is one time instant
    for i=1:height(phi)
        phi_interpolated(i,:)=interp1(time,phi(i,:),time_phi);
    end
    phi_collection(:,j,:) = phi_interpolated(:,:);

end

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

    
for n=1:N_harmonics
    phi_1mode1state (n) = phi_jkn (state,mode,n);
end
phi_1mode1state = phi_1mode1state';

end


function dphi_dt=odefun (t,phi,state_matrix_A)
% used to simulate phi
dphi_dt =  state_matrix_A(t)*phi;

end
