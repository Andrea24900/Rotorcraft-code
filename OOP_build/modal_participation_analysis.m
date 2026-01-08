classdef modal_participation_analysis

    properties
        modal_participation
        modal_solution
        rotor_build
    end


    methods

        function obj = modal_participation_analysis(modes)
            % Call the superclass constructor
            obj.modal_solution=modes.modal_solution;
            obj.rotor_build=modes.rotor_build;
            %call the full range analysis
            
                obj = modal_part_all(obj);
            
        end


            function obj = modal_part_all(obj)

                for i = 1:max(size(obj.modal_solution))
                    if obj.rotor_build.problem.required_solver == "LTP"
                        phi_jkn_mat = modal_participation_analysis.mod_part_single_LTP(obj,i);
                    elseif obj.rotor_build.problem.required_solver == "HD"
                        phi_jkn_mat = modal_participation_analysis.mod_part_single_HD(obj,i);
                    end
                    
                    obj.modal_participation(i).OMEGA = obj.modal_solution(i).OMEGA;
                    obj.modal_participation(i).OMEGA_RPM = obj.modal_solution(i).OMEGA_RPM;
                    obj.modal_participation(i).phi_jkn=phi_jkn_mat;
                    
                end

            end

    end




        methods (Static)

        function phi_jkn_mat = mod_part_single_LTP(obj,index) 

            % determination of the period T
            T = 2*pi/obj.modal_solution(index).OMEGA;
            
            A = @(t) obj.rotor_build.state_matrix_A_handles(t,obj.modal_solution(index).OMEGA);
            

            % determination of the size of A
            size_A = width(A(0));

            % retrieve characteristic exponents
            eta_M = obj.modal_solution(index).damping + 1i*obj.modal_solution(index).frequency;

            % creation of exponential matrix function handle

            exp_mat = @(t) expm(-diag(eta_M).*t);


            % compute phi(t) for a set of time instants, as in monodromy_computer

            phi_0_matrix = eye(size_A);

            if T > 1
                ode_options=odeset('RelTol',1e-7,'AbsTol',1e-9);
            else
                ode_options = odeset('RelTol',1e-5,'AbsTol',1e-7);
            end

            time_phi = linspace(0,T,500);
            %for t = 2:length(time_phi)
            for j=1:size_A

                [time,phi] = ode78 (@(t,phi) ...
                    modal_participation_analysis.odefun(t,phi,A),[0 T],phi_0_matrix(:,j),ode_options);

                phi = phi'; % phi is composed by number_states rows and number_time columns

                % we want to create a 3d matrix where each row is a state, each columns
                % is one of the integrations, each layer is one time instant

                
                for i=1:size(phi,1)
                    phi_interpolated=interp1(time,phi(i,:),time_phi);
                    for t = 1:length(time_phi)
                        phi_collection(i,j,t) = phi_interpolated(t);
                    end
                end

            end

            for i=1:length(time_phi)

                V_collection(:,:,i) = phi_collection(:,:,i) * ...
                    obj.modal_solution(index).eigensolution.eigenvectors * exp_mat(time_phi(i));

            end
            %% computation of the fourier series coefficients of V(t)

            for j=1:size_A %j-th state, follows the rows of V(t)
                for k=1:size_A % k-th mode, follows the columns of V(t)

                    c_jk = modal_participation_analysis.complex_coeff ...
                        (V_collection(j,k,:),obj.rotor_build.problem.number_harmonics*2);

                    c_jk_matrix (j,k,:) = c_jk;


                end
            end


            %% normalisation and computation of modal participation factors

            c_jk_matrix = abs(c_jk_matrix);

            N_harmonics = 2*obj.rotor_build.problem.number_harmonics*2+1;

            for j=1:size_A
                for k=1:size_A
                    for n=1:N_harmonics
                        phi_jkn_mat (j,k,n) = c_jk_matrix (j,k,n) / sum(c_jk_matrix(j,k,:));
                    end
                end
            end
        end

            function dphi_dt=odefun (t,phi,state_matrix_A)
                % used to simulate phi
                dphi_dt =  state_matrix_A(t)*phi;

            end

            function complex_coeffs = complex_coeff (V_jk_t,N_harmonics)
                % this function uses the fft to determine the complex coefficients of the
                % fourier series for V(t)
                % VALIDATED
                % the fft is computed, then centered
                all_coeff = fftshift(fft(V_jk_t)/length(V_jk_t));
                mid_index = 1 + floor(length(V_jk_t)/2);
                % the required N_harmonics coefficients are assigned
                complex_coeffs = all_coeff(mid_index-N_harmonics:mid_index+N_harmonics);

            end

            function phi_jkn_mat = mod_part_single_HD(obj,index)

            % this process must be performed for all k eigenvectors of the
            % expanded eigenvector matrix
            
            eigenvector_matrix = obj.modal_solution(index).eigensolution.eigenvectors;
            block_size = max(size(eigenvector_matrix))/(2*obj.rotor_build.problem.number_harmonics+1);
            
            % now I want to work on the k-th eigenvector, while cycing for
            % all k
            
            N_harmonics = 2*obj.rotor_build.problem.number_harmonics+1;

            for k=1:max(size(eigenvector_matrix))

                % now I need to extract the components of x_0,x_1c,x_1s, etc
                  X_0_k = eigenvector_matrix(1:block_size,k);
                  
                  % work the remaining components until they are all
                  % extracted

                  j=block_size+1;
                  n = 1;
                  while n<=obj.rotor_build.problem.number_harmonics

                      X_nc_k(:,n) = eigenvector_matrix(j:j+block_size-1, k);  
                      X_ns_k(:,n) = eigenvector_matrix(j+block_size:j+2*block_size-1, k);  

                      j=j+2*block_size;
                      

                      C_k_nplus(:,n) = 1/2*(X_nc_k(:,n)-1i*X_ns_k(:,n));
                      C_k_nminus(:,n) = 1/2*(X_nc_k(:,n)+1i*X_ns_k(:,n));

                      n=n+1;
                  end

                  C_k_nminus = fliplr(C_k_nminus);

                  for i=1:obj.rotor_build.problem.number_harmonics

                  c_jk_matrix(:,k,i)=C_k_nminus(:,i);
                  
                  end

                  c_jk_matrix(:,k,obj.rotor_build.problem.number_harmonics+1)=X_0_k;
                  j=1;
                  for i=obj.rotor_build.problem.number_harmonics+2:N_harmonics

                  c_jk_matrix(:,k,i)=C_k_nplus(:,j);
                  j=j+1;
                  
                  end



                  


            end

            
            %% normalisation and computation of modal participation factors

            c_jk_matrix = abs(c_jk_matrix);

            

            for j=1:block_size
                for k=1:max(size(eigenvector_matrix))
                    for n=1:N_harmonics
                        phi_jkn_mat (j,k,n) = c_jk_matrix (j,k,n) / sum(c_jk_matrix(j,k,:));
                    end
                end
            end

            end

        end
    
end