classdef stability_analysis

    % this class contains the main functions of the various stability
    % analyses
    properties
        rotor_build
        modal_solution
    end

    methods (Static)

        %method to run the full stability analysis according to the file
        %definition, and also to build the constructor for the subclass
        function obj = run_stability(rotor_build)


            if rotor_build.problem.required_solver == "LTI"

                if rotor_build.problem.continuation == "YES"

                    obj = continuation_analysis(rotor_build);

                    obj = continuation(obj);

                else

                    obj = LTI_stability(rotor_build);

                    obj = eigen_full_range(obj);

                end

            elseif rotor_build.problem.required_solver == "LTP"

                if rotor_build.problem.continuation == "YES"

                    obj = continuation_analysis(rotor_build);

                    obj = continuation(obj);

                else

                    obj = LTP_stability(rotor_build);

                    obj = LTP_full_range(obj);

                end

            elseif rotor_build.problem.required_solver == "HD"

                if rotor_build.problem.continuation == "YES"

                    obj = continuation_analysis(rotor_build);

                    obj = continuation(obj);

                else

                    obj = HD_stability(rotor_build);

                    obj = HD_full_range(obj);

                end

            end

        end

        function monodromy_matrix_M = monodromy_computer (A_handle_time,period_T)
            %% this function computes the monodromy matrix for a given LTP matrix A
            % inputs: - A_handle_time: function handle of A, as a function of time
            %         - period T of the oscillations in A
            % outputs: monodromy matrix M

            %the column phi(T), corresponding to the initial condition x_O(k) = 1 for
            %k=1 and 0 elsewhere, is equal to the i-th column of M

            phi_0_matrix = eye(width(A_handle_time(0)));
            monodromy_matrix_M = zeros(width(A_handle_time(0)));
            % increased precision for greater accuracy in longer simulations
            %if period_T > 1
              %  ode_options=odeset('RelTol',1e-6,'AbsTol',1e-7);
            %else
             %   ode_options = odeset('RelTol',1e-4,'AbsTol',1e-5);
            %end

            for j=1:width(phi_0_matrix)
                % phi is integrated over one period to find the monodromy matrix =
                % phi(T)
                [~,phi] = ode45(@(t,phi) stability_analysis.odefun_monodromy(t,phi,A_handle_time),[0 period_T],phi_0_matrix(:,j));

                phi=phi';
                % the final value of the time realisation is extracted from phi
                x_T = phi(:,end);
                % this corresponds to the monodromy matrix
                monodromy_matrix_M(:,j)=x_T;
            end
        end

        function [dphi_dt] = odefun_monodromy(t,phi,A)
            % used locally to find phi
            dphi_dt =  A(t)*phi;

        end

        function A_HD = HD_computer (A_handle_time,time,number_harmonics,OMEGA)

            state_number = width(A_handle_time(0));
            size_A_HD = (1+2*number_harmonics)*state_number;
            A_HD = zeros(size_A_HD);



            %% compute fft and determine coefficients for all elements of A
            for instants=1:length(time)

                time_realisation_A(:,:,instants) = A_handle_time (time(instants));

            end


            for j = 1:height(time_realisation_A)
                for k=1:width(time_realisation_A)

                    full_transform(j,k,:) = fft(time_realisation_A(j,k,:))/length(time_realisation_A);

                    A_coeff(j,k,1) = full_transform(j,k,1);

                    index_transform = 2;
                    index_A_coeff = 2;

                    while index_transform <= size(full_transform,3)

                        A_coeff(j,k,index_A_coeff) = 2*real(full_transform(j,k,index_transform));

                        A_coeff(j,k,index_A_coeff+1) = -2*imag(full_transform(j,k,index_transform));

                        index_transform = index_transform + 1;

                        index_A_coeff = index_A_coeff + 2;
                    end

                end
            end



            %assignment of first diagonal block
            A_HD(1:state_number,1:state_number) = stability_analysis.H0M_assign (A_coeff);

            %assignment of first row
            col_A_HD = state_number+1;
            i=1;
            while col_A_HD <= size_A_HD

                A_HD (1:state_number,col_A_HD:col_A_HD+state_number-1) = stability_analysis.H0MiC_assign (A_coeff,i);
                A_HD (1:state_number,col_A_HD+state_number:col_A_HD+2*state_number-1) = stability_analysis.H0MiS_assign (A_coeff,i);
                i = i+1;
                col_A_HD = col_A_HD+2*state_number; %IMPORTANT: multiplication by 2 is motivated by the presence of
                % cos and sin coefficients, which represent two blocks

            end

            %assignment of first column
            row_A_HD = state_number+1;
            i=1;
            while row_A_HD <= size_A_HD

                A_HD (row_A_HD:row_A_HD+state_number-1,1:state_number) = stability_analysis.HiCM_assign (A_coeff,i);
                A_HD (row_A_HD+state_number:row_A_HD+2*state_number-1,1:state_number) = stability_analysis.HiSM_assign (A_coeff,i);
                i=i+1;
                row_A_HD = row_A_HD +2*state_number;
            end

            % assignment of diagonal terms
            col_A_HD = state_number+1;
            row_A_HD = state_number+1;
            i=1;

            while row_A_HD <= size_A_HD % verified for rows, it's the same for columns

                % upper left corner element: HiCMiC
                A_HD (row_A_HD:row_A_HD+state_number-1,col_A_HD:col_A_HD+state_number-1) = stability_analysis.HiCMjC_assign (A_coeff,i,i);
                % upper right corner element: HiCMiS - i*OMEGA
                A_HD (row_A_HD:row_A_HD+state_number-1,col_A_HD+state_number:col_A_HD+2*state_number-1) = ...
                    stability_analysis.HiCMjS_assign (A_coeff,i,i) - eye(state_number)*i*OMEGA;
                % lower left corner element: HiSFiC + i*OMEGA
                A_HD (row_A_HD+state_number:row_A_HD+2*state_number-1,col_A_HD:col_A_HD+state_number-1) =...
                    stability_analysis.HiSMjC_assign (A_coeff,i,i) + eye(state_number)*i*OMEGA;
                % lower right corner element: HiSFiS
                A_HD (row_A_HD+state_number:row_A_HD+2*state_number-1,col_A_HD+state_number:col_A_HD+2*state_number-1) =...
                    stability_analysis.HiSMjS_assign (A_coeff,i,i);
                i=i+1;
                row_A_HD = row_A_HD +2*state_number;
                col_A_HD = col_A_HD + 2*state_number;


            end

            %assignment of extra-diagonal terms (matrix scanned by columns, then by
            % rows)

            col_A_HD = 3 * state_number+1; % directly skip first diagonal block
            row_A_HD = state_number+1;
            i=1; % index for "block rows"
            j=2; % index for "block columns"

            while row_A_HD <= size_A_HD % verified for rows (external loop)

                while col_A_HD <=size_A_HD % internal column loop
                    if i~=j
                        % upper left corner element: HiCMjC
                        A_HD (row_A_HD:row_A_HD+state_number-1,col_A_HD:col_A_HD+state_number-1) = ...
                            stability_analysis.HiCMjC_assign (A_coeff,i,j);
                        % upper right corner element: HiCMjS
                        A_HD (row_A_HD:row_A_HD+state_number-1,col_A_HD+state_number:col_A_HD+2*state_number-1) = ...
                            stability_analysis.HiCMjS_assign (A_coeff,i,j);
                        % lower left corner element: HiSFjC
                        A_HD (row_A_HD+state_number:row_A_HD+2*state_number-1,col_A_HD:col_A_HD+state_number-1) =...
                            stability_analysis.HiSMjC_assign (A_coeff,i,j);
                        % lower right corner element: HiSFjS
                        A_HD (row_A_HD+state_number:row_A_HD+2*state_number-1,col_A_HD+state_number:col_A_HD+2*state_number-1) =...
                            stability_analysis.HiSMjS_assign (A_coeff,i,j);
                    end
                    j = j+1;
                    col_A_HD = col_A_HD + 2*state_number;
                end
                col_A_HD = state_number+1; % restart from first block column in the second block row
                j=1; %restart from first set of coefficients in the next row
                i=i+1;
                row_A_HD = row_A_HD +2*state_number;


            end

        function A_HB = HD_computer_complex (A_handle_time,time,number_harmonics,OMEGA)

        state_number = width(A_handle_time(0));
        size_A_HB = (1+2*number_harmonics)*state_number;
        A_HB = zeros(size_A_HB);

        %% compute fft
            for instants=1:length(time)

                time_realisation_A(:,:,instants) = A_handle_time (time(instants));

            end


            for j = 1:height(time_realisation_A)
                for k=1:width(time_realisation_A)

                    A_coeff(j,k,:) = fft(time_realisation_A(j,k,:))/length(time_realisation_A);

                 end
            end

            % k is used to cycle the harmonics from -N to N (rows)
            % j is used to cycle the blocks in the A matrix (columns)

            k = -number_harmonics;
            j = -number_harmonics;

            for block_row = 0:2*number_harmonics

              for block_column = 0:2*number_harmonics

                m = j-k;

                if m<0
                  A_m = A_coeff(:,:,end+1+m);
                elseif m==0
                  A_m = A_coeff(:,:,1)-1i*k*OMEGA*eye(state_number);
                elseif m>0
                  A_m = A_coeff(:,:,m+1);
                end

                A_HB(block_row*state_number+1:(block_row+1)*state_number,block_column*state_number+1:(block_column+1)*state_number) = A_m;

                % cycle the column index
                k = k+1;

              endfor

               % cycle the row index
                j = j+1;
                %reset the column index
                k=-number_harmonics;
            end

        end

        end




        %% H assignment functions
        % H0 terms
        function H_0_M = H0M_assign (M_coeff)

            H_0_M = M_coeff(:,:,1);

        end

        function H_0_M_iC = H0MiC_assign (M_coeff,i)

            H_0_M_iC = M_coeff(:,:,2*i)./2;

        end

        function H_0_M_iS = H0MiS_assign (M_coeff,i)

            H_0_M_iS = M_coeff(:,:,2*i+1)./2;

        end
        %HiC/SM terms
        function H_iC_M = HiCM_assign (M_coeff,i)

            H_iC_M = M_coeff(:,:,2*i);

        end

        function H_iS_M = HiSM_assign (M_coeff,i)

            H_iS_M = M_coeff(:,:,2*i+1);

        end
        %HiCMjC
        function H_iC_M_jC = HiCMjC_assign (M_coeff,i,j)

            if i==j
                k=i+j;
                H_iC_M_jC = M_coeff(:,:,1) + M_coeff(:,:,k*2)./2;
            elseif i>j
                k=i+j;
                l=i-j;
                H_iC_M_jC = (M_coeff(:,:,l*2) + M_coeff(:,:,k*2))./2;
            elseif j>i
                k=i+j;
                m=j-i;
                H_iC_M_jC = (M_coeff(:,:,m*2) + M_coeff(:,:,k*2))./2;
            end
        end
        %HiCMjS
        function H_iC_M_jS = HiCMjS_assign (M_coeff,i,j)

            if i==j
                k=i+j;
                H_iC_M_jS = M_coeff(:,:,k*2+1)./2;
            elseif i>j
                k=i+j;
                l=i-j;
                H_iC_M_jS = (M_coeff(:,:,k*2+1) - M_coeff(:,:,l*2+1))./2;
            elseif j>i
                k=i+j;
                m=j-i;
                H_iC_M_jS = (M_coeff(:,:,m*2+1) + M_coeff(:,:,k*2+1))./2;
            end
        end

        %HiSMjC
        function H_iS_M_jC = HiSMjC_assign (M_coeff,i,j)

            if i==j
                k=i+j;
                H_iS_M_jC = M_coeff(:,:,k*2+1)./2;
            elseif i>j
                k=i+j;
                l=i-j;
                H_iS_M_jC = (M_coeff(:,:,k*2+1) + M_coeff(:,:,l*2+1))./2;
            elseif j>i
                k=i+j;
                m=j-i;
                H_iS_M_jC = (M_coeff(:,:,k*2+1) - M_coeff(:,:,m*2+1))./2;
            end
        end

        %HiSMjS
        function H_iS_M_jS = HiSMjS_assign (M_coeff,i,j)

            if i==j
                k=i+j;
                H_iS_M_jS = M_coeff(:,:,1) - M_coeff(:,:,k*2)./2;
            elseif i>j
                k=i+j;
                l=i-j;
                H_iS_M_jS = (M_coeff(:,:,l*2) - M_coeff(:,:,k*2))./2;
            elseif j>i
                k=i+j;
                m=j-i;
                H_iS_M_jS = (M_coeff(:,:,m*2) - M_coeff(:,:,k*2))./2;
            end
        end

    end

    methods

        function obj = stability_analysis(rotor_build)
            % this is the constructor of the class
            obj.rotor_build = rotor_build;
        end

        function obj = assign_range_OMEGA(obj)

            rotor_OMEGA_vector=linspace(obj.rotor_build.problem.lower_rotor_RPM/60*2*pi,...
                obj.rotor_build.problem.higher_rotor_RPM/60*2*pi,obj.rotor_build.problem.number_points);
            rotor_OMEGA_vector_RPM = linspace(obj.rotor_build.problem.lower_rotor_RPM,...
                obj.rotor_build.problem.higher_rotor_RPM,obj.rotor_build.problem.number_points);

            if obj.rotor_build.problem.solution_direction == "R2L"

                rotor_OMEGA_vector = fliplr (rotor_OMEGA_vector);
                rotor_OMEGA_vector_RPM = fliplr (rotor_OMEGA_vector_RPM);
            end
            for i = 1:obj.rotor_build.problem.number_points
                obj.modal_solution(i).OMEGA = rotor_OMEGA_vector(i);
                obj.modal_solution(i).OMEGA_RPM = rotor_OMEGA_vector_RPM(i);
            end

        end

        function obj = assign_period_T(obj)
            for i = 1:obj.rotor_build.problem.number_points
                obj.modal_solution(i).T = 2*pi/obj.modal_solution(i).OMEGA;
            end
        end


    end

end
