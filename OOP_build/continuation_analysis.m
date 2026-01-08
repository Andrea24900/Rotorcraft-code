classdef continuation_analysis < stability_analysis

    methods

        function obj = continuation_analysis(rotor_build)
            % Call the superclass constructor
            obj@stability_analysis(rotor_build);
        end


        function obj=continuation(obj)

            obj = assign_range_OMEGA(obj);

            switch obj.rotor_build.problem.required_solver

                case 'LTI'

                    A_OMEGA = @(OMEGA) obj.rotor_build.state_matrix_A_handles(0,OMEGA);

                    A = A_OMEGA(obj.modal_solution(1).OMEGA);

                    problem_size = size(A,1);

                    d_M_d_OMEGA = continuation_analysis.numerical_sensitivity...
                        ('SS','CEN',A_OMEGA,obj.rotor_build.problem.step_h,obj.modal_solution(1).OMEGA);

                    [eigenvectors, eigenvalues] = eig(A,'vector');

                case 'LTP'

                    T = 2*pi/obj.modal_solution(1).OMEGA;

                    A_time_OMEGA = @(t,OMEGA) obj.rotor_build.state_matrix_A_handles(t,OMEGA);

                    A_time =@(t) A_time_OMEGA(t,obj.modal_solution(1).OMEGA);

                    A = continuation_analysis.monodromy_computer (A_time,T);

                    problem_size = size(A,1);

                    d_M_d_OMEGA = continuation_analysis.numerical_sensitivity...
                        ('MON','CEN',A_time_OMEGA,obj.rotor_build.problem.step_h,...
                        obj.modal_solution(1).OMEGA);

                    [eigenvectors, eigenvalues] = eig(A,'vector');

                case 'HD'

                    T = 2*pi/obj.modal_solution(1).OMEGA;
                    
                    time = linspace(0,T,obj.rotor_build.problem.time_samples);

                    A_time_OMEGA = @(t,OMEGA) obj.rotor_build.state_matrix_A_handles(t,OMEGA);

                    A_time =@(t) A_time_OMEGA(t,obj.modal_solution(1).OMEGA);

                    A = continuation_analysis.HD_computer...
                        (A_time,time,obj.rotor_build.problem.number_harmonics,...
                        obj.modal_solution(1).OMEGA);

                    problem_size = size(A,1);

                    d_M_d_OMEGA = continuation_analysis.numerical_sensitivity...
                        ('HB','CEN',A_time_OMEGA,obj.rotor_build.problem.step_h,...
                        obj.modal_solution(1).OMEGA,1000,...
                        [obj.rotor_build.problem.time_samples obj.rotor_build.problem.number_harmonics]);

                    [eigenvectors, eigenvalues] = eig(A,'vector');


            end

            [eigenvalues,indices]=sort(eigenvalues);

            eigenvectors=eigenvectors(:,indices);

            obj.modal_solution(1).eigensolution.eigenvalues =  eigenvalues;

            obj.modal_solution(1).eigensolution.eigenvectors = eigenvectors;

            if obj.rotor_build.problem.required_solver == "LTP"

            obj.modal_solution(1).damping = 1/(2*T).*log...
                (real(eigenvalues).^2+imag(eigenvalues).^2);

            obj.modal_solution(1).frequency = 1/T.*atan2...
                (imag(eigenvalues),real(eigenvalues));


            else

            obj.modal_solution(1).damping =  real(eigenvalues);

            obj.modal_solution(1).frequency =  imag(eigenvalues);
            
            end
            

            % STEP 1.5: find next eigenvalues through SENS
            delta_OMEGA = abs(obj.modal_solution(2).OMEGA-obj.modal_solution(1).OMEGA);


            for j=1:problem_size

                current_eigenvalue=obj.modal_solution(1).eigensolution.eigenvalues(j,1);

                current_eigenvector=obj.modal_solution(1).eigensolution.eigenvectors(:,j);

                A_sens = [current_eigenvector eye(problem_size).*current_eigenvalue-A;...
                    0 2*current_eigenvector'];
                b_sens = [d_M_d_OMEGA*current_eigenvector;0];

                if rcond(A_sens) < eps
                    x_sens=pinv(A_sens) * b_sens;
                else
                    x_sens = A_sens\b_sens;
                end

                next_eigenvalue(j,1) = x_sens(1).*delta_OMEGA+current_eigenvalue;

                next_eigenvector(:,j) = x_sens(2:end).*delta_OMEGA+current_eigenvector;

            end
            
            %%
            % STEP 2: use iteration on all omegas
            for i = 2:obj.rotor_build.problem.number_points
                i
                
                switch obj.rotor_build.problem.required_solver

                    case 'LTI'

                        A = A_OMEGA(obj.modal_solution(i).OMEGA);

                        d_M_d_OMEGA = continuation_analysis.numerical_sensitivity...
                            ('SS','CEN',A_OMEGA,obj.rotor_build.problem.step_h,obj.modal_solution(i).OMEGA);
                        
                    case 'LTP'

                    T = 2*pi/obj.modal_solution(i).OMEGA;

                    A_time_OMEGA = @(t,OMEGA) obj.rotor_build.state_matrix_A_handles(t,OMEGA);

                    A_time =@(t) A_time_OMEGA(t,obj.modal_solution(i).OMEGA);

                    A = continuation_analysis.monodromy_computer (A_time,T);

                    d_M_d_OMEGA = continuation_analysis.numerical_sensitivity...
                        ('MON','CEN',A_time_OMEGA,obj.rotor_build.problem.step_h,...
                        obj.modal_solution(i).OMEGA);

                    case 'HD'
                    T = 2*pi/obj.modal_solution(i).OMEGA;

                    time = linspace(0,T,obj.rotor_build.problem.time_samples);

                    A_time_OMEGA = @(t,OMEGA) obj.rotor_build.state_matrix_A_handles(t,OMEGA);

                    A_time =@(t) A_time_OMEGA(t,obj.modal_solution(i).OMEGA);

                    A = continuation_analysis.HD_computer...
                        (A_time,time,obj.rotor_build.problem.number_harmonics,...
                        obj.modal_solution(i).OMEGA);

                    problem_size = size(A,1);

                    d_M_d_OMEGA = continuation_analysis.numerical_sensitivity...
                        ('HB','CEN',A_time_OMEGA,obj.rotor_build.problem.step_h,...
                        obj.modal_solution(i).OMEGA,1000,...
                        [obj.rotor_build.problem.time_samples obj.rotor_build.problem.number_harmonics]);

                end
                % use the CON matrices to iterate from the current eigenvalue to its
                % precise value

                for j = 1:problem_size

                    % set iterative values as initial current guesses, indicated as
                    % NEXT from previous step

                    iter_eigenvalue = next_eigenvalue(j);

                    iter_eigenvector = next_eigenvector(:,j);

                    %determine residual of perturbation equation (right hand side)

                    residual_b_con=[-iter_eigenvalue*iter_eigenvector+A*iter_eigenvector;...
                        1 - (iter_eigenvector'*iter_eigenvector)];

                    iter_norm = norm([iter_eigenvalue ;iter_eigenvector]);

                    iter_number=0;

                    % start iterating the continuation equation in the same omega until
                    % no further perturbation is needed

                    while norm(residual_b_con)/iter_norm>obj.rotor_build.problem.continuation_tolerance...
                            &&iter_number<obj.rotor_build.problem.continuation_max_iter

                        %find specific A_con
                        A_con = [iter_eigenvector eye(problem_size).*iter_eigenvalue-A;...
                            0 2*iter_eigenvector'];

                        if rcond(A_con) < eps
                            x_con=pinv(A_con) * residual_b_con;
                        else
                            x_con = A_con\residual_b_con;
                        end
                        %update the iterative values
                        iter_eigenvalue = x_con(1)+iter_eigenvalue;

                        iter_eigenvector = x_con(2:end)+iter_eigenvector;

                        iter_number=iter_number+1;

                        residual_b_con=[-iter_eigenvalue*iter_eigenvector+A*iter_eigenvector;...
                            1 - (iter_eigenvector'*iter_eigenvector)];

                    end
                    %save current (i-th omega value) updated eigenvalues in the main
                    %vectors
                    
                    obj.modal_solution(i).eigensolution.eigenvalues(j,1) =  iter_eigenvalue;
                    
                    obj.modal_solution(i).eigensolution.eigenvectors(:,j) = iter_eigenvector;

                    if obj.rotor_build.problem.required_solver == "LTP"

                        obj.modal_solution(i).damping(j,1) = 1/(2*T).*log...
                            (real(iter_eigenvalue).^2+imag(iter_eigenvalue).^2);

                        obj.modal_solution(i).frequency(j,1) = 1/T.*atan2...
                            (imag(iter_eigenvalue),real(iter_eigenvalue));


                    else

                        obj.modal_solution(i).damping(j,1) =  real(iter_eigenvalue);
 
                        obj.modal_solution(i).frequency(j,1) =  imag(iter_eigenvalue);

                    end 

                    %compute the updated term for the next iteration through the SENS
                    %equation
                    if i<obj.rotor_build.problem.number_points

                        current_eigenvalue=iter_eigenvalue;

                        current_eigenvector=iter_eigenvector;

                        A_sens = [current_eigenvector eye(problem_size).*current_eigenvalue-A;...
                            0 current_eigenvector'];
                        b_sens = [d_M_d_OMEGA*current_eigenvector;0];

                        if rcond(A_sens) < eps
                            x_sens=pinv(A_sens) * b_sens;
                        else
                            x_sens = A_sens\b_sens;
                        end

                        next_eigenvalue(j) = x_sens(1).*delta_OMEGA+current_eigenvalue;

                        next_eigenvector(:,j) = x_sens(2:end).*delta_OMEGA+current_eigenvector;

                    end
                end
            end
        end
    end



    methods (Static)
        function sensitivity_M_on_par = numerical_sensitivity(...
                type,type_der,A_handle_all,step_size_h,OMEGA,par,HB_paramaters,varargin)
            % this function computes the numerical sensitivity matrix for the
            % continuation algorithm
            % inputs:
            % - type, choose between SS (state space), HB (harmonic balance), MON
            % (monodromy matrix)
            % - A_handle, the full handle (t,par,OMEGA) or (t,OMEGA) of the LTI or LTP
            % matrix
            % - step_size_h, the h value for the finite forward difference
            % - OMEGA, the frequency of the oscillations
            % - par (optional), the value of the parameter for the perturbation
            % - HB_parameters (only needed if type is "HB"): a vector [NP NH], where NP is the number of points for
            % the harmonic balance computation (lenght of the transform), while NH is
            % the number of harmonics

            switch type
                case 'SS'
                    if nargin == 5 % no par is assigned

                        state_matrix_A_plus_h=A_handle_all(OMEGA+step_size_h);

                        switch type_der
                            case 'FW'
                                state_matrix_A = A_handle_all(OMEGA);
                            case 'CEN'
                                state_matrix_A_min_h = A_handle_all(OMEGA-step_size_h);
                        end
                    elseif nargin == 6 %  par is assigned

                        state_matrix_A_plus_h=A_handle_all(par+step_size_h,OMEGA);
                        switch type_der
                            case 'FW'
                                state_matrix_A = A_handle_all(par,OMEGA);
                            case 'CEN'
                                state_matrix_A_min_h = A_handle_all(par-step_size_h,OMEGA);
                        end
                    end

                    switch type_der
                        case 'FW'
                            sensitivity_M_on_par = 1/(step_size_h).*...
                                (state_matrix_A_plus_h-state_matrix_A);
                        case 'CEN'
                            sensitivity_M_on_par = 1/(2*step_size_h).*...
                                (state_matrix_A_plus_h-state_matrix_A_min_h);
                    end
                case 'HB'
                    number_time_instants=HB_paramaters(1);
                    number_harmonics = HB_paramaters(2);


                    if nargin == 7&&par==1000 % no par is assigned
                        switch type_der
                            case 'FW'
                                %normal sistem
                                T = 2*pi/OMEGA;
                                time = linspace(0,T,number_time_instants);
                                %generation of A(t) from A_handle_all(t,OMEGA)
                                A =@(t) A_handle_all(t,OMEGA);
                                %computation of A_HB
                                A_HB = continuation_analysis.HD_computer(A,time,number_harmonics,OMEGA);
                            case 'CEN'
                                time_min_h = linspace(0,2*pi/(OMEGA-step_size_h),number_time_instants);
                                %generation of A(t) from A_handle_all(t,OMEGA)
                                A_min_h =@(t) A_handle_all(t,OMEGA);
                                %computation of A_HB
                                A_HB_min_h = continuation_analysis.HD_computer(A_min_h,time_min_h,number_harmonics,OMEGA);
                        end
                        % same procedure with OMEGA+h
                        time_plus_h = linspace(0,2*pi/(OMEGA+step_size_h),number_time_instants);
                        A_plus_h =@(t) A_handle_all(t,OMEGA+step_size_h);
                        A_HB_plus_h = continuation_analysis.HD_computer(A_plus_h,time_plus_h,number_harmonics,OMEGA+step_size_h);

                    elseif nargin == 7 && par~=1000%  par is assigned
                        % OMEGA is employed, the variation is par+h
                        T = 2*pi/OMEGA;
                        time = linspace(0,T,number_time_instants);
                        A_plus_h =@(t) A_handle_all(t,par+step_size_h,OMEGA);
                        A_HB_plus_h = continuation_analysis.HD_computer(A_plus_h,time,number_harmonics,OMEGA);

                        switch type_der
                            case 'FW'
                                A =@(t) A_handle_all(t,par,OMEGA);
                                A_HB = continuation_analysis.HD_computer(A,time,number_harmonics,OMEGA);
                            case 'CEN'
                                A_min_h =@(t) A_handle_all(t,par-step_size_h,OMEGA);
                                A_HB_min_h = continuation_analysis.HD_computer(A_min_h,time,number_harmonics,OMEGA);
                        end
                    end
                    %final sens computation
                    switch type_der
                        case 'FW'
                            sensitivity_M_on_par = 1/(step_size_h).*(A_HB_plus_h-A_HB);

                        case 'CEN'
                            sensitivity_M_on_par = 1/(2*step_size_h).*...
                                (A_HB_plus_h-A_HB_min_h);
                    end
                case 'MON'

                    if nargin == 5 % variation in OMEGA
                        % generation of T and T(OMEGA+h)
                        T_plus_h = 2*pi/(OMEGA+step_size_h);
                        state_matrix_A_plus_h=@(t)A_handle_all(t,OMEGA+step_size_h);
                        monodromy_matrix_M_plus_h = continuation_analysis.monodromy_computer (state_matrix_A_plus_h,T_plus_h);

                        switch type_der
                            case 'FW'
                                T = 2*pi/OMEGA;
                                % generation of A(t)
                                state_matrix_A = @(t) A_handle_all(t,OMEGA);
                                % computation of the two monodromy matrices
                                monodromy_matrix_M =continuation_analysis.monodromy_computer (state_matrix_A,T);
                            case 'CEN'
                                T_min_h = 2*pi/(OMEGA-step_size_h);
                                % generation of A(t)
                                state_matrix_A_min_h = @(t) A_handle_all(t,OMEGA-step_size_h);
                                % computation of the two monodromy matrices
                                monodromy_matrix_M_min_h = continuation_analysis.monodromy_computer (state_matrix_A_min_h,T_min_h);
                        end
                    elseif nargin == 6 %same procedure for par
                        T = 2*pi/OMEGA;
                        state_matrix_A_plus_h=@(t)A_handle_all(t,par+step_size_h,OMEGA);
                        monodromy_matrix_M_plus_h = monodromy_computer (state_matrix_A_plus_h,T);

                        switch type_der
                            case 'FW'
                                state_matrix_A = @(t) A_handle_all(t,par,OMEGA);
                                monodromy_matrix_M = monodromy_computer (state_matrix_A,T);
                            case 'CEN'
                                state_matrix_A_min_h = @(t) A_handle_all(t,par-step_size_h,OMEGA);
                                monodromy_matrix_M_min_h = monodromy_computer (state_matrix_A_min_h,T);
                        end
                    end
                    switch type_der
                        case 'FW'
                            sensitivity_M_on_par = 1/(step_size_h).*(monodromy_matrix_M_plus_h-monodromy_matrix_M);
                        case 'CEN'
                            sensitivity_M_on_par = 1/(2*step_size_h).*(monodromy_matrix_M_plus_h-monodromy_matrix_M_min_h);

                    end

            end
        end

    end
end
