classdef HD_stability < LTP_stability

    methods

        function obj = HD_stability(rotor_build)
            % Call the superclass constructor
            obj@LTP_stability(rotor_build);
        end

        function eigensolution = HD_single_point(obj,OMEGA,T)

            A = @(t) obj.rotor_build.state_matrix_A_handles(t,OMEGA);

            time = linspace(0,T,200); %move the number of points in the settings

            disp('ok up to here')

            A_HD = stability_analysis.HD_computer(A,time,obj.rotor_build.problem.number_harmonics,OMEGA);

            [eigensolution.eigevectors,eigensolution.eigenvalues] = eig(A_HD,'vector');

        end

        function obj = HD_full_range(obj)

            obj = assign_range_OMEGA(obj);
            obj = assign_period_T(obj);

            for i = 1:obj.rotor_build.problem.number_points

                eigensolution=HD_single_point(obj,obj.modal_solution(i).OMEGA,obj.modal_solution(i).T);

                obj.modal_solution(i).damping=real(eigensolution.eigenvalues);

                obj.modal_solution(i).frequency=imag(eigensolution.eigenvalues);

                obj.modal_solution(i).eigensolution=eigensolution;

            end
        end
    end




end
