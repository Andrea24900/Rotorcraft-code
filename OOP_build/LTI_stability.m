classdef LTI_stability < stability_analysis

    % all methods for the continuation analysis of classic LTI problems, or
    % eventually the analysis of a single OMEGA value
    methods
        
        function obj = LTI_stability(rotor_build)
            % Call the superclass constructor
            obj@stability_analysis(rotor_build);
        end

        function eigensolution = eigen_single_point(obj,OMEGA)

            A = obj.rotor_build.state_matrix_A_handles(1,OMEGA);
            [eigensolution.eigenvectors,eigensolution.eigenvalues] = eig(A,'vector');

        end

        function obj = eigen_full_range(obj)

            obj = assign_range_OMEGA(obj);

            for i = 1:obj.rotor_build.problem.number_points 

                eigensolution=eigen_single_point(obj,obj.modal_solution(i).OMEGA);

                obj.modal_solution(i).damping=real(eigensolution.eigenvalues);

                obj.modal_solution(i).frequency=imag(eigensolution.eigenvalues);
            
                obj.modal_solution(i).eigensolution=eigensolution;

            end

        end
       




    end

end