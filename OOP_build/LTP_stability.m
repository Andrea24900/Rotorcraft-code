classdef LTP_stability < stability_analysis
    
    methods
        
        function obj = LTP_stability(rotor_build)
            % Call the superclass constructor
            obj@stability_analysis(rotor_build);
        end

        function char_solutions = LTP_single_point(obj,OMEGA,T)

            A = @(t) obj.rotor_build.state_matrix_A_handles(t,OMEGA);
            
            monodromy_matrix_M = stability_analysis.monodromy_computer (A,T);

            [char_solutions.char_multipliers] = eig(monodromy_matrix_M,'vector');

            char_solutions.real_char_exp = 1/(2*T).*log...
                (real(char_solutions.char_multipliers).^2+imag(char_solutions.char_multipliers).^2);

            char_solutions.imag_char_exp = 1/T.*atan2...
                (imag(char_solutions.char_multipliers),real(char_solutions.char_multipliers));

        end

        function obj = LTP_full_range(obj)

            obj = assign_range_OMEGA(obj);

            obj = assign_period_T(obj);

            for i = 1:obj.rotor_build.problem.number_points 
                char_solution=LTP_single_point(obj,obj.modal_solution(i).OMEGA,obj.modal_solution(i).T);

                obj.modal_solution(i).char_solution=char_solution;

                obj.modal_solution(i).damping=char_solution.real_char_exp;

                obj.modal_solution(i).frequency=char_solution.imag_char_exp;
            end
        end
    end

        
end

