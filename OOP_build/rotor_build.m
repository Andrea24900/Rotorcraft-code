classdef rotor_build
    %ROTOR This class contains all methods to build the state-space matrix 
    % for given rotor characteristics
    
    properties
        problem
        state_matrix_A_handles
    end
    
    methods (Static)
        
            function obj=build_all
                run("problem_definition.m")
                % USE THIS TO BUILD THE ROTOR FROM THE PROBLEM DEFINITION
                obj = rotor_build(problem);
                obj = add_A_to_rotor(obj);
            end

    end

    methods
        function obj = rotor_build(rotorStruct)
            %ROTOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.problem = rotorStruct;
        end
        
        % build M,C,K matrices from predetermined symbolic formulation

        function mass_matrix_M = build_mass_matrix(obj)
             %METHOD: build_mass_matrix: this method builds the mass matrix
             % from the rotor characteristics 
             % and according to predetermined symbolic formulations
             mass_matrix_M = mass_matrix_repository(obj);
        end

        function [damping_matrix_C,stiffness_matrix_K] = build_damping_matrix(obj)
             %METHOD: build_damping_matrix: this method builds the damping matrix
             % as a function handle with OMEGA and time dependence from the rotor characteristics 
             % and according to predetermined symbolic formulations
             [damping_matrix_C,stiffness_matrix_K] = damping_matrix_repository(obj);
        end

        % build A matrix
        function state_matrix_A = build_state_matrix(obj)
            %this method produces the state matrix with all possible
            %dependencies (even if they are not needed)
            mass_matrix_M = mass_matrix_repository(obj);
            [damping_matrix_C,stiffness_matrix_K] = damping_matrix_repository(obj);
            state_matrix_A =@(t,rotor_OMEGA) [-mass_matrix_M\eye(width(mass_matrix_M))*damping_matrix_C(t,rotor_OMEGA) -mass_matrix_M\eye(width(mass_matrix_M))*stiffness_matrix_K(t,rotor_OMEGA);
            eye(width(mass_matrix_M)) zeros(width(mass_matrix_M))];
        end

        % add A matrix to problem
        function obj=add_A_to_rotor(obj)
            state_matrix_A = build_state_matrix(obj);
            obj.state_matrix_A_handles = state_matrix_A;
        end
    end
        
        


end

