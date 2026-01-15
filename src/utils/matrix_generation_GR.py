"""Matrix generation functions for mass, damping, and stiffness matrices of a Hammond's
rotor in ground resonance

Based on Hammond's work (1974) and the conventions in Bassi (2024).

Damper activation modes:
- ALL: All dampers active (produces LTI system)
- ODI: One Damper Inoperative (produces LTP system with periodic coefficients)
- MDI: More Dampers Inoperative (User choice through vector of binaries)

Damper attachments: H2B (hub-to-blade) and B2B (blade-to-blade or interblade)

Number of blades: from 3 (minimum number for GR prone rotor) to arbitrary Nb
"""

import numpy as np
from typing import Callable, Tuple
from ..config.problem_definition import ProblemDefinition

def build_mass_matrix(
    problem: ProblemDefinition,mbc_matrices = None
) -> np.ndarray:
    """Build mass matrix for the rotor system. If the MBC matrices are not created, build them.

    Args:
        problem: Problem definition containing rotor characteristics

    Returns:
        mass_matrix, mass_matrix_rot (function handles in Omega,t)

    Note:
        Matrix dimensions depend on number of blades (Nb blade states and 2 hub ones):
        - 3 blades: 5x5 matrix
        - 4 blades: 6x6 matrix
        - 5 blades: 7x7 matrix
    """
    num_blades = problem.number_blades
    rotor_chars = problem.rotor_characteristics

    # Extract rotor characteristics for clarity
    blade_inertia = rotor_chars.blade_inertia_Ib
    blade_mass = rotor_chars.blade_mass_Mb
    blade_static_moment = rotor_chars.blade_static_Sb
    hub_mass_x = rotor_chars.hub_mass_Mx
    hub_mass_y = rotor_chars.hub_mass_My

    # Total hub mass including blades
    total_hub_mass_x = hub_mass_x + blade_mass * num_blades
    total_hub_mass_y = hub_mass_y + blade_mass * num_blades

    def mass_matrix_rot(time: float, rotor_omega:float) -> np.ndarray:
        # Define the four blocks which compose the mass matrix
        # Upper Left
        diagonal = blade_inertia*np.ones(num_blades)
        upper_left = np.diag(diagonal)
        # Lower right
        lower_right = np.array([
            [total_hub_mass_x,0],
            [0,total_hub_mass_y]
            ])
        # Time varying terms (which use a psi azimuth angle time function)
        # the upper right term is a Nb*2 array containing periodic functions
        left_vec=np.zeros(num_blades)
        right_vec=np.zeros(num_blades)
        for blade_index in range(1,num_blades+1):
            psi = azimuth_angle(blade_index,num_blades)
            left_vec[blade_index-1]=np.sin(psi(rotor_omega,time))
            right_vec[blade_index-1]=-np.cos(psi(rotor_omega,time))
        upper_right = blade_static_moment*np.column_stack([left_vec,right_vec])
        # the lower left term is a 2*Nb array
        upper_vec=left_vec
        lower_vec=right_vec
        lower_left = blade_static_moment*np.vstack([upper_vec,lower_vec])
        
        #generation of mass matrix (partially in rotating frame)
        mass_matrix_result = np.block([[upper_left,upper_right],[lower_left,lower_right]])

        return mass_matrix_result

    def mass_matrix(time:float, rotor_omega:float):

        # rotation of the matrix in the proper frame
        mbc_matrix = build_mbc_matrix(problem, time, rotor_omega)
        mass_matrix_result = np.transpose(mbc_matrix)@mass_matrix_rot(time, rotor_omega)@mbc_matrix

        return mass_matrix_result
    
    return mass_matrix,mass_matrix_rot


def build_damping_matrix(
    problem: ProblemDefinition
) -> Tuple[Callable[[float, float], np.ndarray], Callable[[float, float], np.ndarray]]:
    """Build damping and stiffness matrices as function handles.

    Args:
        problem: Problem definition containing rotor characteristics

    Returns:
        Tuple of (damping_matrix_C, stiffness_matrix_K) as function handles
        Both functions take (t, rotor_omega) as arguments and return matrices

    """
    num_blades = problem.number_blades
    rotor_chars = problem.rotor_characteristics

    # Extract commonly used parameters for clarity
    blade_static_moment = rotor_chars.blade_static_Sb
    lag_hinge_offset = rotor_chars.lag_hinge_offset_e
    nominal_damping = rotor_chars.nominal_damping_Cd
    nominal_stiffness = rotor_chars.nominal_stiffness_Kd
    hub_damping_x = rotor_chars.hub_damping_Cx
    hub_damping_y = rotor_chars.hub_damping_Cy
    hub_stiffness_x = rotor_chars.hub_stiffness_Kx
    hub_stiffness_y = rotor_chars.hub_stiffness_Ky

    # add a check for the presence of a mass_matrix_rot from previous computation 

    # if no cached matrix: creation of the mass matrix 
    (mass_matrix,mass_matrix_rot) = build_mass_matrix(problem)

    # working principle of this matrix generator
    # 1) generate C and K without dampers
    # 2) sum the damping and stiffness matrices with the chosen damper configuration
    # 3) transform the resulting matrices through the MBC transformation

    
    def damping_matrix_rot(time:float ,rotor_omega:float):
        # 1) C without dampers
        #block generation of matrix
        upper_left = np.diag(np.zeros(num_blades))
        # Lower right
        lower_right = np.array([
            [hub_damping_x,0],
            [0,hub_damping_y]
            ])
        # Time varying terms (which use a psi azimuth angle time function)
        # the upper right term is a Nb*2 array containing zeros
        upper_right = np.zeros((num_blades,2))
        # the lower left term is a 2*Nb array with periodic terms
        upper_vec=np.zeros(num_blades)
        lower_vec=np.zeros(num_blades)
        for blade_index in range(1,num_blades+1):
            psi = azimuth_angle(blade_index,num_blades)
            upper_vec[blade_index-1]=2*rotor_omega*blade_static_moment*np.cos(psi(rotor_omega,time))
            lower_vec[blade_index-1]=2*rotor_omega*blade_static_moment*np.sin(psi(rotor_omega,time))
        lower_left = np.vstack([upper_vec,lower_vec])
        #generation of mass matrix (partially in rotating frame)
        damping_matrix_partial = np.block([[upper_left,upper_right],[lower_left,lower_right]])
        
        # 2) C with the addition of the dampers
        damping_matrix_addition = np.zeros((num_blades+2,num_blades+2))
        if problem.damper_connection == "H2B":
            diagonal = damping_activation_H2B(problem)
            damping_matrix_addition[0:num_blades,0:num_blades] = np.diag(diagonal)
        #elif problem.damper_connection == "B2B":
            # B2B damper connection logic to be implemented
            
        # summation of matrices
        damping_matrix_result = damping_matrix_partial+damping_matrix_addition
        # 3) transformation through MBC
        return damping_matrix_result

    def damping_matrix_C(time,rotor_omega):
        # build the MBC matrices
        mbc_matrix = build_mbc_matrix(problem,time,rotor_omega)
        mbc_matrix_derivative = build_mbc_matrix_derivative(problem,time,rotor_omega)
        # complete the tranformation
        damping_matrix_result = np.transpose(mbc_matrix)@damping_matrix_rot(time,rotor_omega)@mbc_matrix + 2*np.transpose(mbc_matrix)@mass_matrix_rot(time,rotor_omega)@mbc_matrix_derivative
        return damping_matrix_result
    
    # the damping matrix is now ready, we move to the stiffness one    
    def stiffness_matrix_K(time:float ,rotor_omega:float):
        # 1) K without dampers
            #block generation of matrix
        diagonal = lag_hinge_offset*rotor_omega**2*blade_static_moment*np.ones(num_blades)
        upper_left = np.diag(diagonal)
            # Lower right
        lower_right = np.array([
                [hub_stiffness_x,0],
                [0,hub_stiffness_y]
                ])
        # Time varying terms (which use a psi azimuth angle time function)
        # the upper right term is a Nb*2 array containing zeros
        upper_right = np.zeros((num_blades,2))
        # the lower left term is a 2*Nb array with periodic terms
        upper_vec=np.zeros(num_blades)
        lower_vec=np.zeros(num_blades)
        for blade_index in range(1,num_blades+1):
            psi = azimuth_angle(blade_index,num_blades)
            upper_vec[blade_index-1]=-rotor_omega**2*blade_static_moment*np.sin(psi(rotor_omega,time))
            lower_vec[blade_index-1]=rotor_omega**2*blade_static_moment*np.cos(psi(rotor_omega,time))
        lower_left = np.vstack([upper_vec,lower_vec])
        #generation of mass matrix (partially in rotating frame)
        stiffness_matrix_partial = np.block([[upper_left,upper_right],[lower_left,lower_right]])
            
        # 2) K with the addition of the dampers
        stiffness_matrix_addition = np.zeros((num_blades+2,num_blades+2))
        if problem.damper_connection == "H2B":
            diagonal = nominal_stiffness*np.ones(num_blades)
            stiffness_matrix_addition[0:num_blades,0:num_blades] = np.diag(diagonal)
        #elif problem.damper_connection == "B2B":
                # B2B damper connection logic to be implemented
                
            # summation of matrices
        stiffness_matrix_result= stiffness_matrix_partial+stiffness_matrix_addition
            # 3) transformation through MBC
            # build the MBC matrices
        mbc_matrix = build_mbc_matrix(problem,time,rotor_omega)
        mbc_matrix_derivative = build_mbc_matrix_derivative(problem,time,rotor_omega)
        mbc_matrix_derivative_second = build_mbc_matrix_derivative_second(problem,time,rotor_omega)
            # complete the tranformation
        stiffness_matrix_result = np.transpose(mbc_matrix)@stiffness_matrix_result@mbc_matrix+ np.transpose(mbc_matrix)@damping_matrix_rot(time,rotor_omega)@mbc_matrix_derivative + np.transpose(mbc_matrix)@mass_matrix_rot(time,rotor_omega)@mbc_matrix_derivative_second
        return stiffness_matrix_result
        
    return damping_matrix_C, stiffness_matrix_K

def azimuth_angle(
        blade_index: int,
        num_blades: int
        ): 
    """ Function defining the azimth angle psi in the rotor.
    
    Args:
        blade_index: the current blade under consideration
        num_blades: the number of blades in the rotor
        
        
    Returns: 
        azimuth_angle: a function handle dependent on time, 
        rotor_Omega
    """
    azimuth_angle = lambda rotor_omega,time: rotor_omega*time+2*np.pi/num_blades*(blade_index-1)

    return azimuth_angle

def damping_activation_H2B(problem:ProblemDefinition):
    """ This function generates the diagonal term of the additional damping matrix 
    for the H2B problem
    
    Args: ProblemDefinition object
    
    Returns: damping_vector which becomes the diagonal of the additional damping matrix
    """

    num_blades = problem.number_blades
    blade_damping = problem.rotor_characteristics.nominal_damping_Cd
    if problem.damper_activation == "ALL":
        damping_vector = blade_damping*np.ones(num_blades)
    elif problem.damper_activation == "ODI":
        damping_vector = blade_damping*np.ones(num_blades)
        damping_vector[0] = 0
    elif problem.damper_activation == "CUSTOM":
        damping_vector = blade_damping*problem.damper_activation_vector

    return damping_vector

def build_mbc_matrix(problem:ProblemDefinition, time: float, rotor_omega:float)-> np.ndarray:
    """ This function build the mbc_matrix T as a function of time and omega.
    
    Args: ProblemDefinition
    
    Returns: mbc_matrix
    """
    num_blades = problem.number_blades

    num_cyclic = number_mbc(num_blades)

    # definition of the true MBC block
    left_col=np.ones(num_blades)
    right_col=np.zeros(num_blades) #used only when even Nb
    center_block =np.zeros((num_blades,2*num_cyclic))

    # filling the variable terms

    for blade_index in range(1,num_blades+1): #runs along the lines
        if num_blades % 2 == 0: #even blades, we need scissor states
            right_col[blade_index-1]=(-1)**blade_index
        for mbc_index in range(1, num_cyclic + 1):  # mbc_index = 1, 2, ..., num_cyclic
            psi = azimuth_angle(blade_index, num_blades)
            col_cos = 2 * (mbc_index - 1)      # columns 0, 2, 4, ...
            col_sin = 2 * (mbc_index - 1) + 1  # columns 1, 3, 5, ...
            center_block[blade_index - 1, col_cos] = np.cos(mbc_index * psi(rotor_omega, time))
            center_block[blade_index - 1, col_sin] = np.sin(mbc_index * psi(rotor_omega, time))
            #assemble the MBC block
    if num_blades % 2 == 0:
        mbc_matrix = np.column_stack([left_col,center_block,right_col])
    elif num_blades % 2 == 1:
        mbc_matrix = np.column_stack([left_col,center_block])
    
    # definition of the constant blocks
    right_cols = np.zeros((num_blades,2))
    left_rows = np.zeros((2,num_blades))
    right_block =np.eye(2)

    #building the final matrix
    mbc_matrix = np.block([[mbc_matrix,right_cols],[left_rows,right_block]])

    return mbc_matrix

def build_mbc_matrix_derivative(problem:ProblemDefinition, time: float, rotor_omega:float)-> np.ndarray:
    """ This function build the time derivative of mbc_matrix T as a function of time and omega.
    
    Args: ProblemDefinition
    
    Returns: mbc_matrix
    """
    num_blades = problem.number_blades

    num_cyclic = number_mbc(num_blades)

    # definition of the true MBC block
    left_col=np.zeros(num_blades)
    right_col=np.zeros(num_blades) #used only when even Nb
    center_block =np.zeros((num_blades,2*num_cyclic))

    # filling the variable terms

    for blade_index in range(1,num_blades+1): #runs along the lines
        for mbc_index in range(1, num_cyclic + 1):  # mbc_index = 1, 2, ..., num_cyclic
            psi = azimuth_angle(blade_index, num_blades)
            col_sin = 2 * (mbc_index - 1)      # columns 0, 2, 4, ... (derivative of cos)
            col_cos = 2 * (mbc_index - 1) + 1  # columns 1, 3, 5, ... (derivative of sin)
            center_block[blade_index - 1, col_sin] = -rotor_omega*mbc_index*np.sin(mbc_index * psi(rotor_omega, time))
            center_block[blade_index - 1, col_cos] = rotor_omega*mbc_index*np.cos(mbc_index * psi(rotor_omega, time))
            #assemble the MBC block
    if num_blades % 2 == 0:
        mbc_matrix = np.column_stack([left_col,center_block,right_col])
    elif num_blades % 2 == 1:
        mbc_matrix = np.column_stack([left_col,center_block])
    
    # definition of the constant blocks (all zeros)
    right_cols = np.zeros((num_blades,2))
    left_rows = np.zeros((2,num_blades))
    right_block =np.zeros((2,2))

    #building the final matrix
    mbc_matrix_derivative = np.block([[mbc_matrix,right_cols],[left_rows,right_block]])

    return mbc_matrix_derivative

def build_mbc_matrix_derivative_second(problem:ProblemDefinition, time: float, rotor_omega:float)-> np.ndarray:
    """ This function builds the seond time derivative of mbc_matrix T as a function of time and omega.
    
    Args: ProblemDefinition
    
    Returns: mbc_matrix
    """
    num_blades = problem.number_blades

    num_cyclic = number_mbc(num_blades)

    # definition of the true MBC block
    left_col=np.zeros(num_blades)
    right_col=np.zeros(num_blades) #used only when even Nb
    center_block =np.zeros((num_blades,2*num_cyclic))

    # filling the variable terms

    for blade_index in range(1,num_blades+1): #runs along the lines
        for mbc_index in range(1, num_cyclic + 1):  # mbc_index = 1, 2, ..., num_cyclic
            psi = azimuth_angle(blade_index, num_blades)
            col_cos = 2 * (mbc_index - 1)      # columns 0, 2, 4, ... (derivative of sin)
            col_sin = 2 * (mbc_index - 1) + 1  # columns 1, 3, 5, ... (derivative of cos)
            center_block[blade_index - 1, col_cos] = -(rotor_omega*mbc_index)**2*np.cos(mbc_index * psi(rotor_omega, time))
            center_block[blade_index - 1, col_sin] = -(rotor_omega*mbc_index)**2*np.sin(mbc_index * psi(rotor_omega, time))
            #assemble the MBC block
    if num_blades % 2 == 0:
        mbc_matrix = np.column_stack([left_col,center_block,right_col])
    elif num_blades % 2 == 1:
        mbc_matrix = np.column_stack([left_col,center_block])
    
    # definition of the constant blocks (all zeros)
    right_cols = np.zeros((num_blades,2))
    left_rows = np.zeros((2,num_blades))
    right_block =np.zeros((2,2))

    #building the final matrix
    mbc_matrix_derivative_second = np.block([[mbc_matrix,right_cols],[left_rows,right_block]])

    return mbc_matrix_derivative_second

def number_mbc(num_blades) -> int:
    """ Computation of the number of MBC cyclic coordinates (or N*).
    
    Args: num_blades

    Returns: number_cyclic

    """
    if num_blades % 2 == 0: # even number of blades
        number_cyclic = num_blades//2-1
    elif num_blades % 2 == 1: #odd number of blades
        number_cyclic = (num_blades-1)//2

    return number_cyclic

if __name__ == "__main__":
    # Test matrix construction
    from ..config.problem_definition import create_default_problem
    
    problem = create_default_problem()
    
    # Test mass matrix
    M = build_mass_matrix(problem)
    print(f"Mass matrix shape: {M.shape}")
    
    # Test damping and stiffness matrices
    C_func, K_func = build_damping_matrix(problem)
    C = C_func(0.0, 25.0)
    K = K_func(0.0, 25.0)
    print(f"Damping matrix shape: {C.shape}")
    print(f"Stiffness matrix shape: {K.shape}")
