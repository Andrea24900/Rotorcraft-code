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

def build_MCK_matrices(
    problem: ProblemDefinition,mbc_matrices = None
) -> Tuple [np.ndarray,np.ndarray,np.ndarray]:
    """Build M,C,K matrices for the rotor system. If the MBC matrices are not created, build them.

    Args:
        problem: Problem definition containing rotor characteristics

    Returns:
        mass_matrix, damping_matrix, stiffness_matrix (function handles in Omega,t)

    Note:
        Matrix dimensions depend on number of blades (Nb blade states and 2 hub ones):
        - 3 blades: 5x5 matrix
        - 4 blades: 6x6 matrix
        - 5 blades: 7x7 matrix
        - ...
    """
    rotor_chars = problem.rotor_characteristics
    num_blades = rotor_chars.number_blades

    # Extract rotor characteristics for clarity
    blade_mass = rotor_chars.blade_mass_Mb
    blade_inertia = rotor_chars.blade_inertia_Ib
    hub_mass_x = rotor_chars.hub_mass_Mx
    hub_mass_y = rotor_chars.hub_mass_My
    blade_static_moment = rotor_chars.blade_static_Sb
    lag_hinge_offset = rotor_chars.lag_hinge_offset_e
    nominal_damping = rotor_chars.nominal_damping_Cd
    nominal_stiffness = rotor_chars.nominal_stiffness_Kd
    hub_damping_x = rotor_chars.hub_damping_Cx
    hub_damping_y = rotor_chars.hub_damping_Cy
    hub_stiffness_x = rotor_chars.hub_stiffness_Kx
    hub_stiffness_y = rotor_chars.hub_stiffness_Ky
    # Total hub mass including blades
    total_hub_mass_x = hub_mass_x + blade_mass * num_blades
    total_hub_mass_y = hub_mass_y + blade_mass * num_blades

   
    def MCK_matrices(time: float, rotor_omega:float) -> Tuple [np.ndarray,np.ndarray,np.ndarray]:
     
        # define mass_matrix in rotating frame
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
        
        # define damping matrix in rotating frame 
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
            if rotor_chars.damper_connection == "H2B":
                diagonal = damping_activation_H2B(problem)
                damping_matrix_addition[0:num_blades,0:num_blades] = np.diag(diagonal)
            elif rotor_chars.damper_connection == "B2B":
                damping_matrix_addition, _ = damping_activation_B2B(problem)
                
            # summation of matrices
            damping_matrix_result = damping_matrix_partial+damping_matrix_addition
            # 3) transformation through MBC
            return damping_matrix_result
        
        # definition of stiffness matrix in rotating frame 
        def stiffness_matrix_rot(time:float, rotor_omega:float):
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
            if rotor_chars.damper_connection == "H2B":
                diagonal = nominal_stiffness*np.ones(num_blades)
                stiffness_matrix_addition[0:num_blades,0:num_blades] = np.diag(diagonal)
            elif rotor_chars.damper_connection == "B2B":
                _, stiffness_matrix_addition = damping_activation_B2B(problem)
                    
                # summation of matrices
            stiffness_matrix_result= stiffness_matrix_partial+stiffness_matrix_addition
        
            return stiffness_matrix_result
        
        #generation of MBC matrices
        mbc_matrix = build_mbc_matrix(problem,time,rotor_omega)
        mbc_matrix_derivative = build_mbc_matrix_derivative(problem,time,rotor_omega)
        mbc_matrix_derivative_second = build_mbc_matrix_derivative_second(problem,time,rotor_omega)
        mbc_matrix_transpose = np.transpose(mbc_matrix)

        #evaluation of the system matrices

        mass_matrix_rot_eval = mass_matrix_rot(time, rotor_omega)
        damping_matrix_rot_eval = damping_matrix_rot(time,rotor_omega)
        stiffness_matrix_rot_eval = stiffness_matrix_rot(time,rotor_omega)

        # rotation of the matrices
        mass_matrix = mbc_matrix_transpose@mass_matrix_rot_eval@mbc_matrix
        damping_matrix = mbc_matrix_transpose@damping_matrix_rot_eval@mbc_matrix + 2*mbc_matrix_transpose@mass_matrix_rot_eval@mbc_matrix_derivative
        stiffness_matrix = mbc_matrix_transpose@stiffness_matrix_rot_eval@mbc_matrix+ mbc_matrix_transpose@damping_matrix_rot_eval@mbc_matrix_derivative + mbc_matrix_transpose@mass_matrix_rot_eval@mbc_matrix_derivative_second
        
        
        return mass_matrix,damping_matrix,stiffness_matrix
    
    return MCK_matrices

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

    rotor_chars = problem.rotor_characteristics
    num_blades = rotor_chars.number_blades
    blade_damping = rotor_chars.nominal_damping_Cd
    if rotor_chars.damper_activation == "ALL":
        damping_vector = blade_damping*np.ones(num_blades)
    elif rotor_chars.damper_activation == "ODI":
        damping_vector = blade_damping*np.ones(num_blades)
        damping_vector[0] = 0
    elif rotor_chars.damper_activation == "CUSTOM":
        damping_vector = blade_damping*rotor_chars.damper_activation_vector

    return damping_vector

def damping_activation_B2B(problem: ProblemDefinition):
    """Generate the damping and stiffness addition matrices for B2B (blade-to-blade) dampers.

    In B2B configuration, each damper connects blade i to blade i+1 (cyclically).
    This creates a matrix with:
    - Diagonal terms: contribution from dampers on both sides of each blade
    - Extra-diagonal terms: coupling between adjacent blades

    Args:
        problem: ProblemDefinition object containing rotor characteristics

    Returns:
        Tuple of (damping_matrix_addition, stiffness_matrix_addition)
        Both are (num_blades+2) x (num_blades+2) matrices
    """
    rotor_chars = problem.rotor_characteristics
    num_blades = rotor_chars.number_blades

    # Extract damper parameters
    cd = rotor_chars.nominal_damping_Cd
    kd = rotor_chars.nominal_stiffness_Kd
    a = rotor_chars.hub_to_A_a
    b = rotor_chars.hub_to_B_b
    zeta_E = rotor_chars.zeta_E
    phi = rotor_chars.phi
    gamma_E = rotor_chars.gamma_E
    lE = rotor_chars.length_eq_lE
    l0 = rotor_chars.length_undef_l0

    # Precompute trigonometric terms
    sin_minus = np.sin(zeta_E - phi + gamma_E)  # sin(zeta_E - phi + gamma_E)
    sin_plus = np.sin(zeta_E + phi + gamma_E)   # sin(zeta_E + phi + gamma_E)
    cos_minus = np.cos(zeta_E - phi + gamma_E)  # cos(zeta_E - phi + gamma_E)
    cos_plus = np.cos(zeta_E + phi + gamma_E)   # cos(zeta_E + phi + gamma_E)

    # B2B damping coefficients (from symbolic derivation)
    # C_zetaed: extra-diagonal coupling term between adjacent blades
    C_zetaed = a * b * cd * sin_minus * sin_plus
    # C_before: contribution from damper before blade i (damper i-1, connects blade i-1 to i)
    C_before = b**2 * cd * sin_plus**2
    # C_after: contribution from damper after blade i (damper i, connects blade i to i+1)
    C_after = a**2 * cd * sin_minus**2
    # Note: C_zetad = C_before + C_after (full diagonal when both dampers active)

    # B2B stiffness coefficients (from symbolic derivation)
    # K_zetaed: extra-diagonal coupling term between adjacent blades
    K_zetaed = (a * b * kd * sin_minus * sin_plus
                + a * b * kd * (lE - l0) / l0 * cos_minus * cos_plus)
    # K_before: contribution from damper before blade i (damper i-1, connects blade i-1 to i)
    K_before = (b**2 * kd * sin_plus
                - b * kd * (lE - l0) / lE * cos_plus * (lE - b * cos_plus))
    # K_after: contribution from damper after blade i (damper i, connects blade i to i+1)
    K_after = (a**2 * kd * sin_minus
               - a * kd * (lE - l0) / lE * cos_minus * (lE - a * cos_minus))
    # Note: K_zetad = K_before + K_after (full diagonal when both dampers active)

    # Initialize matrices (only blade DOFs, hub DOFs remain zero)
    damping_matrix_addition = np.zeros((num_blades + 2, num_blades + 2))
    stiffness_matrix_addition = np.zeros((num_blades + 2, num_blades + 2))

    # Build activation vector based on damper_activation mode
    # damper_active[i] = 1 if damper i (connecting blade i to blade i+1) is active
    if rotor_chars.damper_activation == "ALL":
        damper_active = np.ones(num_blades)
    elif rotor_chars.damper_activation == "ODI":
        # Damper 1 (between blade 1 and blade 2, index 0) is inoperative
        damper_active = np.ones(num_blades)
        damper_active[0] = 0
    elif rotor_chars.damper_activation == "CUSTOM":
        damper_active = rotor_chars.damper_activation_vector
    else:
        damper_active = np.ones(num_blades)

    # Build the damping matrix for blade DOFs
    # Each blade i is affected by damper i-1 (connecting blade i-1 to blade i)
    # and damper i (connecting blade i to blade i+1)
    for i in range(num_blades):
        # Damper indices (cyclic): damper connecting blade i to blade i+1 has index i
        damper_prev = (i - 1) % num_blades  # Damper connecting blade i-1 to blade i
        damper_curr = i  # Damper connecting blade i to blade i+1

        # Diagonal term: sum of contributions from both adjacent dampers
        # Scale by damper activation ratio (0.0 to 1.0)
        diag_contrib = (damper_active[damper_prev] * C_before +
                        damper_active[damper_curr] * C_after)
        damping_matrix_addition[i, i] = diag_contrib

        # Extra-diagonal terms: coupling with adjacent blades
        i_next = (i + 1) % num_blades
        i_prev = (i - 1) % num_blades

        # Coupling with next blade (through damper i)
        damping_matrix_addition[i, i_next] = damper_active[damper_curr] * C_zetaed

        # Coupling with previous blade (through damper i-1)
        damping_matrix_addition[i, i_prev] = damper_active[damper_prev] * C_zetaed

    # Build the stiffness matrix for blade DOFs (same structure as damping)
    for i in range(num_blades):
        damper_prev = (i - 1) % num_blades
        damper_curr = i

        # Diagonal term: sum of contributions from both adjacent dampers
        # Scale by damper activation ratio (0.0 to 1.0)
        diag_contrib = (damper_active[damper_prev] * K_before +
                        damper_active[damper_curr] * K_after)
        stiffness_matrix_addition[i, i] = diag_contrib

        # Extra-diagonal terms
        i_next = (i + 1) % num_blades
        i_prev = (i - 1) % num_blades

        stiffness_matrix_addition[i, i_next] = damper_active[damper_curr] * K_zetaed
        stiffness_matrix_addition[i, i_prev] = damper_active[damper_prev] * K_zetaed

    return damping_matrix_addition, stiffness_matrix_addition

def build_mbc_matrix(problem:ProblemDefinition, time: float, rotor_omega:float)-> np.ndarray:
    """ This function build the mbc_matrix T as a function of time and omega.
    
    Args: ProblemDefinition
    
    Returns: mbc_matrix
    """
    num_blades = problem.rotor_characteristics.number_blades

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
    num_blades = problem.rotor_characteristics.number_blades

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
    num_blades = problem.rotor_characteristics.number_blades

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
    MCK_func = build_MCK_matrices(problem)
    M,C,K = MCK_func(time=0.0, rotor_omega=1.0)
    print(f"Mass matrix: {M}")
    print(f"Damping matrix: {C}")
    print(f"Stiffness matrix: {K}")
