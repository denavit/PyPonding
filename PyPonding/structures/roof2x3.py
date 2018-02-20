from PyPonding.structures import basic_structure

class roof2x3(basic_structure.basic_structure):

    # Bay Widths
    L_AB = 480
    L_BC = 480
    L_CD = 480
    L_12 = 480
    L_23 = 480

    # Roof Elevations at Columns
    z_A1 =   0
    z_A2 =  -8
    z_A3 =   0
    z_B1 =   0
    z_B2 = -10
    z_B3 =   0
    z_C1 =   0
    z_C2 = -10
    z_C3 =   0
    z_D1 =   0
    z_D2 =  -8
    z_D3 =   0

    # Joist and Joist Girder Section Properties
    E = 29000

    A_J   = 100
    Iz_J  = 215.0
    Iy_J  = 100
    GJ_J  = 100

    A_JG  = 100
    Iz_JG = 2029
    Iy_JG = 100
    GJ_JG = 100

    # Strength properties 
    w_J   = 253/1000/12    
    P_JG  = 7.5
    
    # Other Properties
    c_J     = 1
    c_JG    = 1
    nspaces = 8
    
    # Analysis Options
    nele_J  = 20
    dof_types = ('UX','UY','UZ','RX','RY','RZ')
    na      = 4
    nb      = 4    
    
    def __init__(self):
        pass
        
    def lowest_point(self):
        return min([ \
           self.z_A1, self.z_A2, self.z_A3,
           self.z_B1, self.z_B2, self.z_B3,
           self.z_C1, self.z_C2, self.z_C3,
           self.z_D1, self.z_D2, self.z_D3])         