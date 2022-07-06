from PyPonding.structures import IdealizedBay
from math import pi, cos, cosh
import numpy as np

# Define units
inch    = 1.0
kip     = 1.0
lb      = kip/1000.0
ft      = 12.0*inch
in_per_ft = inch/ft
ksi     = kip/inch**2
plf     = lb/ft
psf     = lb/ft**2
pcf     = psf/ft
kipft   = kip*ft


Cp_list = np.linspace(0.01,0.40,4)
Cs_list = np.linspace(0.01,0.40,4)
total_load_amplification_factor = np.zeros((4,4))


for iCp, Cp in enumerate(Cp_list):
    for iCs, Cs in enumerate(Cs_list):
        print(f'Running with {Cp = } and {Cs = }')
    
        # Bay Geometry
        Lp = 40*ft
        Ls = 40*ft
        num_spaces = 8
        S  = Lp/num_spaces

        # Loads
        qD = 10.0*psf
        γw = 62.4*pcf   
        zw = 2*inch

        # Top of roof elevation
        zTL = 0.0*inch
        zTR = 0.0*inch
        zBL = 0.0*inch
        zBR = 0.0*inch

        # Material properties
        E = 29000.0*ksi

        # Member stiffness
        Ip = (γw*Ls*Lp**4)/(pi**4*E*Cp)
        Is = (γw*S*Ls**4)/(pi**4*E*Cs)

        # Define IdealizedBay object
        input = {
            'primary_member_span': Lp,
            'secondary_member_span': Ls,
            'number_of_joist_spaces': num_spaces,
            'dead_load_uniform': qD,
            'dead_load_primary_member': 0*plf,
            'water_density': γw,  
            'snow_density': 0.0*pcf,
            'snow_height': 0.0*inch,
            'alpha': 1.0,
            'load_factor_dead':    1.0,
            'load_factor_ponding': 1.0,
            'load_factor_snow':    1.0,
            'z_TL': zTL,
            'z_TR': zTR,
            'z_BL': zBL,
            'z_BR': zBR,
            'secondary_member_camber': 0.000*inch,
            'primary_member_camber_T': 0.000*inch,
            'primary_member_camber_B': 0.000*inch,
            'edge_condition_L': 'mirrored',
            'edge_condition_R': 'mirrored',
            'edge_condition_T': 'mirrored',
            'edge_condition_B': 'mirrored',
            'E':  E,
            'As': 100*inch**2,
            'Ap': 100*inch**2,
            'Is': Is,
            'Ip': Ip,
            'analsis_engine': 'FE',
        }
        bay = IdealizedBay(**input)


        bay.include_ponding_effect = False
        results0 = bay.Run_Analysis(zw)
        bay.include_ponding_effect = True
        resultsP = bay.Run_Analysis(zw)

        # Store Results
        total_load_amplification_factor[iCp,iCs] = resultsP.total_factored_load/results0.total_factored_load


# Save Results to File
np.savetxt('results_Cp.csv', Cp_list, delimiter=',')
np.savetxt('results_Cs.csv', Cs_list, delimiter=',')
np.savetxt('results_total_load_amplification_factor.csv', total_load_amplification_factor, delimiter=',')


