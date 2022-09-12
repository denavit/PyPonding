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


def run_single_analysis(case,Cp,Cs,roof_slope,zw_over_zj,qD):

    # Bay Geometry
    Lp = 40*ft
    Ls = 40*ft
    num_spaces = 16
    S  = Lp/num_spaces

    # Loads
    γw = 62.4*pcf
    set_primary_member_load_zero = False

    # Material properties
    E = 29000.0*ksi

    # Member stiffness
    Ip = (γw*Ls*Lp**4)/(pi**4*E*Cp)
    Is = (γw*S*Ls**4)/(pi**4*E*Cs)

    if case == 'A':
        
        # Roof Slope
        zj = roof_slope*Lp
        zw = zw_over_zj*zj
        
        # Primary Member Load
        wD_T = ((384*E*Ip)/(5*360*Lp**3))*plf  
        wD_B = 0
        
        z_TL = roof_slope*Ls
        z_TR = roof_slope*Ls
        z_BL = 0.0*inch
        z_BR = 0.0*inch
        
        edge_condition_L = 'mirrored'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'none'
        edge_condition_B = 'rigid'
        
    elif case == 'B':
        
        # Roof Slope
        zj = roof_slope*Lp
        zw = zw_over_zj*zj
        
        # Primary Member Load
        wD_T = ((384*E*Ip)/(5*360*Lp**3))*plf   
        wD_B = 0
        
        z_TL = roof_slope*Ls
        z_TR = roof_slope*Ls
        z_BL = 0.0*inch
        z_BR = 0.0*inch
        
        edge_condition_L = 'mirrored'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'none'
        edge_condition_B = 'mirrored'
        
    elif case == 'C':
        
        # Roof Slope
        zj = roof_slope*Ls
        zw = zw_over_zj*zj
        
        # Primary Member Load
        wD_T = 0
        wD_B = 0
        
        z_TL = 0.0*inch
        z_TR = roof_slope*Lp
        z_BL = 0.0*inch
        z_BR = roof_slope*Lp
        
        edge_condition_L = 'rigid'
        edge_condition_R = 'none'
        edge_condition_T = 'mirrored'
        edge_condition_B = 'mirrored'
        
    elif case == 'D':
        
        # Roof Slope
        zj = roof_slope*Ls
        zw = zw_over_zj*zj
        
        # Primary Member Load
        wD_T = 0
        wD_B = 0
        
        z_TL = 0.0*inch
        z_TR = roof_slope*Lp
        z_BL = 0.0*inch
        z_BR = roof_slope*Lp
        
        edge_condition_L = 'mirrored'
        edge_condition_R = 'none'
        edge_condition_T = 'mirrored'
        edge_condition_B = 'mirrored'
        
    elif case == 'E':
        
        # Roof Slope
        zj = roof_slope*Lp
        zw = zw_over_zj*zj
        
        # Primary Member Load
        wD_T = 0
        wD_B = 0
        
        z_TL = 0.0*inch
        z_TR = roof_slope*Ls
        z_BL = 0.0*inch
        z_BR = 0.0*inch
        
        edge_condition_L = 'rigid'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'mirrored'
        edge_condition_B = 'rigid'
        
    elif case == 'F':
        
        # Roof Slope
        zj = roof_slope*Lp
        zw = zw_over_zj*zj
        
        # Primary Member Load
        wD_T = 0
        wD_B = 0
        
        z_TL = roof_slope*Ls
        z_TR = 0.0*inch
        z_BL = roof_slope*Ls
        z_BR = roof_slope*Ls
        
        edge_condition_L = 'none'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'mirrored'
        edge_condition_B = 'none'
        
    else:
        raise ValueError(f'Unknown case {case}')

    if set_primary_member_load_zero == True:
        wD_T = 0
        wD_B = 0

    # Define IdealizedBay object
    bay_input = {
        'primary_member_span': Lp,
        'secondary_member_span': Ls,
        'number_of_joist_spaces': num_spaces,
        'dead_load_uniform': qD,
        'dead_load_on_top_primary_member': wD_T,
        'dead_load_on_bottom_primary_member': wD_B,
        'water_density': γw,  
        'snow_density': 0.0*pcf,
        'snow_height': 0.0*inch,
        'alpha': 1.0,
        'load_factor_dead':    1.0,
        'load_factor_ponding': 1.0,
        'load_factor_snow':    1.0,
        'z_TL': z_TL,
        'z_TR': z_TR,
        'z_BL': z_BL,
        'z_BR': z_BR,
        'secondary_member_camber': 0.000*inch,
        'primary_member_camber_T': 0.000*inch,
        'primary_member_camber_B': 0.000*inch,
        'edge_condition_L': edge_condition_L,
        'edge_condition_R': edge_condition_R,
        'edge_condition_T': edge_condition_T,
        'edge_condition_B': edge_condition_B,
        'E':  E,
        'As': 100*inch**2,
        'Ap': 100*inch**2,
        'Is': Is,
        'Ip': Ip,
        'analsis_engine': 'FE',
    }
    bay = IdealizedBay(**bay_input)

    bay.include_ponding_effect = False
    results0 = bay.Run_Analysis(zw)
    bay.include_ponding_effect = True
    resultsP = bay.Run_Analysis(zw)

    amplification_factors = dict()

    AF = 1/(1-1.15*Cp-Cs)
    print('\nBasic Amplification 1/(1-1.15Cp-Cs) = %.3f' % (AF))
    amplification_factors['basic'] = AF
    
    print('==== Total Load ====')
    print(f'No Ponding Effect:    {results0.total_factored_load:.3f}')
    print(f'With Ponding Effect:  {resultsP.total_factored_load:.3f}')
    print(f'Amplification Factor: {resultsP.total_factored_load/results0.total_factored_load:.3f}')
    amplification_factors['total_load'] = resultsP.total_factored_load/results0.total_factored_load

    print('==== Member Output ====')
    print('  Member         Max Shear      Max Moment         Max Shear      Max Moment        Shear AF         Moment AF')
    print('                  (kips)         (kip-ft)           (kips)         (kip-ft)           (---)            (---)  ')
    
    # Secondary Members
    secondary_amplification_factors = []
    for i in range(bay.number_of_joist_spaces+1):
        if i == 0 and bay.edge_condition_L == 'rigid':
            secondary_amplification_factors.append(None)
            continue
        if i == bay.number_of_joist_spaces and bay.edge_condition_R == 'rigid':
            secondary_amplification_factors.append(None)
            continue  
        Vmax0 = max(results0.secondary_members_bot_reaction[i],results0.secondary_members_top_reaction[i])
        Mmax0 = max(results0.secondary_members_moment[i,:])/12
        VmaxP = max(resultsP.secondary_members_bot_reaction[i],resultsP.secondary_members_top_reaction[i])
        MmaxP = max(resultsP.secondary_members_moment[i,:])/12
        if Vmax0 == 0.0:
            if VmaxP == 0.0:
                V_AF_str = '  Undef.'
            else:
                V_AF_str = '   Inf  '
        else:
            V_AF_str = f'{VmaxP/Vmax0:>8.3f}'
        if Mmax0 == 0.0:
            if MmaxP == 0.0:
                M_AF_str = '  Undef.'
                secondary_amplification_factors.append('Undefined')
            else:
                M_AF_str = '   Inf  '
                secondary_amplification_factors.append('Infinite')
        else:
            M_AF_str = f'{MmaxP/Mmax0:>8.3f}'
            secondary_amplification_factors.append(max(VmaxP/Vmax0,MmaxP/Mmax0))
        
        print(f'Secondary {i+1:<2d}   {Vmax0:>8.2f}         {Mmax0:>8.2f}         {VmaxP:>8.2f}         {MmaxP:>8.2f}         {V_AF_str}         {M_AF_str}')
    
    amplification_factors['secondary'] = secondary_amplification_factors
        
        
    # Top Primary Member
    if bay.edge_condition_T == 'rigid':
        amplification_factors['top'] = None
    else:
        Vmax0 = max(abs(results0.top_primary_member_shear))
        Mmax0 = max(results0.top_primary_member_moment)/12
        VmaxP = max(abs(resultsP.top_primary_member_shear))
        MmaxP = max(resultsP.top_primary_member_moment)/12
        print(f'Primary Top    {Vmax0:>8.2f}         {Mmax0:>8.2f}         {VmaxP:>8.2f}         {MmaxP:>8.2f}         {VmaxP/Vmax0:>8.3f}         {MmaxP/Mmax0:>8.3f}')
        amplification_factors['top'] = max(VmaxP/Vmax0,MmaxP/Mmax0)
        
    # Bottom Primary Member
    if bay.edge_condition_B == 'rigid':
        amplification_factors['bottom'] = None
    else:
        Vmax0 = max(abs(results0.bot_primary_member_shear))
        Mmax0 = max(results0.bot_primary_member_moment)/12
        VmaxP = max(abs(resultsP.bot_primary_member_shear))
        MmaxP = max(resultsP.bot_primary_member_moment)/12
        print(f'Primary Bottom {Vmax0:>8.2f}         {Mmax0:>8.2f}         {VmaxP:>8.2f}         {MmaxP:>8.2f}         {VmaxP/Vmax0:>8.3f}         {MmaxP/Mmax0:>8.3f}')
        amplification_factors['bottom'] = max(VmaxP/Vmax0,MmaxP/Mmax0)
  
    return amplification_factors

    
if __name__ == "__main__":
    case = 'C'
    Cp = 0.3
    Cs = 0.3
    roof_slope = 0.5*in_per_ft
    zw_over_zj = 0.5
    qD = 0.0*psf
    amplification_factors = run_single_analysis(case,Cp,Cs,roof_slope,zw_over_zj,qD)
    print(amplification_factors)
    