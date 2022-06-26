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
Cp = 0.1
Cs = 0.1
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

AF = 1/(1-1.15*Cp-Cs)
print('\nBasic Amplification 1/(1-1.15Cp-Cs) = %.3f' % (AF))

print('==== Total Load ====')
print(f'No Ponding Effect:    {results0.total_factored_load:.3f}')
print(f'With Ponding Effect:  {resultsP.total_factored_load:.3f}')
print(f'Amplification Factor: {resultsP.total_factored_load/results0.total_factored_load:.3f}')

print('==== Member Output ====')
print('  Member         Max Shear      Max Moment         Max Shear      Max Moment        Shear AF         Moment AF')
print('                  (kips)         (kip-ft)           (kips)         (kip-ft)           (---)            (---)  ')
for i in range(bay.number_of_joist_spaces+1):
    Vmax0 = max(results0.secondary_members_bot_reaction[i],results0.secondary_members_top_reaction[i])
    Mmax0 = max(results0.secondary_members_moment[i,:])/12
    VmaxP = max(resultsP.secondary_members_bot_reaction[i],resultsP.secondary_members_top_reaction[i])
    MmaxP = max(resultsP.secondary_members_moment[i,:])/12
    print(f'Secondary {i+1:<2d}   {Vmax0:>8.2f}         {Mmax0:>8.2f}         {VmaxP:>8.2f}         {MmaxP:>8.2f}         {VmaxP/Vmax0:>8.3f}         {MmaxP/Mmax0:>8.3f}')
Vmax0 = max(abs(results0.top_primary_member_shear))
Mmax0 = max(results0.top_primary_member_moment)/12
VmaxP = max(abs(resultsP.top_primary_member_shear))
MmaxP = max(resultsP.top_primary_member_moment)/12
print(f'Primary Top    {Vmax0:>8.2f}         {Mmax0:>8.2f}         {VmaxP:>8.2f}         {MmaxP:>8.2f}         {VmaxP/Vmax0:>8.3f}         {MmaxP/Mmax0:>8.3f}')
Vmax0 = max(abs(results0.bot_primary_member_shear))
Mmax0 = max(results0.bot_primary_member_moment)/12
VmaxP = max(abs(resultsP.bot_primary_member_shear))
MmaxP = max(resultsP.bot_primary_member_moment)/12
print(f'Primary Bot    {Vmax0:>8.2f}         {Mmax0:>8.2f}         {VmaxP:>8.2f}         {MmaxP:>8.2f}         {VmaxP/Vmax0:>8.3f}         {MmaxP/Mmax0:>8.3f}')





