from PyPonding.structures import bay
from PyPonding import FE
from math import pi, cos, cosh
import numpy as np

# Define units
inch    = 1.0
kip     = 1.0
lb      = kip/1000.0
ft      = 12.0*inch
ksi     = kip/inch**2
psf     = lb/ft**2
pcf     = psf/ft

water_level = 2*inch

# Define Steel Beam Object
myBay = bay();

# Geometry
myBay.primary_member_span         = 40*ft
myBay.secondary_member_span       = 40*ft
myBay.num_spaces                  = 8

# Loads and load factors
myBay.dead_load_uniform           = 10.0*psf
myBay.dead_load_primary_member    = 0.0*lb/ft           # Self-weight of joist girder
myBay.water_density               = 62.4*pcf   

# Top of roof elevation
myBay.z_TL                        = 0.0*inch
myBay.z_TR                        = 0.0*inch
myBay.z_BL                        = 0.0*inch
myBay.z_BR                        = 0.0*inch

# Camber
myBay.primary_member_camber_T     = 0.0*inch
myBay.primary_member_camber_B     = 0.0*inch
myBay.secondary_member_camber     = 0.0*inch

# Edge conditions ('mirrored', 'rigid', or '')
myBay.edge_condition_L            = 'mirrored'
myBay.edge_condition_R            = 'mirrored'
myBay.edge_condition_T            = 'mirrored'
myBay.edge_condition_B            = 'mirrored'

# Material and section properties
myBay.E                           = 29000.0*ksi
myBay.Ap                          = 100*inch**2
myBay.As                          = 100*inch**2
myBay.Ip                          = 1820*inch**4
myBay.Is                          = 182*inch**4

# Analysis options
myBay.MAX_ITER                    = 50                # Maximum number of ponding analysis iterations

myBay.alpha                       = 1
myBay.load_factor_dead            = 1.0               # Dead
myBay.load_factor_ponding         = 1.0               # Impounded Water
myBay.load_factor_snow1           = 0.0               # Snow in Ponding Load Cell
myBay.load_factor_snow2           = 0.0               # Snow as Simple Load

myBay.snow_density                = 0.0
myBay.snow_height                 = 0.0




# Run Ponding Analysis
myBay.include_ponding_effect = True
print('Running Ponding Analysis')
myBay.Run_Analysis(water_level)
Mmax_TP = np.amax(myBay.results_top_primary_M)
Mmax_BP = np.amax(myBay.results_bot_primary_M)
Mmax_J = np.zeros(myBay.num_spaces+1)
for i in range(myBay.num_spaces+1):
    Mmax_J[i] = np.amax(myBay.results_secondary_M[i,:])
Vmax_TP = np.amax(myBay.results_top_primary_V)
Vmax_BP = np.amax(myBay.results_bot_primary_V)
Vmax_J = np.zeros(myBay.num_spaces+1)
for i in range(myBay.num_spaces+1):
    Vmax_J[i] = max(myBay.results_secondary_R_top[i],myBay.results_secondary_R_bot[i])
    
# Run First-Order Analysis
myBay.include_ponding_effect = False
myBay.Run_Analysis(water_level)
Mmaxo_TP = np.amax(myBay.results_top_primary_M)
Mmaxo_BP = np.amax(myBay.results_bot_primary_M)
Mmaxo_J = np.zeros(myBay.num_spaces+1)
for i in range(myBay.num_spaces+1):
    Mmaxo_J[i] = np.amax(myBay.results_secondary_M[i,:])
Vmaxo_TP = np.amax(myBay.results_top_primary_V)
Vmaxo_BP = np.amax(myBay.results_bot_primary_V)
Vmaxo_J = np.zeros(myBay.num_spaces+1)
for i in range(myBay.num_spaces+1):
    Vmaxo_J[i] = max(myBay.results_secondary_R_top[i],myBay.results_secondary_R_bot[i])
    
    
    
gamma = myBay.alpha*myBay.load_factor_ponding*myBay.water_density
Cp = gamma*myBay.secondary_member_span*myBay.primary_member_span**4/(pi**4*myBay.E*myBay.Ip)
Cs = gamma*(myBay.primary_member_span/myBay.num_spaces)*myBay.secondary_member_span**4/(pi**4*myBay.E*myBay.Is)
AF = 1/(1-1.15*Cp-Cs)

# Print Results
print('\nPrinting Results')
print('Cp = %.3f' % Cp)
print('Cs = %.3f' % Cs)
print('Basic Amplification 1/(1-1.15Cp-Cs) = %.3f' % (AF))
print('                                        with ponding           no ponding          amplification')
print('Max Moment in Top Primary Member      %9.3f kip-in      %9.3f kip-in       %9.3f' % (Mmax_TP,Mmaxo_TP,Mmax_TP/Mmaxo_TP))
print('Max Moment in Bot Primary Member      %9.3f kip-in      %9.3f kip-in       %9.3f' % (Mmax_BP,Mmaxo_BP,Mmax_BP/Mmaxo_BP))
maxAFm = max(Mmax_TP/Mmaxo_TP,Mmax_BP/Mmaxo_BP)
for i in range(myBay.num_spaces+1):
    print('Max Moment in Secondary Member %2i     %9.3f kip-in      %9.3f kip-in       %9.3f' % (i+1,Mmax_J[i],Mmaxo_J[i],Mmax_J[i]/Mmaxo_J[i]))
    maxAFm = max(maxAFm,Mmax_J[i]/Mmaxo_J[i])
print('Max Shear in Top Primary Member       %9.3f kip         %9.3f kip          %9.3f' % (Vmax_TP,Vmaxo_TP,Vmax_TP/Vmaxo_TP))
print('Max Shear in Bot Primary Member       %9.3f kip         %9.3f kip          %9.3f' % (Vmax_BP,Vmaxo_BP,Vmax_BP/Vmaxo_BP))
maxAFv = max(Vmax_TP/Vmaxo_TP,Vmax_BP/Vmaxo_BP)
for i in range(myBay.num_spaces+1):
    print('Max Shear in Secondary Member %2i      %9.3f kip         %9.3f kin          %9.3f' % (i+1,Vmax_J[i],Vmaxo_J[i],Vmax_J[i]/Vmaxo_J[i]))
    maxAFv = max(maxAFv,Vmax_J[i]/Vmaxo_J[i])
    
print('Maximum Moment Amplification = %9.3f' % (maxAFm))
print('Maximum Shear Amplification  = %9.3f' % (maxAFv))
print('Maximum Amplification        = %9.3f' % (max(maxAFm,maxAFv)))




