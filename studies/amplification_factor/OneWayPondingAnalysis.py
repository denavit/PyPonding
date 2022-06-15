from PyPonding.structures import ElasticBeam2d
from math import pi, cos, cosh
import numpy as np
import matplotlib.pyplot as plt

# Define units
inch    = 1.0
kip     = 1.0
lb      = kip/1000.0
ft      = 12.0*inch
ksi     = kip/inch**2
psf     = lb/ft**2
pcf     = psf/ft
in_per_ft = inch/ft


# Input parameters
L       = 40*ft         # Beam span
S       = 5*ft          # Tributary width
E       = 29000.0*ksi   # Modulus of elasticity
I       = 75*inch**4    # Moment of inertia
gamma   = 62.4*pcf      # Unit weight of water
qD      = 0*psf         # Dead load
slope   = 0*in_per_ft # Beam slope
zw      = 2*inch        # Elevation of water


# Build beam object
beam = ElasticBeam2d(L,S,E,I,gamma,zj=slope*L,qD=qD)
beam.num_elements = 60

# Run Ponding Analysis
resultsP = beam.run_analysis_OPS('IterativeLevel',target_zw=zw)
VmaxP = np.amax(np.absolute(resultsP.shear_along_length))
MmaxP = np.amax(resultsP.bending_moment_along_length)
x_at_MmaxP = resultsP.position_along_length[np.argmax(resultsP.bending_moment_along_length)]
TotalLoadP = resultsP.shear_along_length[0]-resultsP.shear_along_length[-1]

beam.plot_deformed(zw=zw)

# Run First-Order Analysis
beam.include_ponding_effect = False
results1 = beam.run_analysis_OPS('IterativeLevel',target_zw=zw)
Vmax1 = np.amax(np.absolute(results1.shear_along_length))
Mmax1 = np.amax(results1.bending_moment_along_length)
x_at_Mmax1 = results1.position_along_length[np.argmax(results1.bending_moment_along_length)]
TotalLoad1 = results1.shear_along_length[0]-results1.shear_along_length[-1]


# Print Results
Cs = (gamma*S*L**4)/(pi**4*E*I)
print('\nCs = %.3f' % Cs)
print('Basic Amplification 1/(1-Cs) = %.3f' % (1/(1-Cs)))
x = pi/2*Cs**0.25
print('More Exact Amplification     = %.3f\n' % ((1/cos(x)-1/cosh(x))/x**2)) # from Silver (2010)
print('                      with ponding           no ponding          amplification')
print('Total Load          %9.3f kip         %9.3f kip          %9.3f' % (TotalLoadP,TotalLoad1,TotalLoadP/TotalLoad1))
print('Max Shear           %9.3f kip         %9.3f kip          %9.3f' % (VmaxP,Vmax1,VmaxP/Vmax1))
print('Max Moment          %9.3f kip-in      %9.3f kip-in       %9.3f' % (MmaxP,Mmax1,MmaxP/Mmax1))
print('Pos. of Mmax      %5.1f in (%4.3f*L)    %5.1f in (%4.3f*L)\n' % (x_at_MmaxP,x_at_MmaxP/L,x_at_Mmax1,x_at_Mmax1/L))


# Plot Results
plt.figure()
plt.plot(resultsP.position_along_length,resultsP.bending_moment_along_length,label='With Ponding')
plt.plot(results1.position_along_length,results1.bending_moment_along_length,label='No Ponding')
plt.legend()
plt.xlabel('Position along length of beam (in.)')
plt.ylabel('Bending Moment (kip-in.)')
plt.xlim([0,L])

plt.figure()
plt.plot(resultsP.position_along_length,resultsP.shear_along_length,label='With Ponding')
plt.plot(results1.position_along_length,results1.shear_along_length,label='No Ponding')
plt.legend()
plt.xlabel('Position along length of beam (in.)')
plt.ylabel('Shear (kips)')
plt.xlim([0,L])

plt.figure()
plt.plot(resultsP.position_along_length,resultsP.bending_moment_along_length/results1.bending_moment_along_length,label='Bending Moment')
plt.plot(results1.position_along_length,resultsP.shear_along_length/results1.shear_along_length,label='Shear')
plt.legend()
plt.xlabel('Position along length of beam (in.)')
plt.ylabel('Amplification Factor')
plt.xlim([0,L])
#plt.ylim([0,2*(1/(1-Cs))])

plt.show()