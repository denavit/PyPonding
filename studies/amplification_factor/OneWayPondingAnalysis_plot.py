import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos, cosh
from PyPonding.structures import ElasticBeam2d

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
C       = 0.25          # Flexibility coefficient (ponding factor)
gamma   = 62.4*pcf      # Unit weight of water
qD      = 0*psf         # Dead load
slope   = 0.25*in_per_ft   # Beam slope
zw      = 2*inch        # Elevation of water


# Define loading cases
zw_over_zj_max = 1.6
n = 49 # number of analyses

zj = slope*L
zw = np.linspace(0,zw_over_zj_max*zj,n)
Bp_V = np.zeros(n)
Bp_M = np.zeros(n)
Bp_TL = np.zeros(n)
Bpo = 1/(1-C)

for i in range(n):
    # Build beam object
    I = (gamma*S*L**4)/(pi**4*E*C) # Moment of inertia
    beam = ElasticBeam2d(L,S,E,I,gamma,zj=zj,qD=qD)
    beam.num_elements = 60

    # Run Ponding Analysis
    resultsP = beam.run_analysis_OPS('IterativeLevel',target_zw=zw[i])
    VmaxP = np.amax(np.absolute(resultsP.shear_along_length))
    MmaxP = np.amax(resultsP.bending_moment_along_length)
    x_at_MmaxP = resultsP.position_along_length[np.argmax(resultsP.bending_moment_along_length)]
    TotalLoadP = resultsP.shear_along_length[0]-resultsP.shear_along_length[-1]    
    
    # Run First-Order Analysis
    beam.include_ponding_effect = False
    results1 = beam.run_analysis_OPS('IterativeLevel',target_zw=zw[i])
    Vmax1 = np.amax(np.absolute(results1.shear_along_length))
    Mmax1 = np.amax(results1.bending_moment_along_length)
    x_at_Mmax1 = results1.position_along_length[np.argmax(results1.bending_moment_along_length)]
    TotalLoad1 = results1.shear_along_length[0]-results1.shear_along_length[-1]
    
    Bp_V[i] = VmaxP/Vmax1
    Bp_M[i] = MmaxP/Mmax1
    Bp_TL[i] = TotalLoadP/TotalLoad1
   

# Make Plot
fig = plt.figure()
ax = fig.add_axes([0.15,0.15,0.80,0.80])
line0, = plt.plot([0.0,0.2,0.8,zw_over_zj_max],[1,1,Bpo,Bpo],'--', label='Idealized')
line1, = plt.plot(zw/zj,Bp_V,'x-',label='Analysis (Shear)')
line2, = plt.plot(zw/zj,Bp_M,'x-',label='Analysis (Moment)')
line3, = plt.plot(zw/zj,Bp_TL,'x-',label='Analysis (Total Load)')
plt.xlabel('zw/zj')
plt.ylabel('Bp')
plt.xlim(0.0,zw_over_zj_max)
plt.legend(handles=[line0,line1,line2,line3],frameon=False)
plt.show()
