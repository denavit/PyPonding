import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos, cosh
from PyPonding.structures.steel_beam import steel_beam
from PyPonding import FE

# Define units
inch    = 1.0
kip     = 1.0
lb      = kip/1000.0
ft      = 12.0*inch
ksi     = kip/inch**2
psf     = lb/ft**2
pcf     = psf/ft

# Input parameters
E       = 29000.0*ksi   # Modulus of elasticity
I       = 50*inch**4   # Moment of inertia
L       = 30*ft         # Beam span
TW      = 5*ft          # Tributary width
c       = 0*inch        # Camber
qD      = 10.0*psf      # Dead load (force per area)
gamma   = 62.4*pcf      # Unit weight of water
slope   = 1.0*inch/ft
zi      = 0*inch        # Elevation of left end
zj      = slope*L       # Elevation of right end

zw_over_zj_max = 1.6
n = 49 # number of analyses

zw = np.linspace(0,zw_over_zj_max*zj,n)
Bp = np.zeros(n)

Cs  = (gamma*TW*L**4)/(pi**4*E*I)
Bpo = 1/(1-Cs)
print(Cs)
print(Bpo)

for i in range(n):
    # Define Steel Beam Object
    beam = steel_beam();

    beam.L  = L
    beam.tw = TW
    beam.zi = zi
    beam.zj = zj
    beam.c  = c

    beam.E  = E
    beam.A  = 1000.0
    beam.I  = I

    beam.alpha  = 1.0
    beam.LF_D   = 1.0
    beam.wd     = qD
    beam.LF_P   = 1.0
    beam.gamma  = gamma
    beam.LF_S1  = 0.0
    beam.LF_S2  = 0.0
    beam.gammas = 0.0
    beam.hs     = 0.0

    beam.nele   = 40
    beam.BuildModel();

    # Run Ponding Analysis
    PA = FE.PondingAnalysis(beam.model,'Constant_Level')
    PA.max_iterations_z = 60
    res = PA.run({'DEAD':1.0},zw[i])
    if res != 0:
        print('Not converged')
    Vmax = beam.Maximum_Shear(PA)
    Mmax = beam.Maximum_Moment(PA)
    
    # Run First-Order Analysis
    PAo = FE.PondingAnalysis(beam.model,'No_Ponding_Effect')
    res = PAo.run({'DEAD':1.0},zw[i])
    if res != 0:
        print('Not converged')
    Vmaxo = beam.Maximum_Shear(PAo)    
    Mmaxo = beam.Maximum_Moment(PAo)
    
    Bp[i] = max(Vmax/Vmaxo,Mmax/Mmaxo)
   

# Make Plot
fig = plt.figure()
ax = fig.add_axes()
line0, = plt.plot([0.0,0.2,0.8,zw_over_zj_max],[1/Bpo,1/Bpo,1.0,1.0],'--', label='Idealized')
line1, = plt.plot(zw/zj,Bp/Bpo,'x-',label='Analysis')
plt.xlabel('zw/zj')
plt.ylabel('Bp/Bpo')
plt.xlim(0.0,zw_over_zj_max)
plt.legend(handles=[line0,line1],frameon=False)
plt.show()
