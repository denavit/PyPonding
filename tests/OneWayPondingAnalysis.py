#import sys
#sys.path.append('/home/mhscott/OpenSees/SRC/interpreter')
from PyPonding.structures.steel_beam import steel_beam
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


# Input parameters
E       = 29000.0*ksi   # Modulus of elasticity
I       = 75*inch**4   # Moment of inertia

L       = 40*ft         # Beam span
TW      = 5*ft          # Tributary width
c       = 2*inch        # Camber

qD      = 0.0*psf      # Dead load (force per area)
gamma   = 62.4*pcf      # Unit weight of water

zi      = 0*inch        # Elevation of left end
zj      = 0*inch        # Elevation of right end
zw      = 2*inch        # Elevation of water



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
res = PA.run({'DEAD':1.0},zw)
if res != 0:
    print('Not converged')
Vmax = beam.Maximum_Shear(PA)
Mmax = beam.Maximum_Moment(PA)
(x,M,V) = beam.Moment_and_Shear(PA)
x_at_Mmax = x[np.argmax(np.absolute(M))]


# Run First-Order Analysis
PAo = FE.PondingAnalysis(beam.model,'No_Ponding_Effect')
res = PAo.run({'DEAD':1.0},zw)
if res != 0:
    print('Not converged')
Vmaxo = beam.Maximum_Shear(PAo)    
Mmaxo = beam.Maximum_Moment(PAo)
(x,M,V) = beam.Moment_and_Shear(PAo)
x_at_Mmaxo = x[np.argmax(np.absolute(M))]


# Print Results
Cs = (gamma*TW*L**4)/(pi**4*E*I)
print('\nCs = %.3f' % Cs)
print('Basic Amplification 1/(1-Cs) = %.3f' % (1/(1-Cs)))
x = pi/2*Cs**0.25
print('More Exact Amplification     = %.3f\n' % ((1/cos(x)-1/cosh(x))/x**2)) # from Silver (2010)
print('                      with ponding           no ponding          amplification')
print('Max Shear           %9.3f kip         %9.3f kip          %9.3f' % (Vmax,Vmaxo,Vmax/Vmaxo))
print('Max Moment          %9.3f kip-in      %9.3f kip-in       %9.3f' % (Mmax,Mmaxo,Mmax/Mmaxo))
print('Pos. of Mmax      %5.1f in (%4.3f*L)    %5.1f in (%4.3f*L)\n' % (x_at_Mmax,x_at_Mmax/L,x_at_Mmaxo,x_at_Mmaxo/L))


