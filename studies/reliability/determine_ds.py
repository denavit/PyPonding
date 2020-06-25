from math import pi,cos,cosh,ceil,sqrt
import numpy as np
import matplotlib.pyplot as plt
from wide_flange import wf,wf_shapes
from PyPonding import FE
  
# Define units
inch = 1.0
kip  = 1.0
sec  = 1.0

lb = kip/1000.0
ft = 12.0*inch

in_per_ft = inch/ft

ksi = kip/inch**2
psf = lb/ft**2
pcf = psf/ft

hour = 60*60*sec
in_per_hour = inch/hour

in_per_ss = inch/sec**2

# Input parameters
Fy      = 50.0*ksi
E       = 29000.0*ksi
Hk      = 1.0e-4*E

TW      = 10*ft

shape_name = 'W14X22';
L       = 30*ft
slope   = 0.0*in_per_ft
qD      = 15.0*psf

  
# Lookup shape data
shape_data = wf_shapes[shape_name]
d  = shape_data['d']*inch
bf = shape_data['bf']*inch
tf = shape_data['tf']*inch
tw = shape_data['tw']*inch

# Create wide-flange object
wf_section = wf(d,tw,bf,tf,Fy,E,Hk)

# Set additional properties
wf_section.L        = L 
wf_section.gamma    = 62.4*pcf
wf_section.TW       = TW
wf_section.zi       = 0.0*inch
wf_section.zj       = L*slope
wf_section.wD       = qD*TW # Dead load per length

# Set design method and compute zw_lim (ds+dh)
#method = 'AISC Appendix 2'
#method = 'DAMP'
method = 'Proposed for ASCE 7'
zw_lim = wf_section.maximum_permitted_zw(method)

# Compute hydraulic head
i  = 3.75*in_per_hour  # Rain intensity
Ss = 40*ft   # Scupper spacing
ws = 6*inch  # Scupper width
A  = Ss*L    # Tributary area (per scupper)
Q  = A*i     # Scupper flow rate
g  = 386.4*in_per_ss
cd = 0.6
dh = ((3*Q)/(2*cd*ws*sqrt(2*g)))**(2/3)

# Print Results
print('Ix = %.3f in^4' % (wf_section.Iz()))
print('Zx = %.3f in^3' % (wf_section.Zz()))
print('Sx = %.3f in^3' % (wf_section.Sz()))

print('dh = %.3f in' % (dh))
print('ds = %.3f in' % (zw_lim-dh))

