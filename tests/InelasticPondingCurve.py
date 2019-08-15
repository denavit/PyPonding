import openseespy.opensees as ops
from math import pi,cos,cosh,ceil
import numpy as np
import matplotlib.pyplot as plt
from wide_flange import wf,wf_shapes

  
# Define units
inch = 1.0
kip = 1.0

lb = kip/1000.0
ft = 12.0*inch

in_per_ft = inch/ft

ksi = kip/inch**2
psf = lb/ft**2
pcf = psf/ft


# Input parameters
Fy      = 50.0*ksi
E       = 29000.0*ksi
Hk      = E/1000.0
TW      = 5*ft
shape_name = 'W14X22';
qD      = 0.0*psf
gamma   = 62.4*pcf


# Lookup shape data
shape_data = wf_shapes[shape_name]
d  = shape_data['d']*inch
bf = shape_data['bf']*inch
tf = shape_data['tf']*inch
tw = shape_data['tw']*inch
Ix = shape_data['Ix']*inch**4
Zx = shape_data['Zx']*inch**3

Mp = Fy*Zx

# Create wide-flange object
wf_section = wf(d,tw,bf,tf,Fy,E,Hk)

# Set additional properties
wf_section.gamma    = gamma
wf_section.TW       = TW
wf_section.wD       = qD*TW # Dead load per length

# Set OpenSees analysis options
wf_section.material_type = 'Hardening'

L_cr = ((pi**4*E*Ix)/(gamma*TW))**0.25

max_length = 1.1*L_cr
L_OPS    = np.arange(10*ft,max_length,2*ft)
zmax_OPS = np.empty(L_OPS.shape)
hp_OPS   = np.empty(L_OPS.shape)

# Run OpenSees analyses 
for iL in range(L_OPS.size):
    wf_section.L = L_OPS[iL]
    
    wf_section.num_steps = 5000
    wf_section.max_volume = (500*inch)*wf_section.L*wf_section.TW
    wf_section.vol_tol = wf_section.max_volume/wf_section.num_steps/100.
    
    (data_volume,data_height) = wf_section.perform_OpenSees_analysis();
    zmax_OPS[iL] = np.max(data_height)
    hp_OPS[iL] = (8*Mp)/(gamma*TW*L_OPS[iL]**2)

#print(zmax_OPS)
#print(hp_OPS)
    
plt.plot([0,L_cr,L_cr],[1,1,0],'-')
plt.plot(L_OPS, zmax_OPS/hp_OPS, 'ro') 
plt.xlim(0,1.1*L_cr)
plt.ylim(0,1.1)
plt.xlabel('Beam length (in)')
plt.ylabel('Normalized height (h/h_p)')
plt.show()