import openseespy.opensees as ops
from math import pi,cos,cosh,ceil
import numpy as np
import matplotlib.pyplot as plt
from wide_flange import wf,wf_shapes
from PyPonding import FE
  
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
Hk      = 1.0e-4*E

TW      = 5*ft

shape_name = 'W14X22';
L       = 30*ft
slope   = 0.5*in_per_ft
qD      = 10.0*psf

    
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

# Set OpenSees analysis options
wf_section.material_type = 'Hardening'
wf_section.max_volume = (500*inch)*wf_section.L*wf_section.TW
wf_section.num_steps = 10000
wf_section.vol_tol = wf_section.max_volume/wf_section.num_steps/10000.

# Run OpenSees Analysis 
(data_volume_1,data_height_1) = wf_section.perform_OpenSees_analysis();

# Run PyPonding Analyses
num_analyses = 50
data_height_2 = np.linspace(0.0,1.1*np.amax(data_height_1),num=num_analyses)
data_volume_2 = np.empty(num_analyses)
elastic_beam = wf_section.steel_beam_object()
PA = FE.PondingAnalysis(elastic_beam.model,'Constant_Level')
PA.max_iterations_z = 60
for i in range(len(data_height_2)):
    res = PA.run({'DEAD':1.0},data_height_2[i])
    if res == 0:
        (V,dVdz) = elastic_beam.model.GetPondingVolume(PA.d,data_height_2[i])
        data_volume_2[i] = V
    else:
        print('Not converged')


# Plot Results
plt.figure(1)  
line1, = plt.plot(data_volume_1, data_height_1, '-')
line2, = plt.plot(data_volume_2, data_height_2, '-')
plt.legend((line1, line2), ('Analysis 1', 'Analysis 2'))
plt.xlabel('Water Volume')
plt.ylabel('Water Height')

plt.show()