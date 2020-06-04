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
Hk      = E/10000.0

TW      = 10.0*ft

shape_name = 'W24X55';
L       = 40*ft
slope   = 0.25*in_per_ft
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

# Run OpenSees analysis 
(data_volume,data_height) = wf_section.perform_OpenSees_analysis();
zmax_inelastic = np.max(data_height)


# Elastic Analyses
elastic_beam = wf_section.steel_beam_object()
zmax_elastic = elastic_beam.Run_To_Strength_Limit(start_level = 1.0,max_level=200.0)
tau = elastic_beam.determine_stiffness_reduction(zmax_inelastic)

# Print Data
print('\nz_max from the inelastic analysis: %.4f in.' % (zmax_inelastic))
print('    (this is the peak height from the volume vs. height graph)')
print('z_max from an elastic analysis with nominal EI: %.4f in.' % (zmax_elastic))
print('    (this is the height that causes M_max = M_p)')
print('The necessary stiffness reduction is %.4f' % (tau))
print('    (this is the factor that when applied to EI in an elastic model')
print('    results in M_max = M_p at z_max from the inelastic analysis)\n')
