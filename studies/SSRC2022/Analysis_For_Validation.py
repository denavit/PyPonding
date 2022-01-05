import matplotlib.pyplot as plt
# import openseespy.opensees as ops
from PyPonding.structures import wf


from aisc import angle_database, double_angle_database, wide_flange_database

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
Fy = 50.0*ksi
E = 29000.0*ksi
Hk = 1.0e-4*E
TW = 10*ft
shape_name = 'W14X22'
L = 40*ft
slope = 0.0*in_per_ft
qD = 20.0*psf
gamma = 62.4*pcf

# Lookup shape data
shape_data = wide_flange_database[shape_name]
d = shape_data['d']*inch
bf = shape_data['bf']*inch
tf = shape_data['tf']*inch
tw = shape_data['tw']*inch

# Create wide-flange object
wf_section = wf(d, tw, bf, tf, Fy, E, Hk)
wf_section.L = L
wf_section.gamma = gamma
wf_section.TW = TW
wf_section.zi = 0.0*inch
wf_section.zj = L*slope
wf_section.wD = qD*TW  # Dead load per length
wf_section.material_type = 'Elastic'
wf_section.num_steps = 2000
wf_section.percent_drop = 1.0

# Run OpenSees analyses
wf_section.max_volume = (6*inch)*L*TW
wf_section.vol_tol = wf_section.max_volume/wf_section.num_steps/10000.

(data_volume, data_height) = wf_section.perform_OpenSees_analysis()


# Plot Results
fig = plt.figure()
plt.plot(data_volume/(L*TW), data_height, 'k-')
plt.xlabel('Normalized Water Volume, V/(L*TW) (in)')
plt.ylabel('Water Level (in)')

plt.show()
