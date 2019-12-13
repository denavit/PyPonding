import openseespy.opensees as ops
from math import pi,cos,cosh,ceil,tan,tanh
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
Hk      = 1.0e-4*E
TW      = 10*ft
shape_name = 'W14X22';
L1      = 40*ft
L2      = 49.57*ft
L3      = 54.85*ft
slope   = 0.0*in_per_ft
qD      = 0.0*psf
 
# Lookup shape data
shape_data = wf_shapes[shape_name]
d  = shape_data['d']*inch
bf = shape_data['bf']*inch
tf = shape_data['tf']*inch
tw = shape_data['tw']*inch

# Create wide-flange object
wf_section = wf(d,tw,bf,tf,Fy,E,Hk)

# Set additional properties
wf_section.L        = L1
wf_section.gamma    = 62.4*pcf
wf_section.TW       = TW
wf_section.zi       = 0.0*inch
wf_section.zj       = L1*slope
wf_section.wD       = qD*TW # Dead load per length

# Set OpenSees analysis options
wf_section.material_type = 'Elastic'
wf_section.num_steps = 2000


# Run OpenSees analyses 
# L1
wf_section.L  = L1
C1 = wf_section.C()
wf_section.max_volume = (14*inch)*L1*TW
wf_section.vol_tol = wf_section.max_volume/wf_section.num_steps/10000.

wf_section.zj = (0.00)*L1
(data_volume_L1,data_height_L1) = wf_section.perform_OpenSees_analysis();

wf_section.zj = (0.25*in_per_ft)*L1
(data_volume_L1s,data_height_L1s) = wf_section.perform_OpenSees_analysis();

wf_section.material_type = 'Hardening'
(data_volume_L1si,data_height_L1si) = wf_section.perform_OpenSees_analysis();

wf_section.frc = -0.3*Fy
(data_volume_L1sirs,data_height_L1sirs) = wf_section.perform_OpenSees_analysis();

# L2
wf_section.material_type = 'Elastic'
wf_section.frc = 0.0
wf_section.L  = L2
C2 = wf_section.C()
wf_section.max_volume = (10*inch)*L2*TW
wf_section.vol_tol = wf_section.max_volume/wf_section.num_steps/10000.

wf_section.zj = (0.00)*L2
(data_volume_L2,data_height_L2) = wf_section.perform_OpenSees_analysis();

wf_section.zj = (0.25*in_per_ft)*L2
(data_volume_L2s,data_height_L2s) = wf_section.perform_OpenSees_analysis();

# L3
wf_section.L  = L3
C3 = wf_section.C()
wf_section.max_volume = (10*inch)*L3*TW
wf_section.vol_tol = wf_section.max_volume/wf_section.num_steps/10000.
wf_section.percent_drop = 10000.0

wf_section.zj = (0.00)*L3
(data_volume_L3,data_height_L3) = wf_section.perform_OpenSees_analysis();

wf_section.zj = (0.25*in_per_ft)*L3
(data_volume_L3s,data_height_L3s) = wf_section.perform_OpenSees_analysis();


slope1 = (pi*C1**0.25)/(TW*L1*(tan(0.5*pi*C1**0.25)+tanh(0.5*pi*C1**0.25)))
slope2 = (pi*C2**0.25)/(TW*L2*(tan(0.5*pi*C2**0.25)+tanh(0.5*pi*C2**0.25)))
slope3 = (pi*C3**0.25)/(TW*L3*(tan(0.5*pi*C3**0.25)+tanh(0.5*pi*C3**0.25)))

# Plot Results
plt.figure(1)  
line1, = plt.plot(data_volume_L1/(TW*L1), data_height_L1, '-')
line2, = plt.plot(data_volume_L2/(TW*L2), data_height_L2, '-')
line3, = plt.plot(data_volume_L3/(TW*L3), data_height_L3, '-')

line4, = plt.plot(np.array([0,10]), np.array([0,10*TW*L1*slope1]), 'o')
line5, = plt.plot(np.array([0,10]), np.array([0,10*TW*L2*slope2]), 'o')
line6, = plt.plot(np.array([0,10]), np.array([0,10*TW*L3*slope3]), 'o')

print(C1)
print(C2)
print(C3)

plt.legend((line1, line2, line3), ('C = 0.424', 'C = 1.0', 'C = 1.5'))
plt.xlabel('Normalized Water Volume (V/SL)')
plt.ylabel('Water Height')
plt.xlim( 0,10)
plt.ylim(-7,7)



plt.figure(2)  
line1, = plt.plot(data_volume_L1s/(TW*L1), data_height_L1s, '-')
line2, = plt.plot(data_volume_L2s/(TW*L2), data_height_L2s, '-')
line3, = plt.plot(data_volume_L3s/(TW*L3), data_height_L3s, '-')

plt.legend((line1, line2, line3), ('C = 0.424', 'C = 1.0', 'C = 1.5'))
plt.xlabel('Normalized Water Volume (V/SL)')
plt.ylabel('Water Height')
plt.xlim(0,10)
plt.ylim(0,12)



plt.figure(3)  
line1, = plt.plot(   data_volume_L1s/(TW*L1),    data_height_L1s, '-')
line2, = plt.plot(  data_volume_L1si/(TW*L1),   data_height_L1si, '-')
line3, = plt.plot(data_volume_L1sirs/(TW*L1), data_height_L1sirs, '-')

plt.legend((line1, line2, line3), ('Elastic', 'Inelastic', 'Inelastic with RS'))
plt.xlabel('Normalized Water Volume (V/SL)')
plt.ylabel('Water Height')
plt.xlim(0,14)
plt.ylim(0,14)

plt.show()
