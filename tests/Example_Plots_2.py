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
mm  = 1/25.4*inch
m   = 1000*mm
kN  = 1/4.448222*kip
kNm = kN*m

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
gamma   = 62.4*pcf

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
wf_section.gamma    = gamma
wf_section.TW       = TW
wf_section.zi       = 0.0*inch
wf_section.zj       = L1*slope
wf_section.wD       = qD*TW # Dead load per length

# Set OpenSees analysis options
wf_section.material_type = 'Elastic'
wf_section.num_steps = 2000
wf_section.percent_drop = 1.0

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
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=8)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)

fig = plt.figure(figsize=(3.25,2.50))
ax = fig.add_axes([0.18,0.18,0.78,0.78])
line1, = plt.plot(data_volume_L1/(TW*L1)/mm, data_height_L1/mm, 'k-')
line2, = plt.plot(data_volume_L2/(TW*L2)/mm, data_height_L2/mm, 'k--')
line3, = plt.plot(data_volume_L3/(TW*L3)/mm, data_height_L3/mm, 'k-.')
x_closedform = np.linspace(0.0, 250.0, num=10)
line4, = plt.plot(x_closedform,x_closedform*TW*L1*slope1,'ko', markersize=4)
line5, = plt.plot(x_closedform,x_closedform*TW*L2*slope2,'ko', markersize=4)
#line6, = plt.plot(x_closedform,x_closedform*TW*L3*slope3,'ko', markersize=4)
print(C1)
print(C2)
print(C3)
plt.legend((line4,line1,line2,line3), ('Closed-form Solution','C = 0.424','C = 1.0','C = 1.5'),frameon=False)
plt.xlabel('Normalized Water Volume (V/SL, mm)')
plt.ylabel('Water Level (mm)')
plt.xlim( 0,250)
plt.ylim(-250,250)
plt.savefig('Example_Plot_1.png',dpi=300)
plt.savefig('Example_Plot_1.pdf')

fig = plt.figure(figsize=(3.25,2.50))
ax = fig.add_axes([0.15,0.18,0.80,0.78])
line1, = plt.plot(data_volume_L1s/(TW*L1)/mm, data_height_L1s/mm, 'k-')
line2, = plt.plot(data_volume_L2s/(TW*L2)/mm, data_height_L2s/mm, 'k--')
line3, = plt.plot(data_volume_L3s/(TW*L3)/mm, data_height_L3s/mm, 'k-.')
plt.legend((line1, line2, line3), ('C = 0.424', 'C = 1.0', 'C = 1.5'),frameon=False)
plt.xlabel('Normalized Water Volume (V/SL, mm)')
plt.ylabel('Water Level (mm)')
plt.xlim(0,250)
plt.ylim(0,300)
plt.savefig('Example_Plot_2.png',dpi=300)
plt.savefig('Example_Plot_2.pdf')

fig = plt.figure(figsize=(7.00,2.50))
ax = fig.add_axes([0.075,0.21,0.40,0.75])
line1, = plt.plot(   data_volume_L1s/(TW*L1)/mm,    data_height_L1s/mm, 'k-')
line2, = plt.plot(  data_volume_L1si/(TW*L1)/mm,   data_height_L1si/mm, 'k--')
line3, = plt.plot(data_volume_L1sirs/(TW*L1)/mm, data_height_L1sirs/mm, 'r-.')
plt.legend((line1, line2, line3), ('Elastic', 'Inelastic (no residual stress)', 'Inelastic'),frameon=False)
plt.xlabel('Normalized Water Volume (V/SL, mm)\n(a)')
plt.ylabel('Water Level (mm)')
plt.xlim(0,350)
plt.ylim(0,350)

ax = fig.add_axes([0.575,0.21,0.40,0.75])
line1, = plt.plot(   data_volume_L1s/(TW*L1)/mm,    data_height_L1s/mm, 'k-')
line2, = plt.plot(  data_volume_L1si/(TW*L1)/mm,   data_height_L1si/mm, 'k--')
line3, = plt.plot(data_volume_L1sirs/(TW*L1)/mm, data_height_L1sirs/mm, 'r-.')
plt.legend((line1, line2, line3), ('Elastic', 'Inelastic (no residual stress)', 'Inelastic'),frameon=False)
plt.xlabel('Normalized Water Volume (V/SL, mm)\n(b)')
plt.ylabel('Water Level (mm)')
plt.xlim(210,290)
plt.ylim(250,320)

plt.savefig('Example_Plot_3.png',dpi=300)
plt.savefig('Example_Plot_3.pdf')

plt.show()
