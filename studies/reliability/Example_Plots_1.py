import openseespy.opensees as ops
from math import pi,ceil
import numpy as np
from numpy import cos,cosh
import matplotlib.pyplot as plt
from PyPonding.structures import wf,wf_shapes

# Define units
inch = 1.0
kip = 1.0
lb  = kip/1000.0
ft  = 12.0*inch
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
L       = 40*ft
qD      = 10.0*psf
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
wf_section.L        = L
wf_section.gamma    = gamma
wf_section.TW       = TW
wf_section.zi       = 0.0*inch
wf_section.zj       = 0.0*inch
wf_section.wD       = qD*TW # Dead load per length

# Set OpenSees analysis options
wf_section.material_type = 'Elastic'
wf_section.num_steps     = 2000
wf_section.max_volume    = 249450
wf_section.vol_tol       = wf_section.max_volume/wf_section.num_steps/10000.
wf_section.extra_results = True

# Run OpenSees analyses 
(data_volume,data_height,extra_results) = wf_section.perform_OpenSees_analysis();
x_OPS = extra_results[0]
y_OPS = extra_results[1]
M_OPS = extra_results[2]

# Compute Closed-form Results
C  = wf_section.C()
I  = wf_section.Iz()
wi = (qD + gamma*(2*inch))*TW

x  = np.linspace(0.0, L, num=100)

y1 = -wi*x*(L**3-2*L*x**2+x**3)/(24*E*I)
M1 = 0.5*wi*x*(L-x)

fact = pi*C**0.25
y2 = -wi/(2*gamma*TW)*( (cos(fact*(0.5-x/L))/cos(fact*0.5) - 1) - (1 - cosh(fact*(0.5-x/L))/cosh(fact*0.5)) )
M2 = (wi*L**2)/(2*pi**2*C**0.5)*( cos(fact*(0.5-x/L))/cos(0.5*fact) - cosh(fact*(0.5-x/L))/cosh(0.5*fact) )


# Plot Results
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=8)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)

fig = plt.figure(figsize=(7.00,2.50))
ax = fig.add_axes([0.09,0.21,0.40,0.75])
line2, = plt.plot(x/m,y2/mm,'k-')
line3, = plt.plot(x/m,y1/mm,'k--')
line1, = plt.plot(x_OPS/m,y_OPS/mm,'ko')
plt.legend((line1, line2, line3), ('Numerical Results', 'Closed-form Solution', 'No Ponding Effect'),frameon=False)
plt.xlabel('Position Along Length of Beam (m)\n(a)')
plt.ylabel('Displacement (mm)')
plt.xlim(0,L/m)
plt.ylim(-100,0)

ax = fig.add_axes([0.59,0.21,0.40,0.75])
line2, = plt.plot(x/m,M2/kNm,'k-')
line3, = plt.plot(x/m,M1/kNm,'k--')
line1, = plt.plot(x_OPS/m,M_OPS/kNm,'ko')
plt.legend((line1, line2, line3), ('Numerical Results', 'Closed-form Solution', 'No Ponding Effect'),frameon=False)
plt.xlabel('Position Along Length of Beam (m)\n(b)')
plt.ylabel('Bending Moment (kN-m)')
plt.xlim(0,L/m)
plt.ylim(0,100)
plt.savefig('Example_Defl_and_Moment.png',dpi=300)
plt.savefig('Example_Defl_and_Moment.pdf')

plt.show()
