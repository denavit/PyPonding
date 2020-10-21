import matplotlib.pyplot as plt
import numpy as np
from structures import ExampleBeam
from PyPonding.structures import wf_shapes
from math import pi,tan,tanh

# Define beam
shape_data = wf_shapes['W14X22']
d  = shape_data['d']
tw = shape_data['tw']
bf = shape_data['bf']
tf = shape_data['tf']
Fy = 50
E  = 29000
L  = 40*12
S  = 10*12
qD = 0./1000/12**2    # Uniform dead load (forcer per unit area, downward positive)
zw = 6

beam = ExampleBeam(d,tw,bf,tf,Fy,E,L,S,qD)
beam.material_type = 'ElasticSection'
beam.tol_load_level = 1.0e-6
beam.num_elements = 1
beam.num_divisions_per_element = 20
beam.element_type = 'forceBeamColumn'

# Analytical solution
I = beam.Iz
gamma = beam.gamma
C = (gamma*S*L**4)/(pi**4*E*I)
#amp = 1/(1-C)
amp = (tan(0.5*pi*C**0.25)+tanh(0.5*pi*C**0.25))/(pi*C**0.25)

# Initilize results
num_ip = np.array([2,3,4,5,6,7,8])
water_volume_1 = np.zeros(np.size(num_ip))
water_volume_2 = np.zeros(np.size(num_ip))
water_volume_3 = np.zeros(np.size(num_ip))
water_volume_4 = np.zeros(np.size(num_ip))

# Run analyses
for i in range(len(num_ip)):
    beam.nIP = int(num_ip[i])
    
    # Legendre
    beam.beamint_type = 'Legendre'
    results = beam.RunAnalysis('IterativeLevel' ,target_zw=zw)
    water_volume_1[i] = results.water_volume
    
    # Lobatto
    beam.beamint_type = 'Lobatto'
    results = beam.RunAnalysis('IterativeLevel' ,target_zw=zw)
    water_volume_2[i] = results.water_volume

    # NewtonCotes (Open)
    beam.beamint_type = 'OpenNewtonCotes'
    results = beam.RunAnalysis('IterativeLevel' ,target_zw=zw)
    water_volume_3[i] = results.water_volume

    # Mid-Distance
    beam.beamint_type = 'MidDistance'
    results = beam.RunAnalysis('IterativeLevel' ,target_zw=zw)
    water_volume_4[i] = results.water_volume

print(amp)
    
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
line0, = plt.plot([0,35], [amp,amp],'k-',linewidth=1,label='Analytical')
line1, = plt.plot(num_ip, water_volume_1/(S*L*zw),'b.',label='Legendre')
line2, = plt.plot(num_ip, water_volume_2/(S*L*zw),'rx',label='Lobotto')
line3, = plt.plot(num_ip, water_volume_3/(S*L*zw),'g*',label='Newton-Cotes (Open)')
line4, = plt.plot(num_ip, water_volume_4/(S*L*zw),'+c',label='Mid-Distance')
plt.legend(handles=[line0,line1,line2,line3,line4],frameon=False)
plt.xlabel('Number of Integration Points')
plt.ylabel('Normalized Water Volume ($V/SLz_w$)')
plt.xlim(1.5,8.5)
#plt.ylim(0.95,1.65)
plt.savefig('Beam_Numerical_Integration_1.png',dpi=300)
plt.savefig('Beam_Numerical_Integration_1.pdf')

plt.show()
