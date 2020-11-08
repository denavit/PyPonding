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

inch = 1
mm = 1/25.4*inch

Vw = L*S*350*mm

beam = ExampleBeam(d,tw,bf,tf,Fy,E,L,S,qD)
beam.material_type = 'Hardening'
beam.num_elements = 1
beam.num_divisions_per_element = 20
beam.element_type = 'forceBeamColumn'
beam.yj = L/48.

# Initilize results
num_ip = np.array([2,3,4,5,6,7,8,9,10])
maximum_water_level_1 = np.zeros(np.size(num_ip))
maximum_water_level_2 = np.zeros(np.size(num_ip))
maximum_water_level_3 = np.zeros(np.size(num_ip))
maximum_water_level_4 = np.zeros(np.size(num_ip))

# Run analyses
for i in range(len(num_ip)):
    beam.nIP = int(num_ip[i])
    
    # Legendre
    beam.beamint_type = 'Legendre'
    results = beam.RunAnalysis('SimpleStepVolume',target_Vw=Vw,num_steps=1000)
    maximum_water_level_1[i] = max(results.water_level)
    
    # Lobatto
    beam.beamint_type = 'Lobatto'
    results = beam.RunAnalysis('SimpleStepVolume',target_Vw=Vw,num_steps=1000)
    maximum_water_level_2[i] = max(results.water_level)

    # NewtonCotes (Open)
    beam.beamint_type = 'OpenNewtonCotes'
    results = beam.RunAnalysis('SimpleStepVolume',target_Vw=Vw,num_steps=1000)
    maximum_water_level_3[i] = max(results.water_level)

    # Mid-Distance
    beam.beamint_type = 'MidDistance'
    results = beam.RunAnalysis('SimpleStepVolume',target_Vw=Vw,num_steps=1000)
    maximum_water_level_4[i] = max(results.water_level)

   
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
line1, = plt.plot(num_ip, maximum_water_level_1/mm,'b.--',label='Legendre')
line2, = plt.plot(num_ip, maximum_water_level_2/mm,'rx',label='Lobotto')
line3, = plt.plot(num_ip, maximum_water_level_3/mm,'g*',label='Newton-Cotes (Open)')
line4, = plt.plot(num_ip, maximum_water_level_4/mm,'+c',label='Mid-Distance')
plt.legend(handles=[line1,line2,line3,line4],frameon=False)
plt.xlabel('Number of Integration Points')
plt.ylabel('Maximum Water Level (mm)')
plt.xlim(1.5,10.5)
#plt.ylim(0.95,1.65)
plt.savefig('figure_beam_numerical_integration_2.png',dpi=300)
plt.savefig('figure_beam_numerical_integration_2.pdf')

plt.show()
