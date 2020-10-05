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

I = beam.Iz
gamma = beam.gamma

C = (gamma*S*L**4)/(pi**4*E*I)
#amp = 1/(1-C)
amp = (tan(0.5*pi*C**0.25)+tanh(0.5*pi*C**0.25))/(pi*C**0.25)

# Initilize results
num_divisions = np.array([1,2,4,8,12,16,20,24,28,32])
water_volume_1 = np.zeros(np.size(num_divisions))
water_volume_2 = np.zeros(np.size(num_divisions))
water_volume_3 = np.zeros(np.size(num_divisions))

# Run analyses
for i in range(len(num_divisions)):
    # num_elements = num_divisions
    beam.num_elements = int(num_divisions[i])
    beam.num_divisions_per_element = 1
    beam.element_type = 'dispBeamColumn'
    beam.nIP = 3
    results = beam.RunAnalysis('IterativeLevel' ,target_zw=zw)
    water_volume_1[i] = results.water_volume
    
    # num_elements = 1
    beam.num_elements = 1
    beam.num_divisions_per_element = int(num_divisions[i])
    beam.element_type = 'forceBeamColumn'
    beam.nIP = 6
    results = beam.RunAnalysis('IterativeLevel' ,target_zw=zw)
    water_volume_2[i] = results.water_volume

    # num_elements = 1
    beam.num_elements = 1
    beam.num_divisions_per_element = int(num_divisions[i])
    beam.element_type = 'forceBeamColumn'
    beam.nIP = 4
    results = beam.RunAnalysis('IterativeLevel' ,target_zw=zw)
    water_volume_3[i] = results.water_volume


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
line1, = plt.plot(num_divisions, water_volume_1/(S*L*zw),'b.',label='$n_{ele} = n_{div}$')
line2, = plt.plot(num_divisions, water_volume_2/(S*L*zw),'rx',label='$n_{ele} = 1, n_{ip} = 6$')
line3, = plt.plot(num_divisions, water_volume_3/(S*L*zw),'g*',label='$n_{ele} = 1, n_{ip} = 4$')
plt.legend(handles=[line0,line1,line2,line3],frameon=False)
plt.xlabel('Number of Divisions')
plt.ylabel('Normalized Water Volume ($V/SLz_w$)')
plt.xlim(0,35)
plt.ylim(0.95,1.65)
plt.savefig('Beam_Mesh_Refinement_1.png',dpi=300)
plt.savefig('Beam_Mesh_Refinement_1.pdf')

plt.show()
