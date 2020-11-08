from structures import ExampleBeam
from PyPonding.structures import wf_shapes
import matplotlib.pyplot as plt
import numpy as np

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

inch = 1
mm = 1/25.4*inch

beam = ExampleBeam(d,tw,bf,tf,Fy,E,L,S,qD)
beam.material_type = 'Hardening'
beam.yj = L/48.
beam.kyi = 5.0
beam.kyj = 5.0
beam.test_flag = 0

# Run analyses
beam.num_elements = 20
beam.num_divisions_per_element = 1
beam.element_type = 'dispBeamColumn'
beam.nIP = 3

beam.c = 0.0
results1 = beam.RunAnalysis('SimpleStepLevel',target_zw=zw,num_steps=1000)
beam.c = 1.0
results2 = beam.RunAnalysis('SimpleStepLevel',target_zw=zw,num_steps=1000)
beam.c = 0.0
beam.yi_fixed = False
beam.yj_fixed = False

results3 = beam.RunAnalysis('SimpleStepLevel',target_zw=zw,num_steps=1000)

beam.num_elements = 1
beam.num_divisions_per_element = 20
beam.element_type = 'forceBeamColumn'
beam.nIP = 6
beam.yi_fixed = True
beam.yj_fixed = True

beam.c = 0.0
results4 = beam.RunAnalysis('SimpleStepLevel',target_zw=zw,num_steps=1000)
beam.c = 1.0
results5 = beam.RunAnalysis('SimpleStepLevel',target_zw=zw,num_steps=1000)
beam.c = 0.0
beam.yi_fixed = False
beam.yj_fixed = False
results6 = beam.RunAnalysis('SimpleStepLevel',target_zw=zw,num_steps=1000)

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
line1, = plt.plot(results1.water_volume/(S*L)/mm, results1.water_level/mm, label='$n_{ele} = 20$')
line2, = plt.plot(results2.water_volume/(S*L)/mm, results2.water_level/mm, label='$n_{ele} = 20$ with camber')
line3, = plt.plot(results3.water_volume/(S*L)/mm, results3.water_level/mm, label='$n_{ele} = 20$ with flexible support')
line4, = plt.plot(results4.water_volume/(S*L)/mm, results4.water_level/mm, '--', label='$n_{ele} = 1$')
line5, = plt.plot(results5.water_volume/(S*L)/mm, results5.water_level/mm, '--', label='$n_{ele} = 1$ with camber')
line6, = plt.plot(results6.water_volume/(S*L)/mm, results6.water_level/mm, '--', label='$n_{ele} = 1$ with flexible support')
plt.xlabel('Normalized Water Volume ($V/SL$, mm)')
plt.ylabel('Water Level (mm)')
plt.legend(handles=[line1,line2,line3,line4,line5,line6],frameon=False)
#plt.xlim(1.5,8.5)
#plt.ylim(0.95,1.65)
plt.savefig('figure_beam_camber_and_flexible_supports.png',dpi=300)
plt.savefig('figure_beam_camber_and_flexible_supports.pdf')

plt.show()