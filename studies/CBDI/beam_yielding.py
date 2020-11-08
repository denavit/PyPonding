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

inch = 1
mm = 1/25.4*inch

Vw = L*S*300*mm

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
beam.nIP = 6

beam.residual_stress_ratio = 0.0
results1 = beam.RunAnalysis('SimpleStepVolume',target_Vw=Vw,num_steps=1000)
beam.residual_stress_ratio = -0.3
results2 = beam.RunAnalysis('SimpleStepVolume',target_Vw=Vw,num_steps=1000)

beam.num_elements = 1
beam.num_divisions_per_element = 20
beam.element_type = 'forceBeamColumn'
beam.nIP = 3

beam.residual_stress_ratio = 0.0
results3 = beam.RunAnalysis('SimpleStepVolume',target_Vw=Vw,num_steps=1000)
beam.residual_stress_ratio = -0.3
results4 = beam.RunAnalysis('SimpleStepVolume',target_Vw=Vw,num_steps=1000)

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
line2, = plt.plot(results2.water_volume/(S*L)/mm, results2.water_level/mm, label='$n_{ele} = 20$ with residual stress')
line3, = plt.plot(results3.water_volume/(S*L)/mm, results3.water_level/mm, '--', label='$n_{ele} = 1$')
line4, = plt.plot(results4.water_volume/(S*L)/mm, results4.water_level/mm, '--', label='$n_{ele} = 1$ with residual stress')
plt.xlabel('Normalized Water Volume ($V/SL$, mm)')
plt.ylabel('Water Level (mm)')
plt.legend(handles=[line1,line2,line3,line4],frameon=False)
plt.xlim(0.0,300.0)
plt.ylim(0.0,350.0)
plt.savefig('figure_beam_yielding.png',dpi=300)
plt.savefig('figure_beam_yielding.pdf')

plt.show()