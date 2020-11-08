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
beam.test_flag = 0
beam.use_CBDI = False

# Run analyses
results1 = beam.RunAnalysis('SimpleStepLevel',target_zw=zw,num_steps=10)
results2 = beam.RunAnalysis('SimpleStepLevel',target_zw=zw,num_steps=100)
results3 = beam.RunAnalysis('SimpleStepLevel',target_zw=zw,num_steps=1000)

num_steps_iterative = 10
water_volume_iter = np.zeros(num_steps_iterative)
water_level_iter  = np.zeros(num_steps_iterative)
print(water_level_iter)
for i in range(num_steps_iterative):
    izw = zw*(i+1)/num_steps_iterative
    iresults = beam.RunAnalysis('IterativeLevel' ,target_zw=izw)
    water_volume_iter[i] = iresults.water_volume
    water_level_iter[i] = iresults.water_level


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
line1, = plt.plot(results1.water_volume/(S*L)/mm, results1.water_level/mm, label='$10$ Steps')
line2, = plt.plot(results2.water_volume/(S*L)/mm, results2.water_level/mm, label='$100$ Steps')
line3, = plt.plot(results3.water_volume/(S*L)/mm, results3.water_level/mm, label='$1000$ Steps')
linei, = plt.plot(water_volume_iter/(S*L)/mm, water_level_iter/mm, 'x', label='Iterative')
plt.xlabel('Normalized Water Volume ($V/SL$, mm)')
plt.ylabel('Water Level (mm)')
plt.legend(handles=[line1,line2,line3,linei],frameon=False)
#plt.xlim(1.5,8.5)
#plt.ylim(0.95,1.65)
plt.savefig('figure_beam_number_of_simple_steps_1.png',dpi=300)
plt.savefig('figure_beam_number_of_simple_steps_1.pdf')


plt.show()