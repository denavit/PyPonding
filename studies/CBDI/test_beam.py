from structures import ExampleBeam
from PyPonding.structures import wf_shapes
import matplotlib.pyplot as plt

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

beam = ExampleBeam(d,tw,bf,tf,Fy,E,L,S,qD)
beam.material_type = 'Hardening'
#beam.include_ponding_effect = False
beam.yj = L/48.

# Run analyses
beam.use_CBDI = False
results1 = beam.RunAnalysis('SimpleStepVolume',target_Vw=20*L*S,num_steps=1000)
beam.use_CBDI = True
results2 = beam.RunAnalysis('SimpleStepVolume',target_Vw=20*L*S,num_steps=1000)

# Plot results
fig = plt.figure()
ax = plt.axes()
ax.plot(results1.water_volume, results1.water_level, label='Nodal')
ax.plot(results2.water_volume, results2.water_level, '--', label='CBDI')
plt.xlabel('Water Volume (in^3)')
plt.ylabel('Water Level (in)')
ax.legend()

fig = plt.figure()
ax = plt.axes()
ax.plot(results1.water_volume, results1.Ryi, label='Nodal')
ax.plot(results2.water_volume, results2.Ryi, '--', label='CBDI')
plt.xlabel('Water Volume (in^3)')
plt.ylabel('Left End Vertical Reaction (kips)')
ax.legend()

fig = plt.figure()
ax = plt.axes()
ax.plot(results1.water_volume, results1.Rxi, label='Nodal')
ax.plot(results2.water_volume, results2.Rxi, '--', label='CBDI')
plt.xlabel('Water Volume (in^3)')
plt.ylabel('Left End Horizontal Reaction (kips)')
ax.legend()

plt.show()