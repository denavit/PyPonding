from structures import ExampleBeam
import matplotlib.pyplot as plt

# Define beam
beam = ExampleBeam()
beam.include_ponding_effect = False
beam.yi = 0.
beam.yj = 0.
beam.c = 0.
beam.xj_fixed = True
beam.transf_type = 'Corotational'

# Run analyses
beam.use_CBDI = False
results1 = beam.RunAnalysis()
beam.use_CBDI = True
results2 = beam.RunAnalysis()

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