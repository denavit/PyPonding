from structures import ExampleRoof
import matplotlib.pyplot as plt

# Define roof
roof = ExampleRoof()
roof.num_steps_zw = 25
roof.include_ponding_effect = False
#roof.plot_load_cells = True

# Run analyses
roof.use_CBDI = False
results1 = roof.RunAnalysis()
roof.use_CBDI = True
results2 = roof.RunAnalysis()

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
ax.plot(results1.water_volume, results1.col_react_B2, label='Nodal')
ax.plot(results2.water_volume, results2.col_react_B2, '--', label='CBDI')
plt.xlabel('Water Volume (in^3)')
plt.ylabel('Column Reaction (kips)')
ax.legend()

plt.show()