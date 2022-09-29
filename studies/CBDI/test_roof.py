from structures import ExampleRoof
import matplotlib.pyplot as plt
import time

analysis_type = 'IterativeLevel'
#analysis_type = 'SimpleStepLevel'

# Define roof
roof = ExampleRoof()
#roof.include_ponding_effect = False
#roof.plot_load_cells = True

# Run analyses
print('Without CBDI')
roof.use_CBDI = False
tic = time.perf_counter()
results1 = roof.RunAnalysis(analysis_type,target_zw=6,num_steps=100)
toc = time.perf_counter()
print(f"Total time for RunAnalysis: {(toc-tic):0.4f} seconds")

print('With CBDI')
roof.use_CBDI = True
tic = time.perf_counter()
results2 = roof.RunAnalysis(analysis_type,target_zw=6,num_steps=100)
toc = time.perf_counter()
print(f"Total time for RunAnalysis: {(toc-tic):0.4f} seconds")

if analysis_type == 'SimpleStepLevel':
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