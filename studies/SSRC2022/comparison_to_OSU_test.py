import matplotlib.pyplot as plt
import numpy as np

# Read Data
# Experimental data digitized from Figure 68 in report
# Analytical (PyPonding) data generated from truss_model.py
exp_flat    = np.genfromtxt('OSU_data\Figure_68_flat.csv', delimiter=',', skip_header=1)
exp_pitched = np.genfromtxt('OSU_data\Figure_68_pitched.csv', delimiter=',', skip_header=1)
ana_flat    = np.genfromtxt('OSU_data\PyPonding_flat.csv', delimiter=',', skip_header=1)
ana_pitched = np.genfromtxt('OSU_data\PyPonding_pitched.csv', delimiter=',', skip_header=1)
ana_flat_linear    = np.genfromtxt('OSU_data\PyPonding_flat_linear.csv', delimiter=',', skip_header=1)
ana_pitched_linear = np.genfromtxt('OSU_data\PyPonding_pitched_linear.csv', delimiter=',', skip_header=1)

# Plot settings
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=8)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)

# Make plot
fig = plt.figure(figsize=(6.50,3.00))
ax = fig.add_axes([0.07,0.13,0.90,0.85])
line1, = plt.plot(exp_flat[:,0],exp_flat[:,1],'k-')
line3, = plt.plot(ana_flat_linear[:,0],ana_flat_linear[:,1],'g-')
line2, = plt.plot(ana_flat[:,0],ana_flat[:,1],'r-')
line4, = plt.plot(exp_pitched[:,0],exp_pitched[:,1],'k--')
line6, = plt.plot(ana_pitched_linear[:,0],ana_pitched_linear[:,1],'g--')
line5, = plt.plot(ana_pitched[:,0],ana_pitched[:,1],'r--')
plt.legend((line1, line2, line3, line4, line5, line6), 
    ('Experiment (flat)','PyPonding (flat)','PyPonding (flat,linear)',
     'Experiment (pitched)','PyPonding (pitched)','PyPonding (pitched,linear)'),
    frameon=False,loc=2)
plt.xlabel('Total water load (kips)')
plt.ylabel('Water level (in.)')
plt.xlim(0.0,50.0)
plt.ylim(0.00,18.00)
plt.savefig('Figure_11.png',dpi=300)
plt.savefig('Figure_11.pdf')

plt.show()