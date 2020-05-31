import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from wide_flange import wf_shapes

with open('StrengthEvaluationOutput.csv', 'r') as f:
    reader = csv.reader(f)
    your_list = list(reader)

num_data = len(your_list)-1

shape_name              = [None]*num_data
d_in                    = np.zeros(num_data)
L_ft                    = np.zeros(num_data)
slope                   = np.zeros(num_data)
dead_load_psf           = np.zeros(num_data)
zmax_inelastic          = np.zeros(num_data)
zmax_inelastic_noRS     = np.zeros(num_data)
zmax_elastic            = np.zeros(num_data)
tau                     = np.zeros(num_data)
zmax_elastic_reduced    = np.zeros(num_data)

for i in range(num_data):
    shape_name[i] = your_list[i+1][0]
    shape_data = wf_shapes[shape_name[i]]
    d_in[i] = shape_data['d']
    
    L_ft[i]                 = float(your_list[i+1][1])
    slope[i]                = float(your_list[i+1][2])
    dead_load_psf[i]        = float(your_list[i+1][3])
    zmax_inelastic[i]       = float(your_list[i+1][4])
    zmax_inelastic_noRS[i]  = float(your_list[i+1][5])
    zmax_elastic[i]         = float(your_list[i+1][6])
    tau[i]                  = float(your_list[i+1][7])
    zmax_elastic_reduced[i] = float(your_list[i+1][8])


plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=8)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)

fig = plt.figure(figsize=(3.25,5.50))

ax1 = fig.add_axes([0.18,0.76,0.80,0.22])
ax1.plot([15,32],[1,1],'k-',linewidth=0.5)
ax1.scatter(12*L_ft/d_in,zmax_elastic/zmax_inelastic,4,color='k')
ax1.set(xlabel='Span-to-Depth Ratio ($L/d$)', ylabel='Max Water Level Ratio\n(Nominal Stiffness)', xlim=(15,32), ylim=(0.95,1.3))

ax2 = fig.add_axes([0.18,0.43,0.80,0.22])
ax2.scatter(12*L_ft/d_in,tau,4,color='k')
ax2.set(xlabel='Span-to-Depth Ratio ($L/d$)', ylabel='Computed Stiffness\nReduction Factor', xlim=(15,32), ylim=(0.50,0.90))

ax3 = fig.add_axes([0.18,0.10,0.80,0.22])
ax3.plot([15,32],[1,1],'k-',linewidth=0.5)
ax3.scatter(12*L_ft/d_in,zmax_elastic_reduced/zmax_inelastic,4,color='k')
ax3.set(xlabel='Span-to-Depth Ratio ($L/d$)', ylabel='Max Water Level Ratio\n(Reduced Stiffness)', xlim=(15,32), ylim=(0.95,1.3))

fig.text(0.56,0.67,'(a)',fontsize=8)
fig.text(0.56,0.34,'(b)',fontsize=8)
fig.text(0.56,0.01,'(c)',fontsize=8)

plt.savefig('StrengthEvaluationOutput.png',dpi=300)
plt.savefig('StrengthEvaluationOutput.pdf')
plt.show()
