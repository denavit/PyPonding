import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from math import pi
from PyPonding.structures import wf_shapes

with open('StrengthEvaluationOutput.csv', 'r') as f:
    reader = csv.reader(f)
    your_list = list(reader)

num_data = len(your_list)-1

shape_name              = [None]*num_data
d_in                    = np.zeros(num_data)
Sx_in3                  = np.zeros(num_data)
Ix_in4                  = np.zeros(num_data)
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
    Sx_in3[i] = shape_data['Sx']
    Ix_in4[i] = shape_data['Ix']
    
    L_ft[i]                 = float(your_list[i+1][1])
    slope[i]                = float(your_list[i+1][2])
    dead_load_psf[i]        = float(your_list[i+1][3])
    zmax_inelastic[i]       = float(your_list[i+1][4])
    zmax_inelastic_noRS[i]  = float(your_list[i+1][5])
    zmax_elastic[i]         = float(your_list[i+1][6])
    tau[i]                  = float(your_list[i+1][7])
    zmax_elastic_reduced[i] = float(your_list[i+1][8])



# x-axis options

xlabel = 'Span-to-Depth Ratio ($L/d$)'
x = 12*L_ft/d_in
xlim = (15,32)
filename = 'StrengthEvaluationOutput.pdf'

# xlabel = '$L/S_x$ (1/in.$^2$)'
# x = 12*L_ft/Sx_in3
# xlim = (0,34)
# filename = 'StrengthEvaluationOutput_Sx.pdf'

# xlabel = 'Flexibility coefficient ($C$)'
# gamma = 0.0624/12**3
# S = 120
# E = 29000
# x = gamma*S*(12*L_ft)**4/(pi**4*E*Ix_in4)
# xlim = (0,0.28)
# filename = 'StrengthEvaluationOutput_C.pdf'

print(min(x))
print(max(x))


# Make Plot
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=8)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)

fig = plt.figure(figsize=(3.25,5.50))

ax1 = fig.add_axes([0.18,0.76,0.80,0.22])
ax1.plot([xlim[0],xlim[1]],[1,1],'k-',linewidth=0.5)
ax1.scatter(x,zmax_elastic/zmax_inelastic,4,color='k')
ax1.set(xlabel=xlabel, ylabel='Max Water Level Ratio\n(Nominal Stiffness)', xlim=xlim, ylim=(0.95,1.3))

ax2 = fig.add_axes([0.18,0.43,0.80,0.22])
ax2.scatter(x,tau,4,color='k')
ax2.set(xlabel=xlabel, ylabel='Computed Stiffness\nReduction Factor', xlim=xlim, ylim=(0.50,0.90))

ax3 = fig.add_axes([0.18,0.10,0.80,0.22])
ax3.plot([xlim[0],xlim[1]],[1,1],'k-',linewidth=0.5)
ax3.scatter(x,zmax_elastic_reduced/zmax_inelastic,4,color='k')
ax3.set(xlabel=xlabel, ylabel='Max Water Level Ratio\n(Reduced Stiffness)', xlim=xlim, ylim=(0.95,1.3))

fig.text(0.56,0.67,'(a)',fontsize=8)
fig.text(0.56,0.34,'(b)',fontsize=8)
fig.text(0.56,0.01,'(c)',fontsize=8)

plt.savefig(filename)
plt.show()
