import numpy as np
import matplotlib.pyplot as plt
from math import pi
from PyPonding.structures import IdealizedBay
 
# Plot settings
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=6)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
 
# Define units
inch = 1.0
kip = 1.0
lb  = kip/1000.0
ft  = 12.0*inch
in_per_ft = inch/ft
ksi = kip/inch**2
plf = lb/ft
psf = lb/ft**2
pcf = lb/ft**3
kipft = kip*ft

# Input parameters
Lp = 40*ft
Ls = 40*ft
num_spaces = 8
E  = 29000.0*ksi
zw = 2*inch
qD = 20.0*psf
gamma = 62.4*pcf

S  = Lp/num_spaces

##########################################################
## Run analysis with different flexibility coefficients ##
##########################################################
Cp_list = np.array([0.1,0.2,0.3])

# Run OpenSees analysis
Cs_list_OPS = np.array([0.01,0.1,0.2,0.3,0.4])
#Cs_list_OPS = np.array([0.01])
Mp_max_list_OPS = np.zeros((len(Cp_list),len(Cs_list_OPS)))
Ms_max_list_OPS = np.zeros((len(Cp_list),len(Cs_list_OPS)))

for i in range(len(Cp_list)):
    Cp = Cp_list[i]
    Ip = (gamma*Ls*Lp**4)/(pi**4*E*Cp)

    for j in range(len(Cs_list_OPS)):
        Cs = Cs_list_OPS[j]
        Is = (gamma*S*Ls**4)/(pi**4*E*Cs)

        print(f'\n{Cp = },{Cs = }')
        input = {
            'primary_member_span': Lp,
            'secondary_member_span': Ls,
            'number_of_joist_spaces': num_spaces,
            'dead_load_uniform': qD,
            'dead_load_primary_member': 0*plf,
            'water_density': 62.4*pcf,  
            'snow_density': 0*pcf,
            'snow_height': 0*inch,
            'alpha': 1.0,
            'load_factor_dead':    1.0,
            'load_factor_ponding': 1.0,
            'load_factor_snow1':   0.0,
            'load_factor_snow2':   0.0,
            'z_TL': 0*inch,
            'z_TR': 0*inch,
            'z_BL': 0*inch,
            'z_BR': 0*inch,
            'secondary_member_camber': 0*inch,
            'primary_member_camber_T': 0*inch,
            'primary_member_camber_B': 0*inch,
            'edge_condition_L': 'mirrored',
            'edge_condition_R': 'mirrored',
            'edge_condition_T': 'mirrored',
            'edge_condition_B': 'mirrored',
            'E': E,
            'As': 100*inch**2,
            'Ap': 100*inch**2,
            'Is': Is,
            'Ip': Ip,
            'analsis_engine': 'OpenSees',
        }

        bay = IdealizedBay(**input)
        results = bay.Run_Analysis(zw)
        Mp_max_list_OPS[i,j] = np.amax(results.bot_primary_member_moment)
        Ms_max_list_OPS[i,j] = np.amax(results.secondary_members_moment)


# Compute closed-form results (Marino 1966)
Cs_list_Marino = np.linspace(0.01,0.40,40)
Mp_max_list_Marino = np.zeros((len(Cp_list),len(Cs_list_Marino)))
Ms_max_list_Marino = np.zeros((len(Cp_list),len(Cs_list_Marino)))

wp = (qD + gamma*zw)*Ls
Mp_max_no_ponding = wp*Lp**2/8

ws = (qD + gamma*zw)*S
Ms_max_no_ponding = ws*Ls**2/8

for i in range(len(Cp_list)):
    Cp = Cp_list[i]
    αp = Cp/(1-Cp)

    for j in range(len(Cs_list_Marino)):
        Cs = Cs_list_Marino[j]
        αs = Cs/(1-Cs)
        ρ = Cs/Cp

        Mp_max_list_Marino[i,j] = Mp_max_no_ponding*(1+αp*(1+0.25*pi*αs+0.25*pi*ρ*(1+αs))/(1-0.25*pi*αp*αs))
        Ms_max_list_Marino[i,j] = Ms_max_no_ponding*(1+αs*(1+0.03125*pi**3*αp+0.125*pi**2*(1+αp)/ρ+0.185*αp*αs)/(1-0.25*pi*αp*αs))

# Plot Results
fig = plt.figure(figsize=(3.25,2.50))
ax = fig.add_axes([0.19,0.15,0.78,0.82])
line1, = plt.plot(Cs_list_Marino,Mp_max_list_Marino[0,:]/kipft,'r-')
line2, = plt.plot(Cs_list_Marino,Mp_max_list_Marino[1,:]/kipft,'g-')
line3, = plt.plot(Cs_list_Marino,Mp_max_list_Marino[2,:]/kipft,'b-')
line4, = plt.plot(Cs_list_OPS,Mp_max_list_OPS[0,:]/kipft,'ro')
line5, = plt.plot(Cs_list_OPS,Mp_max_list_OPS[1,:]/kipft,'go')
line6, = plt.plot(Cs_list_OPS,Mp_max_list_OPS[2,:]/kipft,'bo')
line7, = plt.plot([0.00,0.41],[Mp_max_no_ponding/kipft]*2,'k--')
plt.legend((line1, line2, line3, line4, line5, line6, line7), 
    ('Marino, $C_p = 0.1$', 'Marino, $C_p = 0.2$', 'Marino, $C_p = 0.3$', 
    'PyPonding, $C_p = 0.1$', 'PyPonding, $C_p = 0.2$', 'PyPonding, $C_p = 0.3$','No Ponding Effect'),
    frameon=False,loc=2)
plt.xlabel('Flexibility coefficient of secondary member, $C_s$')
plt.ylabel('Bending moment at mid-span\nof primary member (kip-ft)')
plt.xlim(0.00,0.41)
plt.savefig('Figure_7a.png',dpi=300)
plt.savefig('Figure_7a.pdf')

fig = plt.figure(figsize=(3.25,2.50))
ax = fig.add_axes([0.19,0.15,0.78,0.82])
line1, = plt.plot(Cs_list_Marino,Ms_max_list_Marino[0,:]/kipft,'r-')
line2, = plt.plot(Cs_list_Marino,Ms_max_list_Marino[1,:]/kipft,'g-')
line3, = plt.plot(Cs_list_Marino,Ms_max_list_Marino[2,:]/kipft,'b-')
line4, = plt.plot(Cs_list_OPS,Ms_max_list_OPS[0,:]/kipft,'ro')
line5, = plt.plot(Cs_list_OPS,Ms_max_list_OPS[1,:]/kipft,'go')
line6, = plt.plot(Cs_list_OPS,Ms_max_list_OPS[2,:]/kipft,'bo')
line7, = plt.plot([0.00,0.41],[Ms_max_no_ponding/kipft]*2,'k--')
plt.legend((line1, line2, line3, line4, line5, line6, line7), 
    ('Marino, $C_p = 0.1$', 'Marino, $C_p = 0.2$', 'Marino, $C_p = 0.3$', 
    'PyPonding, $C_p = 0.1$', 'PyPonding, $C_p = 0.2$', 'PyPonding, $C_p = 0.3$','No Ponding Effect'),
    frameon=False,loc=2)
plt.xlabel('Flexibility coefficient of secondary member, $C_s$')
plt.ylabel('Bending moment at mid-span\nof secondary member (kip-ft)')
plt.xlim(0.00,0.41)
plt.savefig('Figure_7b.png',dpi=300)
plt.savefig('Figure_7b.pdf')

plt.show()