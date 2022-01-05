import numpy as np
import matplotlib.pyplot as plt
from math import pi
from numpy import cos,cosh
from PyPonding.structures import ElasticBeam2d
from libdenavit.section.database import wide_flange_database
 
# Plot settings
plt.rc('text',usetex=False)
plt.rc('font',family='sans-serif')
plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=8)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
 
# Define units
inch = 1.0
kip = 1.0
lb  = kip/1000.0
ft  = 12.0*inch
in_per_ft = inch/ft
ksi = kip/inch**2
psf = lb/ft**2
pcf = psf/ft
kipft = kip*ft

# Input parameters
L  = 40*ft
S  = 10*ft
E  = 29000.0*ksi
zw = 2*inch
qD = 10.0*psf
gamma = 62.4*pcf

######################################
## Run analysis with specific shape ##
######################################
shape_name = 'W14X22'

# Run OpenSees analysis
shape_data = wide_flange_database[shape_name.upper()]
I  = shape_data['Ix']*inch**4
beam = ElasticBeam2d(L,S,E,I,gamma,qD=qD)
results_OPS = beam.run_analysis_OPS('IterativeLevel',target_zw=zw)

# Compute closed-form results (Silver 2010)
C  = beam.C()
wi = (qD + gamma*zw)*S
x  = np.linspace(0.0, L, num=100)
y_no_ponding = -wi*x*(L**3-2*L*x**2+x**3)/(24*E*I)
M_no_ponding = 0.5*wi*x*(L-x)
fact = pi*C**0.25
y_Silver = -wi/(2*gamma*S)*( (cos(fact*(0.5-x/L))/cos(fact*0.5) - 1) - (1 - cosh(fact*(0.5-x/L))/cosh(fact*0.5)) )
M_Silver = (wi*L**2)/(2*pi**2*C**0.5)*( cos(fact*(0.5-x/L))/cos(0.5*fact) - cosh(fact*(0.5-x/L))/cosh(0.5*fact) )

# Plot Results
fig = plt.figure(figsize=(3.50,2.50))
ax = fig.add_axes([0.15,0.21,0.80,0.75])
line2, = plt.plot(x/ft,y_Silver/inch,'k-')
line3, = plt.plot(x/ft,y_no_ponding/inch,'k--')
line1, = plt.plot(results_OPS.position_along_length/ft,results_OPS.deflection_along_length/inch,'ko')
plt.legend((line1, line2, line3), ('Numerical Results', 'Closed-form Solution', 'No Ponding Effect'),frameon=False)
#plt.legend((line2, line3), ('Including Ponding Effect', 'No Ponding Effect'),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Displacement (in.)')
plt.xlim(0,L/ft)
plt.ylim(-5,0)
#plt.savefig('Figure_X1a.png',dpi=300)
#plt.savefig('Figure_X1a.pdf')

fig = plt.figure(figsize=(3.50,2.50))
ax = fig.add_axes([0.15,0.21,0.80,0.75])
line2, = plt.plot(x/ft,M_Silver/kipft,'k-')
line3, = plt.plot(x/ft,M_no_ponding/kipft,'k--')
line1, = plt.plot(results_OPS.position_along_length/ft,results_OPS.bending_moment_along_length/kipft,'ko')
plt.legend((line1, line2, line3), ('Numerical Results', 'Closed-form Solution', 'No Ponding Effect'),frameon=False)
#plt.legend((line2, line3), ('Including Ponding Effect', 'No Ponding Effect'),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Bending Moment (kip-ft)')
plt.xlim(0,L/ft)
plt.ylim(0,100)
#plt.savefig('Figure_X1b.png',dpi=300)
#plt.savefig('Figure_X1b.pdf')


##########################################################
## Run analysis with different flexibility coefficients ##
##########################################################

# Run OpenSees analysis
C_list_OPS = np.array([0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
Mmax_list_OPS = np.zeros_like(C_list_OPS)
ymax_list_OPS = np.zeros_like(C_list_OPS)

for i in range(len(C_list_OPS)):
    C = C_list_OPS[i]
    I = (gamma*S*L**4)/(pi**4*E*C)
    
    beam = ElasticBeam2d(L,S,E,I,gamma,qD=qD)
    beam.maximum_number_of_iterations_level = 100
    results_OPS = beam.run_analysis_OPS('IterativeLevel',target_zw=zw)
    Mmax_list_OPS[i] = max(results_OPS.bending_moment_along_length)
    ymax_list_OPS[i] = min(results_OPS.deflection_along_length)

# Compute closed-form results (Silver 2010)
C_list = np.linspace(0.01,0.80,80)
Mmax_list_Silver = np.zeros_like(C_list)
ymax_list_Silver = np.zeros_like(C_list)
Mmax_list_no_ponding = np.zeros_like(C_list)
ymax_list_no_ponding = np.zeros_like(C_list)

for i in range(len(C_list)):
    C = C_list[i]
    I = (gamma*S*L**4)/(pi**4*E*C)   
    
    wi = (qD + gamma*zw)*S
    
    Mmax_list_no_ponding[i] = wi*L**2/8
    ymax_list_no_ponding[i] = -(5*wi*L**4)/(384*E*I)
    
    x  = L/2
    fact = pi*C**0.25
    ymax_list_Silver[i] = -wi/(2*gamma*S)*( (cos(fact*(0.5-x/L))/cos(fact*0.5) - 1) - (1 - cosh(fact*(0.5-x/L))/cosh(fact*0.5)) )
    Mmax_list_Silver[i] = (wi*L**2)/(2*pi**2*C**0.5)*( cos(fact*(0.5-x/L))/cos(0.5*fact) - cosh(fact*(0.5-x/L))/cosh(0.5*fact) )

# Plot Results
fig = plt.figure(figsize=(3.50,2.50))
ax = fig.add_axes([0.15,0.21,0.80,0.75])
line2, = plt.plot(C_list,ymax_list_Silver/inch,'k-')
line3, = plt.plot(C_list,ymax_list_no_ponding/inch,'k--')
line1, = plt.plot(C_list_OPS,ymax_list_OPS/inch,'ko')
plt.legend((line1, line2, line3), ('Numerical Results', 'Closed-form Solution', 'No Ponding Effect'),frameon=False)
plt.xlabel('Flexibility coefficient, C')
plt.ylabel('Displacement at mid-span (in.)')
plt.xlim(0.00,0.81)
#plt.ylim(-5,0)
#plt.savefig('Figure_X2a.png',dpi=300)
#plt.savefig('Figure_X2a.pdf')

fig = plt.figure(figsize=(3.50,2.50))
ax = fig.add_axes([0.15,0.21,0.80,0.75])
line2, = plt.plot(C_list,Mmax_list_Silver/kipft,'k-')
line3, = plt.plot(C_list,Mmax_list_no_ponding/kipft,'k--')
line1, = plt.plot(C_list_OPS,Mmax_list_OPS/kipft,'ko')
plt.legend((line1, line2, line3), ('Numerical Results', 'Closed-form Solution', 'No Ponding Effect'),frameon=False)
plt.xlabel('Flexibility coefficient, C')
plt.ylabel('Bending moment at mid-span (kip-ft)')
plt.xlim(0.00,0.81)
#plt.ylim(0,100)
#plt.savefig('Figure_X2b.png',dpi=300)
#plt.savefig('Figure_X2b.pdf')

plt.show()