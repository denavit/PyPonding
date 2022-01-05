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
zi = 0*inch
zj = 0*inch

VwM = (4.0*inch)*L*S
Vw1 = (1.0*inch)*L*S
Vw2 = (2.0*inch)*L*S
Vw3 = (3.0*inch)*L*S

######################################
## Run analysis with specific shape ##
######################################
shape_name = 'W14X22'

# Run OpenSees analysis
shape_data = wide_flange_database[shape_name.upper()]
I  = shape_data['Ix']*inch**4
beam = ElasticBeam2d(L,S,E,I,gamma,qD=qD,zi=zi,zj=zj)
results_PyPonding_ramp = beam.run_analysis_OPS('SimpleStepVolume',target_Vw=VwM,num_steps=500)
results_PyPonding_vol1 = beam.run_analysis_OPS('IterativeVolume', target_Vw=Vw1)
results_PyPonding_vol2 = beam.run_analysis_OPS('IterativeVolume', target_Vw=Vw2)
results_PyPonding_vol3 = beam.run_analysis_OPS('IterativeVolume', target_Vw=Vw3)
    
beam.include_ponding_effect = False
results_PyPonding_ramp_NP = beam.run_analysis_OPS('SimpleStepVolume',target_Vw=VwM,num_steps=100)
results_PyPonding_vol1_NP = beam.run_analysis_OPS('IterativeVolume', target_Vw=Vw1)
results_PyPonding_vol2_NP = beam.run_analysis_OPS('IterativeVolume', target_Vw=Vw2)
results_PyPonding_vol3_NP = beam.run_analysis_OPS('IterativeVolume', target_Vw=Vw3)

# Plot Results
fig = plt.figure(figsize=(3.50,2.50))
fig.add_axes([0.15,0.20,0.80,0.75])
line2, = plt.plot(results_PyPonding_ramp_NP.water_volume/(L*S),results_PyPonding_ramp_NP.water_level,'k--')
line1, = plt.plot(results_PyPonding_ramp.water_volume/(L*S),results_PyPonding_ramp.water_level,'k-')
plt.plot(results_PyPonding_vol1.water_volume/(L*S),results_PyPonding_vol1.water_level,'ro')
plt.plot(results_PyPonding_vol2.water_volume/(L*S),results_PyPonding_vol2.water_level,'go')
plt.plot(results_PyPonding_vol3.water_volume/(L*S),results_PyPonding_vol3.water_level,'bo')
plt.legend((line1,line2), ('PyPonding', 'No Ponding Effect'),frameon=False)
plt.xlabel('Normalized Water Volume, V/LS (in.)')
plt.ylabel('Water Level (in.)')
#plt.savefig('Figure_X1a.png',dpi=300)
#plt.savefig('Figure_X1a.pdf')

fig = plt.figure(figsize=(3.50,2.50))
fig.add_axes([0.15,0.20,0.80,0.75])
#line1, = plt.plot(results_PyPonding_vol1_NP.position_along_length/ft,results_PyPonding_vol1_NP.deflection_along_length/inch,'r--')
#line2, = plt.plot(results_PyPonding_vol2_NP.position_along_length/ft,results_PyPonding_vol2_NP.deflection_along_length/inch,'g--')
#line3, = plt.plot(results_PyPonding_vol3_NP.position_along_length/ft,results_PyPonding_vol3_NP.deflection_along_length/inch,'b--')
line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.deflection_along_length/inch,'r-')
line2, = plt.plot(results_PyPonding_vol2.position_along_length/ft,results_PyPonding_vol2.deflection_along_length/inch,'g-')
line3, = plt.plot(results_PyPonding_vol3.position_along_length/ft,results_PyPonding_vol3.deflection_along_length/inch,'b-')
#plt.legend((line1, line2), ('PyPonding', 'No Ponding Effect'),frameon=False)
plt.legend((line1, line2, line3), ('Volume 1', 'Volume 2', 'Volume 3'),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Displacement (in.)')
plt.xlim(0,L/ft)
#plt.ylim(-5,0)
#plt.savefig('Figure_X1a.png',dpi=300)
#plt.savefig('Figure_X1a.pdf')

fig = plt.figure(figsize=(3.50,2.50))
fig.add_axes([0.15,0.20,0.80,0.75])
#line1, = plt.plot(results_PyPonding_vol1_NP.position_along_length/ft,results_PyPonding_vol1_NP.axial_load_along_length/inch,'r--')
#line2, = plt.plot(results_PyPonding_vol2_NP.position_along_length/ft,results_PyPonding_vol2_NP.axial_load_along_length/inch,'g--')
#line3, = plt.plot(results_PyPonding_vol3_NP.position_along_length/ft,results_PyPonding_vol3_NP.axial_load_along_length/inch,'b--')
line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.axial_load_along_length/inch,'r-')
line2, = plt.plot(results_PyPonding_vol2.position_along_length/ft,results_PyPonding_vol2.axial_load_along_length/inch,'g-')
line3, = plt.plot(results_PyPonding_vol3.position_along_length/ft,results_PyPonding_vol3.axial_load_along_length/inch,'b-')
#plt.legend((line1, line2), ('PyPonding', 'No Ponding Effect'),frameon=False)
plt.legend((line1, line2, line3), ('Volume 1', 'Volume 2', 'Volume 3'),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Axial Load (kips)')
plt.xlim(0,L/ft)
#plt.ylim(-5,0)
#plt.savefig('Figure_X1a.png',dpi=300)
#plt.savefig('Figure_X1a.pdf')

fig = plt.figure(figsize=(3.50,2.50))
fig.add_axes([0.15,0.20,0.80,0.75])
#line1, = plt.plot(results_PyPonding_vol1_NP.position_along_length/ft,results_PyPonding_vol1_NP.shear_along_length/inch,'r--')
#line2, = plt.plot(results_PyPonding_vol2_NP.position_along_length/ft,results_PyPonding_vol2_NP.shear_along_length/inch,'g--')
#line3, = plt.plot(results_PyPonding_vol3_NP.position_along_length/ft,results_PyPonding_vol3_NP.shear_along_length/inch,'b--')
line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.shear_along_length/inch,'r-')
line2, = plt.plot(results_PyPonding_vol2.position_along_length/ft,results_PyPonding_vol2.shear_along_length/inch,'g-')
line3, = plt.plot(results_PyPonding_vol3.position_along_length/ft,results_PyPonding_vol3.shear_along_length/inch,'b-')
#plt.legend((line1, line2), ('PyPonding', 'No Ponding Effect'),frameon=False)
plt.legend((line1, line2, line3), ('Volume 1', 'Volume 2', 'Volume 3'),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Shear (kips)')
plt.xlim(0,L/ft)
#plt.ylim(-5,0)
#plt.savefig('Figure_X1a.png',dpi=300)
#plt.savefig('Figure_X1a.pdf')

fig = plt.figure(figsize=(3.50,2.50))
fig.add_axes([0.15,0.20,0.80,0.75])
#line1, = plt.plot(results_PyPonding_vol1_NP.position_along_length/ft,results_PyPonding_vol1_NP.bending_moment_along_length/inch,'r--')
#line2, = plt.plot(results_PyPonding_vol2_NP.position_along_length/ft,results_PyPonding_vol2_NP.bending_moment_along_length/inch,'g--')
#line3, = plt.plot(results_PyPonding_vol3_NP.position_along_length/ft,results_PyPonding_vol3_NP.bending_moment_along_length/inch,'b--')
line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.bending_moment_along_length/inch,'r-')
line2, = plt.plot(results_PyPonding_vol2.position_along_length/ft,results_PyPonding_vol2.bending_moment_along_length/inch,'g-')
line3, = plt.plot(results_PyPonding_vol3.position_along_length/ft,results_PyPonding_vol3.bending_moment_along_length/inch,'b-')
#plt.legend((line1, line2), ('PyPonding', 'No Ponding Effect'),frameon=False)
plt.legend((line1, line2, line3), ('Volume 1', 'Volume 2', 'Volume 3'),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Bending Moment (kip-in.)')
plt.xlim(0,L/ft)
#plt.ylim(-5,0)
#plt.savefig('Figure_X1a.png',dpi=300)
#plt.savefig('Figure_X1a.pdf')

plt.show()