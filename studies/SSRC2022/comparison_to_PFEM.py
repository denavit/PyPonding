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
qD = 20.0*psf
gamma = 62.4*pcf
zi = 0*inch
zj = 20*inch

if zj <= 0:
  height1 = 20.0
  height2 = 30.0
  height3 = 40.0
else:
  height1 = 15.0
  height2 = 20.0
  height3 = 25.0

VwM = (40.0*inch)*L*S
Vw1 = (height1*inch)*L*S
Vw2 = (height2*inch)*L*S
Vw3 = (height3*inch)*L*S

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
if zj <= 0:
  plt.plot(height1,10,'rx', markersize=5)
  plt.plot(height2,20,'gx', markersize=5)
  plt.plot(height3,25,'bx', markersize=5)
else:
  plt.title(f'Sloped Beam: zj = {zj} in')
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
if zj <= 0:
  disps1 = np.loadtxt(f"disp/disp-{height1}in.txt")
  pfem1, = plt.plot(disps1[:,0],disps1[:,1],'r--',lw=3)
  disps2 = np.loadtxt(f"disp/disp-{height2}in.txt")
  pfem2, = plt.plot(disps2[:,0],disps2[:,1],'g--',lw=3)
  disps3 = np.loadtxt(f"disp/disp-{height3}in.txt")
  pfem3, = plt.plot(disps3[:,0],disps3[:,1],'b--',lw=3)
else:
  disps1 = np.loadtxt(f"disp/disp-{height1}in-{zj}in.txt")
  pfem1, = plt.plot(disps1[:,0],disps1[:,1],'r--',lw=3)
  disps2 = np.loadtxt(f"disp/disp-{height2}in-{zj}in.txt")
  pfem2, = plt.plot(disps2[:,0],disps2[:,1],'g--',lw=3)
  disps3 = np.loadtxt(f"disp/disp-{height3}in-{zj}in.txt")
  pfem3, = plt.plot(disps3[:,0],disps3[:,1],'b--',lw=3)
  plt.title(f'Sloped Beam: zj = {zj} in')
#plt.legend((line1, line2), ('PyPonding', 'No Ponding Effect'),frameon=False)
plt.legend((line1, line2, line3, pfem1, pfem2, pfem3), ('Volume 1', 'Volume 2', 'Volume 3', 'PFEM 1', 'PFEM 2', 'PFEM 3'),frameon=False)
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
if zj <= 0:
  axials1 = np.loadtxt(f"axial/axial-{height1}in.txt")
  pfem1, = plt.plot(axials1[:,0],axials1[:,1],'r--',lw=3)
  axials2 = np.loadtxt(f"axial/axial-{height2}in.txt")
  pfem2, = plt.plot(axials2[:,0],axials2[:,1],'g--',lw=3)
  axials3 = np.loadtxt(f"axial/axial-{height3}in.txt")
  pfem3, = plt.plot(axials3[:,0],axials3[:,1],'b--',lw=3)
else:
  axials1 = np.loadtxt(f"axial/axial-{height1}in-{zj}in.txt")
  pfem1, = plt.plot(axials1[:,0],axials1[:,1],'r--',lw=3)
  axials2 = np.loadtxt(f"axial/axial-{height2}in-{zj}in.txt")
  pfem2, = plt.plot(axials2[:,0],axials2[:,1],'g--',lw=3)
  axials3 = np.loadtxt(f"axial/axial-{height3}in-{zj}in.txt")
  pfem3, = plt.plot(axials3[:,0],axials3[:,1],'b--',lw=3)
  plt.title(f'Sloped Beam: zj = {zj} in')
plt.legend((line1, line2, line3, pfem1, pfem2, pfem3), ('Volume 1', 'Volume 2', 'Volume 3', 'PFEM 1', 'PFEM 2', 'PFEM 3'),frameon=False)
#plt.legend((line1, line2), ('PyPonding', 'No Ponding Effect'),frameon=False)
# plt.legend((line1, line2, line3), ('Volume 1', 'Volume 2', 'Volume 3'),frameon=False)
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
if zj <= 0:
  shears1 = np.loadtxt(f"shear/shear-{height1}in.txt")
  pfem1, = plt.plot(shears1[:,0],shears1[:,1],'r--',lw=3)
  shears2 = np.loadtxt(f"shear/shear-{height2}in.txt")
  pfem2, = plt.plot(shears2[:,0],shears2[:,1],'g--',lw=3)
  shears3 = np.loadtxt(f"shear/shear-{height3}in.txt")
  pfem3, = plt.plot(shears3[:,0],shears3[:,1],'b--',lw=3)
else:
  shears1 = np.loadtxt(f"shear/shear-{height1}in-{zj}in.txt")
  pfem1, = plt.plot(shears1[:,0],shears1[:,1],'r--',lw=3)
  shears2 = np.loadtxt(f"shear/shear-{height2}in-{zj}in.txt")
  pfem2, = plt.plot(shears2[:,0],shears2[:,1],'g--',lw=3)
  shears3 = np.loadtxt(f"shear/shear-{height3}in-{zj}in.txt")
  pfem3, = plt.plot(shears3[:,0],shears3[:,1],'b--',lw=3)
  plt.title(f'Sloped Beam: zj = {zj} in')

plt.legend((line1, line2, line3, pfem1, pfem2, pfem3), ('Volume 1', 'Volume 2', 'Volume 3', 'PFEM 1', 'PFEM 2', 'PFEM 3'),frameon=False)
#plt.legend((line1, line2), ('PyPonding', 'No Ponding Effect'),frameon=False)
# plt.legend((line1, line2, line3), ('Volume 1', 'Volume 2', 'Volume 3'),frameon=False)
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
if zj <= 0:
  moments1 = np.loadtxt(f"moment/moment-{height1}in.txt")
  pfem1, = plt.plot(moments1[:,0],moments1[:,1],'r--',lw=3)
  moments2 = np.loadtxt(f"moment/moment-{height2}in.txt")
  pfem2, = plt.plot(moments2[:,0],moments2[:,1],'g--',lw=3)
  moments3 = np.loadtxt(f"moment/moment-{height3}in.txt")
  pfem3, = plt.plot(moments3[:,0],moments3[:,1],'b--',lw=3)
else:
  moments1 = np.loadtxt(f"moment/moment-{height1}in-{zj}in.txt")
  pfem1, = plt.plot(moments1[:,0],moments1[:,1],'r--',lw=3)
  moments2 = np.loadtxt(f"moment/moment-{height2}in-{zj}in.txt")
  pfem2, = plt.plot(moments2[:,0],moments2[:,1],'g--',lw=3)
  moments3 = np.loadtxt(f"moment/moment-{height3}in-{zj}in.txt")
  pfem3, = plt.plot(moments3[:,0],moments3[:,1],'b--',lw=3)
  plt.title(f'Sloped Beam: zj = {zj} in')

plt.legend((line1, line2, line3, pfem1, pfem2, pfem3), ('Volume 1', 'Volume 2', 'Volume 3', 'PFEM 1', 'PFEM 2', 'PFEM 3'),frameon=False)
#plt.legend((line1, line2), ('PyPonding', 'No Ponding Effect'),frameon=False)
# plt.legend((line1, line2, line3), ('Volume 1', 'Volume 2', 'Volume 3'),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Bending Moment (kip-in.)')
plt.xlim(0,L/ft)
#plt.ylim(-5,0)
#plt.savefig('Figure_X1a.png',dpi=300)
#plt.savefig('Figure_X1a.pdf')

plt.show()