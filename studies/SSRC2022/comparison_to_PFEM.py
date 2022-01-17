import numpy as np
import matplotlib.pyplot as plt
from math import pi
from numpy import cos,cosh
from PyPonding.structures import ElasticBeam2d
from libdenavit.section.database import wide_flange_database
 
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
zj = 0*inch 
# PFEM results are available for zj = 0*inch and zj = 20*inch

if zj == 0*inch:
  normalized_volume1 = 20.0
  normalized_volume2 = 30.0
  normalized_volume3 = 40.0
  time = ""
elif zj == 20*inch:
  normalized_volume1 = 15.0
  normalized_volume2 = 30.0
  normalized_volume3 = 40.0
  time1 = "-30.030sec"
else:
  raise Exception(f'PFEM results not available for {zj =}')

VwM = (40.0*inch)*L*S
Vw1 = (normalized_volume1*inch)*L*S
Vw2 = (normalized_volume2*inch)*L*S
Vw3 = (normalized_volume3*inch)*L*S

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
fig = plt.figure(figsize=(3.25,2.50))
ax = fig.add_axes([0.15,0.15,0.82,0.82])
line2, = plt.plot(results_PyPonding_ramp_NP.water_volume/(L*S),results_PyPonding_ramp_NP.water_level,'k--')
line1, = plt.plot(results_PyPonding_ramp.water_volume/(L*S),results_PyPonding_ramp.water_level,'k-')
if zj == 0*inch:
  plt.plot(results_PyPonding_vol1.water_volume/(L*S),results_PyPonding_vol1.water_level,'ro')
  plt.plot(results_PyPonding_vol2.water_volume/(L*S),results_PyPonding_vol2.water_level,'go')
  plt.plot(results_PyPonding_vol3.water_volume/(L*S),results_PyPonding_vol3.water_level,'bo')
  plt.plot(normalized_volume1,10,'rx', markersize=5)
  plt.plot(normalized_volume2,20,'gx', markersize=5)
  plt.plot(normalized_volume3,25,'bx', markersize=5)
  line3, = plt.plot([],[],'o', color='grey', markersize=5)   
  line4, = plt.plot([],[],'x', color='grey', markersize=5) 
elif zj == 20*inch:
  line3, = plt.plot(results_PyPonding_vol1.water_volume/(L*S),results_PyPonding_vol1.water_level,'bo')
  line4, = plt.plot(normalized_volume1,17,'bx', markersize=5)
plt.legend((line1,line3,line4,line2), ('PyPonding', 'PyPonding', 'PFEM', 'No Ponding Effect'),frameon=False)
plt.xlabel('Normalized Water Volume, $V/LS$ (in.)')
plt.ylabel('Water Level (in.)')
if zj == 0*inch:
  plt.savefig('Figure_9a.png',dpi=300)
  plt.savefig('Figure_9a.pdf')
elif zj == 20*inch:
  plt.savefig('Figure_10a.png',dpi=300)
  plt.savefig('Figure_10a.pdf')

fig = plt.figure(figsize=(3.25,2.50))
ax = fig.add_axes([0.17,0.15,0.80,0.82])
#line1, = plt.plot(results_PyPonding_vol1_NP.position_along_length/ft,results_PyPonding_vol1_NP.deflection_along_length/inch,'r--')
#line2, = plt.plot(results_PyPonding_vol2_NP.position_along_length/ft,results_PyPonding_vol2_NP.deflection_along_length/inch,'g--')
#line3, = plt.plot(results_PyPonding_vol3_NP.position_along_length/ft,results_PyPonding_vol3_NP.deflection_along_length/inch,'b--')
if zj == 0*inch:
  line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.deflection_along_length/inch,'r-')
  line2, = plt.plot(results_PyPonding_vol2.position_along_length/ft,results_PyPonding_vol2.deflection_along_length/inch,'g-')
  line3, = plt.plot(results_PyPonding_vol3.position_along_length/ft,results_PyPonding_vol3.deflection_along_length/inch,'b-')
  disps1 = np.loadtxt(f"PFEM_data/disp/disp-{normalized_volume1}in{time}.txt")
  pfem1, = plt.plot(disps1[:,0],disps1[:,1],'r--',lw=3)
  disps2 = np.loadtxt(f"PFEM_data/disp/disp-{normalized_volume2}in{time}.txt")
  pfem2, = plt.plot(disps2[:,0],disps2[:,1],'g--',lw=3)
  disps3 = np.loadtxt(f"PFEM_data/disp/disp-{normalized_volume3}in{time}.txt")
  pfem3, = plt.plot(disps3[:,0],disps3[:,1],'b--',lw=3)
  plt.legend((line1, line2, line3, pfem1, pfem2, pfem3), 
    ("PyPonding, Volume 1", "PyPonding, Volume 2", "PyPonding, Volume 3", 
     "PFEM, Volume 1", "PFEM, Volume 2", "PFEM, Volume 3"),frameon=False)
elif zj == 20*inch:
  line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.deflection_along_length/inch,'b-')
  disps1 = np.loadtxt(f"PFEM_data/disp/disp-{normalized_volume1}in-{zj}in{time1}.txt")
  pfem1, = plt.plot(disps1[:,0],disps1[:,1],'b--',lw=3)
  plt.legend((line1, pfem1), ("PyPonding","PFEM"),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Deflection (in.)')
plt.xlim(0,L/ft)
if zj == 0*inch:
  plt.savefig('Figure_9b.png',dpi=300)
  plt.savefig('Figure_9b.pdf')
elif zj == 20*inch:
  plt.savefig('Figure_10b.png',dpi=300)
  plt.savefig('Figure_10b.pdf')

fig = plt.figure(figsize=(3.25,2.50))
ax = fig.add_axes([0.13,0.15,0.84,0.82])
#line1, = plt.plot(results_PyPonding_vol1_NP.position_along_length/ft,results_PyPonding_vol1_NP.axial_load_along_length/inch,'r--')
#line2, = plt.plot(results_PyPonding_vol2_NP.position_along_length/ft,results_PyPonding_vol2_NP.axial_load_along_length/inch,'g--')
#line3, = plt.plot(results_PyPonding_vol3_NP.position_along_length/ft,results_PyPonding_vol3_NP.axial_load_along_length/inch,'b--')
if zj == 0*inch:
  line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.axial_load_along_length/inch,'r-')
  line2, = plt.plot(results_PyPonding_vol2.position_along_length/ft,results_PyPonding_vol2.axial_load_along_length/inch,'g-')
  line3, = plt.plot(results_PyPonding_vol3.position_along_length/ft,results_PyPonding_vol3.axial_load_along_length/inch,'b-')
  axials1 = np.loadtxt(f"PFEM_data/axial/axial-{normalized_volume1}in{time}.txt")
  pfem1, = plt.plot(axials1[:,0],axials1[:,1],'r--',lw=3)
  axials2 = np.loadtxt(f"PFEM_data/axial/axial-{normalized_volume2}in{time}.txt")
  pfem2, = plt.plot(axials2[:,0],axials2[:,1],'g--',lw=3)
  axials3 = np.loadtxt(f"PFEM_data/axial/axial-{normalized_volume3}in{time}.txt")
  pfem3, = plt.plot(axials3[:,0],axials3[:,1],'b--',lw=3)
  plt.legend((line1, line2, line3, pfem1, pfem2, pfem3), 
    ("PyPonding, Volume 1", "PyPonding, Volume 2", "PyPonding, Volume 3", 
     "PFEM, Volume 1", "PFEM, Volume 2", "PFEM, Volume 3"),frameon=False,ncol=2)
  plt.ylim(-0.2,6)
elif zj == 20*inch:
  line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.axial_load_along_length/inch,'b-')
  axials1 = np.loadtxt(f"PFEM_data/axial/axial-{normalized_volume1}in-{zj}in{time1}.txt")
  pfem1, = plt.plot(axials1[:,0],axials1[:,1],'b--',lw=3)
  plt.yticks([0,1])
  plt.legend((line1, pfem1), ("PyPonding","PFEM"),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Axial Load (kips)')
plt.xlim(0,L/ft)
if zj == 0*inch:
  plt.savefig('Figure_9c.png',dpi=300)
  plt.savefig('Figure_9c.pdf')
elif zj == 20*inch:
  plt.savefig('Figure_10c.png',dpi=300)
  plt.savefig('Figure_10c.pdf')

fig = plt.figure(figsize=(3.25,2.50))
ax = fig.add_axes([0.17,0.15,0.80,0.82])
#line1, = plt.plot(results_PyPonding_vol1_NP.position_along_length/ft,results_PyPonding_vol1_NP.shear_along_length/inch,'r--')
#line2, = plt.plot(results_PyPonding_vol2_NP.position_along_length/ft,results_PyPonding_vol2_NP.shear_along_length/inch,'g--')
#line3, = plt.plot(results_PyPonding_vol3_NP.position_along_length/ft,results_PyPonding_vol3_NP.shear_along_length/inch,'b--')
if zj == 0*inch:
  line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.shear_along_length/inch,'r-')
  line2, = plt.plot(results_PyPonding_vol2.position_along_length/ft,results_PyPonding_vol2.shear_along_length/inch,'g-')
  line3, = plt.plot(results_PyPonding_vol3.position_along_length/ft,results_PyPonding_vol3.shear_along_length/inch,'b-')
  shears1 = np.loadtxt(f"PFEM_data/shear/shear-{normalized_volume1}in{time}.txt")
  pfem1, = plt.plot(shears1[:,0],shears1[:,1],'r--',lw=3)
  shears2 = np.loadtxt(f"PFEM_data/shear/shear-{normalized_volume2}in{time}.txt")
  pfem2, = plt.plot(shears2[:,0],shears2[:,1],'g--',lw=3)
  shears3 = np.loadtxt(f"PFEM_data/shear/shear-{normalized_volume3}in{time}.txt")
  pfem3, = plt.plot(shears3[:,0],shears3[:,1],'b--',lw=3)
  plt.legend((line1, line2, line3, pfem1, pfem2, pfem3), 
    ("PyPonding, Volume 1", "PyPonding, Volume 2", "PyPonding, Volume 3", 
     "PFEM, Volume 1", "PFEM, Volume 2", "PFEM, Volume 3"),frameon=False)
elif zj == 20*inch:
  line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.shear_along_length/inch,'b-')
  shears1 = np.loadtxt(f"PFEM_data/shear/shear-{normalized_volume1}in-{zj}in{time1}.txt")
  pfem1, = plt.plot(shears1[:,0],shears1[:,1],'b--',lw=3)
  plt.legend((line1, pfem1), ("PyPonding","PFEM"),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Shear (kips)')
plt.xlim(0,L/ft)
if zj == 0*inch:
  plt.savefig('Figure_9d.png',dpi=300)
  plt.savefig('Figure_9d.pdf')
elif zj == 20*inch:
  plt.savefig('Figure_10d.png',dpi=300)
  plt.savefig('Figure_10d.pdf')

fig = plt.figure(figsize=(3.25,2.50))
ax = fig.add_axes([0.17,0.15,0.80,0.82])
#line1, = plt.plot(results_PyPonding_vol1_NP.position_along_length/ft,results_PyPonding_vol1_NP.bending_moment_along_length/inch,'r--')
#line2, = plt.plot(results_PyPonding_vol2_NP.position_along_length/ft,results_PyPonding_vol2_NP.bending_moment_along_length/inch,'g--')
#line3, = plt.plot(results_PyPonding_vol3_NP.position_along_length/ft,results_PyPonding_vol3_NP.bending_moment_along_length/inch,'b--')
if zj == 0*inch:
  line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.bending_moment_along_length/inch,'r-')
  line2, = plt.plot(results_PyPonding_vol2.position_along_length/ft,results_PyPonding_vol2.bending_moment_along_length/inch,'g-')
  line3, = plt.plot(results_PyPonding_vol3.position_along_length/ft,results_PyPonding_vol3.bending_moment_along_length/inch,'b-')
  moments1 = np.loadtxt(f"PFEM_data/moment/moment-{normalized_volume1}in{time}.txt")
  pfem1, = plt.plot(moments1[:,0],moments1[:,1],'r--',lw=3)
  moments2 = np.loadtxt(f"PFEM_data/moment/moment-{normalized_volume2}in{time}.txt")
  pfem2, = plt.plot(moments2[:,0],moments2[:,1],'g--',lw=3)
  moments3 = np.loadtxt(f"PFEM_data/moment/moment-{normalized_volume3}in{time}.txt")
  pfem3, = plt.plot(moments3[:,0],moments3[:,1],'b--',lw=3)
  plt.legend((line1, line2, line3, pfem1, pfem2, pfem3), 
    ("PyPonding, Volume 1", "PyPonding, Volume 2", "PyPonding, Volume 3", 
     "PFEM, Volume 1", "PFEM, Volume 2", "PFEM, Volume 3"),frameon=False)
elif zj == 20*inch:
  line1, = plt.plot(results_PyPonding_vol1.position_along_length/ft,results_PyPonding_vol1.bending_moment_along_length/inch,'b-')
  moments1 = np.loadtxt(f"PFEM_data/moment/moment-{normalized_volume1}in-{zj}in{time1}.txt")
  pfem1, = plt.plot(moments1[:,0],moments1[:,1],'b--',lw=3)
  plt.legend((line1, pfem1), ("PyPonding","PFEM"),frameon=False)
plt.xlabel('Position Along Length of Beam (ft)')
plt.ylabel('Bending Moment (kip-in.)')
plt.xlim(0,L/ft)
if zj == 0*inch:
  plt.savefig('Figure_9e.png',dpi=300)
  plt.savefig('Figure_9e.pdf')
elif zj == 20*inch:
  plt.savefig('Figure_10e.png',dpi=300)
  plt.savefig('Figure_10e.pdf')


plt.show()