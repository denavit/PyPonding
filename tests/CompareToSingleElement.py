import openseespy.opensees as ops
from math import pi,sin,cos,ceil,atan2
import numpy as np
import matplotlib.pyplot as plt
import PyPonding.FE as FE
from PyPonding import PondingLoadManager2d
from PyPonding.structures import steel_beam


# Define units
inch = 1.0
kip = 1.0

lb = kip/1000.0
ft = 12.0*inch

in_per_ft = inch/ft

ksi = kip/inch**2
psf = lb/ft**2
pcf = psf/ft

# Input parameters
E       = 29000.0*ksi
A       = 100.*inch**2
Iz      = 100.*inch**4
TW      = 5*ft
L       = 40*ft
slope   = 0.25*in_per_ft
qD      = 10.0*psf
gamma   = 62.4*pcf

# Computed parameters
zi      = 0.0*inch
zj      = L*slope
wD      = qD*TW # Dead load per length

# Set OpenSees analysis options
max_volume  = (20*inch)*L*TW
num_steps   = 1000
vol_tol     = max_volume/num_steps/10000.
nsteps_vol  = 30
percent_drop = 0.05

# ---------------------------------------------------------------------- #
#   Run first OpenSees analysis (with multiple elements along length)    #
# ---------------------------------------------------------------------- #
num_elements = 20

# set modelbuilder
ops.model('basic', '-ndm', 2, '-ndf', 3)

# create nodes
for i in range(num_elements+1):
    x_over_L = float(i)/(num_elements)
    ix = L*x_over_L
    iy = zi+(zj-zi)*x_over_L
    ops.node(i,ix,iy)

# set boundary condition
ops.fix(0,1,1,0)
ops.fix(num_elements,0,1,0)

# define coordinate transformation
ops.geomTransf('Linear',1)

# define cross section
ops.section('Elastic', 1, E, A, Iz) 
ops.beamIntegration('Lobatto', 1, 1, 3)

# Time series for loads
ops.timeSeries("Constant", 1)

# define elements
for i in range(0,num_elements):
    ops.element("dispBeamColumn",i,i,i+1,1,1)            

# Dead load
ops.pattern('Plain',-1,1)
for i in range(0,num_elements):
    iNodes = ops.eleNodes(i)
    xI,yI = ops.nodeCoord(iNodes[0])
    xJ,yJ = ops.nodeCoord(iNodes[1])
    ele_angle = atan2(yJ-yI, xJ-xI)
    ops.eleLoad('-ele',i,'-type','beamUniform',-wD*cos(ele_angle),-wD*sin(ele_angle))
    
# ------------------------------
# Start of analysis generation
# ------------------------------

# create SOE
ops.system("UmfPack")

# create DOF number
ops.numberer("RCM")

# create constraint handler
ops.constraints("Plain")

# create integrator
ops.integrator("LoadControl", 0.0)

# create algorithm
ops.algorithm("Newton")

# create analysis object
ops.analysis("Static")

# Run dead load analysis
ops.analyze(1)
zo = float('inf')
for i in range(num_elements+1):
    iz = ops.nodeCoord(i, 2) + ops.nodeDisp(i, 2)
    if iz < zo:
        zo = iz

# Initilize data        
data_volume_1 = np.zeros((num_steps+1,1))
data_height_1 = np.full((num_steps+1,1),zo)
data_reactIx_1 = np.zeros((num_steps+1,1))
data_reactIy_1 = np.zeros((num_steps+1,1))
data_reactJy_1 = np.zeros((num_steps+1,1))
end_step = num_steps+1

ops.reactions()
data_reactIx_1[0] = ops.nodeReaction(0,1)
data_reactIy_1[0] = ops.nodeReaction(0,2)
data_reactJy_1[0] = ops.nodeReaction(num_elements,2)

# ------------------------------
# Finally perform the analysis
# ------------------------------

# define ponding load cells    
PondingLoadManager = PondingLoadManager2d()
for i in range(0,num_elements):
    PondingLoadManager.add_cell(i,i,i+1,gamma,TW)

# Perform analysis, ramping up volume      
zw = zo+0.1

for iStep in range(0,num_steps):
    
    target_volume = (iStep+1)/num_steps*max_volume
    
    # Update ponding load cells
    PondingLoadManager.update()
    
    # Estimate water height
    for i in range(nsteps_vol):
        (V,dVdz) = PondingLoadManager.get_volume(zw)
        zw = zw - (V-target_volume)/dVdz
        if abs(target_volume-V) <= vol_tol:
            break 
    
    # Compute load vector
    PondingLoadManager.compute_current_load_vector(zw)
        
    # Apply difference to model
    ops.pattern("Plain", iStep, 1)
    PondingLoadManager.apply_load_increment()
    PondingLoadManager.commit_current_load_vector()

    # Run analysis
    ops.analyze(1)
    ops.reactions()
    
    # Store Data
    data_volume_1[iStep+1] = target_volume
    data_height_1[iStep+1] = zw
    data_reactIx_1[iStep+1] = ops.nodeReaction(0,1)
    data_reactIy_1[iStep+1] = ops.nodeReaction(0,2)
    data_reactJy_1[iStep+1] = ops.nodeReaction(num_elements,2)

    # Stop analysis if water level too low
    if (zw-zo) <= (1-percent_drop)*(np.amax(data_height_1)-zo):
        end_step = iStep+1
        break  
    
# Wipe Analysis
ops.wipe()
del PondingLoadManager

# Trim results
data_volume_1 = data_volume_1[:end_step]
data_height_1 = data_height_1[:end_step]


# ---------------------------------------------------------------------- #
#   Run second OpenSees analysis (with one CBDI element along lenth)     #
# ---------------------------------------------------------------------- #
# set modelbuilder
ops.model('basic', '-ndm', 2, '-ndf', 3)

# create nodes
ops.node(0,0,zi)
ops.node(1,L,zj)

# set boundary condition
ops.fix(0,1,1,0)
ops.fix(1,0,1,0)

# define coordinate transformation
ops.geomTransf('Linear',1)

# define cross section
ops.section('Elastic', 1, E, A, Iz) 
ops.beamIntegration('Lobatto', 1, 1, 3)

# Time series for loads
ops.timeSeries("Constant", 1)

# define elements
ops.element("dispBeamColumn",0,0,1,1,1)            

# Dead load
ops.pattern('Plain',-1,1)
iNodes = ops.eleNodes(0)
xI,yI = ops.nodeCoord(iNodes[0])
xJ,yJ = ops.nodeCoord(iNodes[1])
ele_angle = atan2(yJ-yI, xJ-xI)
ops.eleLoad('-ele',0,'-type','beamUniform',-wD*cos(ele_angle),-wD*sin(ele_angle))
    
# ------------------------------
# Start of analysis generation
# ------------------------------

# create SOE
ops.system("UmfPack")

# create DOF number
ops.numberer("RCM")

# create constraint handler
ops.constraints("Plain")

# create integrator
ops.integrator("LoadControl", 0.0)

# create algorithm
ops.algorithm("Newton")

# create analysis object
ops.analysis("Static")

# Run dead load analysis
ops.analyze(1)
zo = float('inf')
for i in range(1):
    iz = ops.nodeCoord(i, 2) + ops.nodeDisp(i, 2)
    if iz < zo:
        zo = iz

# Initilize data        
data_volume_2 = np.zeros((num_steps+1,1))
data_height_2 = np.full((num_steps+1,1),zo)
data_reactIx_2 = np.zeros((num_steps+1,1))
data_reactIy_2 = np.zeros((num_steps+1,1))
data_reactJy_2 = np.zeros((num_steps+1,1))
end_step = num_steps+1

ops.reactions()
data_reactIx_2[0] = ops.nodeReaction(0,1)
data_reactIy_2[0] = ops.nodeReaction(0,2)
data_reactJy_2[0] = ops.nodeReaction(1,2)

# ------------------------------
# Finally perform the analysis
# ------------------------------
num_cells = 20;

# define ponding load cells    
PondingLoadManager = PondingLoadManager2d()
for i in range(0,num_cells):
    if i == 0:
        endI = ('node',0)
    else:
        endI = ('element',0,i/num_cells)
        
    if i == (num_cells-1):
        endJ = ('node',1)
    else:
        endJ = ('element',0,(i+1)/num_cells)

    PondingLoadManager.add_cell(i,endI,endJ,gamma,TW)

# Perform analysis, ramping up volume      
zw = zo+0.1

for iStep in range(0,num_steps):
   
    target_volume = (iStep+1)/num_steps*max_volume
    
    # Update ponding load cells
    PondingLoadManager.update()
    
    # Estimate water height
    for i in range(nsteps_vol):
        (V,dVdz) = PondingLoadManager.get_volume(zw)        
        zw = zw - (V-target_volume)/dVdz
        if abs(target_volume-V) <= vol_tol:
            break 
    
    # Compute load vector
    PondingLoadManager.compute_current_load_vector(zw)
        
    # Apply difference to model
    ops.pattern("Plain", iStep, 1)
    PondingLoadManager.apply_load_increment()
    PondingLoadManager.commit_current_load_vector()

    # Run analysis
    ops.analyze(1)
    ops.reactions()
    
    # Store Data
    data_volume_2[iStep+1] = target_volume
    data_height_2[iStep+1] = zw
    data_reactIx_2[iStep+1] = ops.nodeReaction(0,1)
    data_reactIy_2[iStep+1] = ops.nodeReaction(0,2)
    data_reactJy_2[iStep+1] = ops.nodeReaction(1,2)

    # Stop analysis if water level too low
    if (zw-zo) <= (1-percent_drop)*(np.amax(data_height_2)-zo):
        end_step = iStep+1
        break  
    
# Wipe Analysis
ops.wipe()

# Trim results
data_volume_2 = data_volume_2[:end_step]
data_height_2 = data_height_2[:end_step]


# ---------------------------------------------------------------------- #
#   Run PyPonding analyses (with and without ponding effect)             #
# ---------------------------------------------------------------------- #
num_analyses = 20

elastic_beam = steel_beam();
elastic_beam.L   = L
elastic_beam.tw  = TW
elastic_beam.zi  = zi
elastic_beam.zj  = zj
elastic_beam.c   = 0.0

elastic_beam.E   = E
elastic_beam.A   = A
elastic_beam.I   = Iz

elastic_beam.alpha   = 1.0
elastic_beam.LF_D    = 1.0
elastic_beam.wd      = wD/TW
elastic_beam.LF_P    = 1.0
elastic_beam.gamma   = gamma
elastic_beam.LF_S1   = 0.0
elastic_beam.LF_S2   = 0.0
elastic_beam.gammas  = 0.0
elastic_beam.hs      = 0.0
elastic_beam.BuildModel();

data_height_3 = np.linspace(0.0,np.amax(data_height_1),num=num_analyses)
data_volume_3 = np.empty(num_analyses)
data_reactIx_3 = np.empty(num_analyses)
data_reactIy_3 = np.empty(num_analyses)
data_reactJy_3 = np.empty(num_analyses)

PA = FE.PondingAnalysis(elastic_beam.model,'Constant_Level')
PA.max_iterations_z = 60
for i in range(len(data_height_3)):
    res = PA.run({'DEAD':1.0},data_height_3[i])
    if res == 0:
        (V,dVdz) = elastic_beam.model.GetPondingVolume(PA.d,data_height_3[i])
        data_volume_3[i] = V
        (Ri,Rj) = elastic_beam.Reactions(PA)
        data_reactIx_3[i] = 0.0
        data_reactIy_3[i] = Ri
        data_reactJy_3[i] = Rj
        
    else:
        print('Not converged')

data_height_4 = np.linspace(0.0,np.amax(data_height_1),num=num_analyses)
data_volume_4 = np.empty(num_analyses)
data_reactIx_4 = np.empty(num_analyses)
data_reactIy_4 = np.empty(num_analyses)
data_reactJy_4 = np.empty(num_analyses)

PA = FE.PondingAnalysis(elastic_beam.model,'No_Ponding_Effect')
PA.max_iterations_z = 60
for i in range(len(data_height_4)):
    res = PA.run({'DEAD':1.0},data_height_4[i])
    if res == 0:
        (V,dVdz) = elastic_beam.model.GetPondingVolume(PA.d,data_height_4[i])
        data_volume_4[i] = V
        (Ri,Rj) = elastic_beam.Reactions(PA)
        data_reactIx_4[i] = 0.0
        data_reactIy_4[i] = Ri
        data_reactJy_4[i] = Rj
    else:
        print('Not converged')


# ---------------------------------------------------------------------- #
#   Plot results                                                         #
# ---------------------------------------------------------------------- #
plt.figure(1)  
line1, = plt.plot(data_volume_1, data_height_1, '-')
line2, = plt.plot(data_volume_2, data_height_2, '-')
line3, = plt.plot(data_volume_3, data_height_3, 'o')
line4, = plt.plot(data_volume_4, data_height_4, 'o')
plt.legend((line1, line2, line3, line4), ('OpenSees (multiple elements)', 'OpenSees (one element)', 'PyPonding', 'PyPonding (no ponding effect)'),loc='lower right',frameon=False)
plt.xlabel('Water Volume')
plt.ylabel('Water Height')

plt.figure(2)  
line1, = plt.plot(data_reactIx_1, data_height_1, '-')
line2, = plt.plot(data_reactIx_2, data_height_2, '-')
line3, = plt.plot(data_reactIx_3, data_height_3, 'o')
line4, = plt.plot(data_reactIx_4, data_height_4, 'o')
plt.legend((line1, line2, line3, line4), ('OpenSees (multiple elements)', 'OpenSees (one element)', 'PyPonding', 'PyPonding (no ponding effect)'),loc='lower right',frameon=False)
plt.xlabel('Horizontal Reaction')
plt.ylabel('Water Height')

plt.figure(3)  
line1, = plt.plot(data_reactIy_1, data_height_1, '-')
line2, = plt.plot(data_reactIy_2, data_height_2, '-')
line3, = plt.plot(data_reactIy_3, data_height_3, 'o')
line4, = plt.plot(data_reactIy_4, data_height_4, 'o')
plt.legend((line1, line2, line3, line4), ('OpenSees (multiple elements)', 'OpenSees (one element)', 'PyPonding', 'PyPonding (no ponding effect)'),loc='lower right',frameon=False)
plt.xlabel('Vertical Reaction at I-end')
plt.ylabel('Water Height')

plt.show()