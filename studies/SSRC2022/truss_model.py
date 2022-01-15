import numpy as np
import matplotlib.pyplot as plt
#import openseespy.opensees as ops
#import openseespy.postprocessing.Get_Rendering as opsplt
#import openseespy.postprocessing.ops_vis as opsv
import pandas as pd
from math import sqrt
from PyPonding import PondingLoadManager3d
from PyPonding import opensees as ops

#specimen = 'pitched'
specimen = 'flat'

file_suffix = '_linear'

# Geometric Properties
L1 = 4.0
L2 = 576.0
L3 = 4.0
L4 = 67.0
L5 = 5.0

# Define Vertical Profile
if specimen == 'flat':
    def z(x):
        c = 1.0
        L = 576.0
        if c == 0:
            z = 0.
        else:
            r = (c**2 + (L/2)**2)/(2*c)
            z = sqrt(r**2 - (x - L/2)**2) + (c-r)
        return z
elif specimen == 'pitched':
    def z(x):
        c = 1.0
        L = 576.0
        zj = 12.0
        if c == 0:
            z = zj*x/L
        else:
            r = (c**2 + (L/2)**2)/(2*c)
            z = zj*x/L + sqrt(r**2 - (x - L/2)**2) + (c-r)
        return z
else:
    raise Exception('Unknown specimen: %s' % specimen)


# Material Properties 
material_linear = False
geometric_linear = True
E  = 29000.0 
H_iso = E/200.0
H_kin = 0.0
GJ = 1.0e6

#Iz_center = 309.33 # 24K9
Iy_top_chord = 1000.0
A_top_chord = 1.273
Iz_top_chord = 0.490

Iz_edge = 187.74 # 24KSP
A_edge = 1000.0
Iy_edge = 1000.0

# Loads
gamma = 62.4/1000/12**3
target_zw = {'flat':12.5, 'pitched':19.0}
target_V1 = {'flat':60.0/gamma, 'pitched':60.0/gamma}
target_V2 = {'flat':100.0/gamma, 'pitched':100.0/gamma}
num_steps1 = 60*2
num_steps2 = 40*8
num_steps = num_steps1 + num_steps2
na = 4
nb = 4
qD = 5/1000/12**2
xD_center = 0.625 # percent of dead load on center joist. 

# Analysis Options
analysis_type = 'target_volume'
#analysis_type = 'target_level'

N = 48

tol_volume = 0.1        # Tolerance for volume iterations
max_iter_volume = 30    # Maximum number of volume iterations
tol_load_level = 0.0001 # Tolerance for force in 'IterativeLevel' analyses
max_iter_level = 30     # Maximum number of iterations for 'IterativeLevel' analyses
nIP = 5                 # Number of integration points override
#element_type = 'dispBeamColumn'
element_type = 'mixedBeamColumn'
test_flag = 1
include_ponding_effect = True


# Tags    
transf_tag = {'joist': 1, 'beam': 2}
section_tag = {'center_joist_top_chord': 1, 'edge_joist': 2, 'beam': 3}
beamint_tag = {'center_joist_top_chord': 1, 'edge_joist': 2, 'beam': 3}
dead_load_pattern_tag = 1
dead_load_ts_tag = 1
ponding_load_pattern_tag_start = 2
ponding_load_ts_tag = 2  


# Create OpenSees model
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

###########################################################
# Define Joists
###########################################################
if geometric_linear:
    ops.geomTransf('Linear', transf_tag['joist'], 0.0, -1.0, 0.0)
else:
    ops.geomTransf('Corotational', transf_tag['joist'], 0.0, -1.0, 0.0)

if material_linear:
    ops.section('Elastic', section_tag['center_joist_top_chord'], E, A_top_chord, Iz_top_chord, Iy_top_chord, 1, GJ)
else:
    ops.section('Fiber', section_tag['center_joist_top_chord'], '-GJ', GJ)
    matTag_TC = 11;
    ops.uniaxialMaterial('Hardening', matTag_TC, E, 60.15, H_iso, H_kin)
    ops.patch('rect', matTag_TC, 3,  4,  0.395, -2.250, 0.561, -0.250)
    ops.patch('rect', matTag_TC, 3,  4,  0.395,  0.250, 0.561,  2.250)
    ops.patch('rect', matTag_TC, 20, 1, -1.439, -0.416, 0.395, -0.250)
    ops.patch('rect', matTag_TC, 20, 1, -1.439,  0.250, 0.395,  0.416)   
    
ops.beamIntegration('Lobatto', beamint_tag['center_joist_top_chord'], section_tag['center_joist_top_chord'], nIP)

ops.section('Elastic', section_tag['edge_joist'], E, A_edge, Iz_edge, Iy_edge, 1, GJ)
ops.beamIntegration('Lobatto', beamint_tag['edge_joist'], section_tag['edge_joist'], 3)


ops.node(1000, -L1, -L4-L5, z(-L1))
ops.node(2000, -L1,    -L4, z(-L1))
ops.node(3000, -L1,    0.0, z(-L1))
ops.node(4000, -L1,     L4, z(-L1))
ops.node(5000, -L1,  L4+L5, z(-L1))

ops.fix(1000,1,1,0,1,1,1)
ops.fix(2000,0,1,0,1,0,1)
ops.fix(3000,0,1,0,1,0,1)
ops.fix(4000,0,1,0,1,0,1)
ops.fix(5000,1,1,0,1,1,1)
ops.equalDOF(2000,1000,3)
ops.equalDOF(4000,5000,3)

for i in range(N+1):
    n = 110101+i

    x = (i/N)*L2
    ops.node(1001+i, x, -L4-L5, z(x))
    ops.node(2001+i, x,    -L4, z(x))
    ops.node(3001+i, x,    0.0, z(x))
    ops.node(4001+i, x,     L4, z(x))
    ops.node(5001+i, x,  L4+L5, z(x))

    if i == 0:
        ops.fix(1001+i,1,1,1,1,1,1)
        ops.fix(2001+i,1,1,1,1,0,1)
        ops.fix(3001+i,0,1,0,1,0,1)
        ops.fix(4001+i,1,1,1,1,0,1)
        ops.fix(5001+i,1,1,1,1,1,1)   
    elif i == N:  
        ops.fix(1001+i,1,1,1,1,1,1)
        ops.fix(2001+i,0,1,1,1,0,1)
        ops.fix(3001+i,0,1,0,1,0,1)
        ops.fix(4001+i,0,1,1,1,0,1)
        ops.fix(5001+i,1,1,1,1,1,1)
    else:
        ops.fix(1001+i,1,1,0,1,1,1)
        ops.fix(2001+i,0,1,0,1,0,1)
        ops.fix(3001+i,0,1,0,1,0,1)
        ops.fix(4001+i,0,1,0,1,0,1)
        ops.fix(5001+i,1,1,0,1,1,1)
        ops.equalDOF(2001+i,1001+i,3)
        ops.equalDOF(4001+i,5001+i,3)

ops.node(1002+N, L2+L1, -L4-L5, z(L2+L1))
ops.node(2002+N, L2+L1,    -L4, z(L2+L1))
ops.node(3002+N, L2+L1,    0.0, z(L2+L1))
ops.node(4002+N, L2+L1,     L4, z(L2+L1))
ops.node(5002+N, L2+L1,  L4+L5, z(L2+L1))

ops.fix(1002+N,1,1,0,1,1,1)
ops.fix(2002+N,0,1,0,1,0,1)
ops.fix(3002+N,0,1,0,1,0,1)
ops.fix(4002+N,0,1,0,1,0,1)
ops.fix(5002+N,1,1,0,1,1,1)
ops.equalDOF(2002+N,1002+N,3)
ops.equalDOF(4002+N,5002+N,3)

for i in range(N+2):
    ops.element(element_type, 2000+i, 2000+i, 2001+i, transf_tag['joist'], beamint_tag['edge_joist'])
    ops.element(element_type, 3000+i, 3000+i, 3001+i, transf_tag['joist'], beamint_tag['center_joist_top_chord'])
    ops.element(element_type, 4000+i, 4000+i, 4001+i, transf_tag['joist'], beamint_tag['edge_joist'])


###########################################################
# Define Bottom Chord and Web Members
###########################################################

web_A = {'BR15/16': 0.6903, 'C1': 0.2178, 'C1-3/8': 0.4341, 'U1': 0.3498}
web_matTag = {'BR15/16': 81, 'C1': 82, 'C1-3/8': 83, 'U1': 84}
bottom_chord_A = 1.096
bottom_chord_matTag = 85
deff = 22.886


if material_linear:
    ops.uniaxialMaterial('Elastic', web_matTag['BR15/16'], E)
    ops.uniaxialMaterial('Elastic', web_matTag['C1'], E)
    ops.uniaxialMaterial('Elastic', web_matTag['C1-3/8'], E)
    ops.uniaxialMaterial('Elastic', web_matTag['U1'], E)
    ops.uniaxialMaterial('Elastic', bottom_chord_matTag, E)
else:
    ops.uniaxialMaterial('Hardening', web_matTag['BR15/16'], E, 52.2, H_iso, H_kin)
    ops.uniaxialMaterial('Hardening', web_matTag['C1'], E, 62.3, H_iso, H_kin)
    ops.uniaxialMaterial('Hardening', web_matTag['C1-3/8'], E, 64.4, H_iso, H_kin)
    ops.uniaxialMaterial('Hardening', web_matTag['U1'], E, 58.8, H_iso, H_kin)
    ops.uniaxialMaterial('Hardening', bottom_chord_matTag, E, 59.2, H_iso, H_kin)


ops.node(8001,  48.0, 0.0, z(48.0)-deff)
ops.node(8002,  96.0, 0.0, z(96.0)-deff)
ops.node(8003, 144.0, 0.0, z(144.0)-deff)
ops.node(8004, 192.0, 0.0, z(192.0)-deff)
ops.node(8005, 240.0, 0.0, z(240.0)-deff)
ops.node(8006, 288.0, 0.0, z(288.0)-deff)
ops.node(8007, 336.0, 0.0, z(336.0)-deff)
ops.node(8008, 384.0, 0.0, z(384.0)-deff)
ops.node(8009, 432.0, 0.0, z(432.0)-deff)
ops.node(8010, 480.0, 0.0, z(480.0)-deff)
ops.node(8011, 528.0, 0.0, z(528.0)-deff)

ops.fix(8001,0,1,0,1,1,1)
ops.fix(8002,0,1,0,1,1,1)
ops.fix(8003,0,1,0,1,1,1)
ops.fix(8004,0,1,0,1,1,1)
ops.fix(8005,0,1,0,1,1,1)
ops.fix(8006,0,1,0,1,1,1)
ops.fix(8007,0,1,0,1,1,1)
ops.fix(8008,0,1,0,1,1,1)
ops.fix(8009,0,1,0,1,1,1)
ops.fix(8010,0,1,0,1,1,1)
ops.fix(8011,0,1,0,1,1,1)


if geometric_linear:
    truss_ele_type = 'Truss'
else:
    truss_ele_type = 'corotTruss'

ops.element(truss_ele_type, 8001, 3001, 8001, web_A['BR15/16'], web_matTag['BR15/16'])
ops.element(truss_ele_type, 8002, 3004, 8001, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8003, 3007, 8001, web_A['C1-3/8'], web_matTag['C1-3/8'])
ops.element(truss_ele_type, 8004, 3007, 8002, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8005, 3009, 8002, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8006, 3011, 8002, web_A['U1'], web_matTag['U1'])
ops.element(truss_ele_type, 8007, 3011, 8003, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8008, 3013, 8003, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8009, 3015, 8003, web_A['U1'], web_matTag['U1'])
ops.element(truss_ele_type, 8010, 3015, 8004, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8011, 3017, 8004, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8012, 3019, 8004, web_A['U1'], web_matTag['U1'])
ops.element(truss_ele_type, 8013, 3019, 8005, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8014, 3021, 8005, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8015, 3023, 8005, web_A['U1'], web_matTag['U1'])
ops.element(truss_ele_type, 8016, 3023, 8006, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8017, 3025, 8006, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8018, 3027, 8006, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8019, 3027, 8007, web_A['U1'], web_matTag['U1'])
ops.element(truss_ele_type, 8020, 3029, 8007, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8021, 3031, 8007, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8022, 3031, 8008, web_A['U1'], web_matTag['U1'])
ops.element(truss_ele_type, 8023, 3033, 8008, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8024, 3035, 8008, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8025, 3035, 8009, web_A['U1'], web_matTag['U1'])
ops.element(truss_ele_type, 8026, 3037, 8009, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8027, 3039, 8009, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8028, 3039, 8010, web_A['U1'], web_matTag['U1'])
ops.element(truss_ele_type, 8029, 3041, 8010, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8030, 3043, 8010, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8031, 3043, 8011, web_A['C1-3/8'], web_matTag['C1-3/8'])
ops.element(truss_ele_type, 8032, 3046, 8011, web_A['C1'], web_matTag['C1'])
ops.element(truss_ele_type, 8033, 3049, 8011, web_A['BR15/16'], web_matTag['BR15/16'])

ops.element(truss_ele_type, 8101, 8001, 8002, bottom_chord_A, bottom_chord_matTag)
ops.element(truss_ele_type, 8102, 8002, 8003, bottom_chord_A, bottom_chord_matTag)
ops.element(truss_ele_type, 8103, 8003, 8004, bottom_chord_A, bottom_chord_matTag)
ops.element(truss_ele_type, 8104, 8004, 8005, bottom_chord_A, bottom_chord_matTag)
ops.element(truss_ele_type, 8105, 8005, 8006, bottom_chord_A, bottom_chord_matTag)
ops.element(truss_ele_type, 8106, 8006, 8007, bottom_chord_A, bottom_chord_matTag)
ops.element(truss_ele_type, 8107, 8007, 8008, bottom_chord_A, bottom_chord_matTag)
ops.element(truss_ele_type, 8108, 8008, 8009, bottom_chord_A, bottom_chord_matTag)
ops.element(truss_ele_type, 8109, 8009, 8010, bottom_chord_A, bottom_chord_matTag)
ops.element(truss_ele_type, 8110, 8010, 8011, bottom_chord_A, bottom_chord_matTag)

###########################################################
# Define Beams
###########################################################
if geometric_linear:
    ops.geomTransf('Linear', transf_tag['beam'], 1.0, 0.0, 0.0)
else:
    ops.geomTransf('Corotational', transf_tag['beam'], 1.0, 0.0, 0.0)

ops.section('Elastic', section_tag['beam'], E, 14.4, 272, 93.4, 11200.0, 1.39) # W10x49
ops.beamIntegration('Lobatto', beamint_tag['beam'], section_tag['beam'], 3)

ops.node(6000, 0.0, -L4, z(0.0))
ops.node(6001, 0.0, 0.0, z(0.0))
ops.node(6002, 0.0,  L4, z(0.0))
ops.fix(6000,1,0,1,0,1,0)
ops.fix(6002,1,1,1,0,1,0)
ops.equalDOF(6001,3001,1,3)
ops.element(element_type, 6000, 6000, 6001, transf_tag['beam'], beamint_tag['beam'])
ops.element(element_type, 6001, 6001, 6002, transf_tag['beam'], beamint_tag['beam'])

ops.node(7000, L2, -L4, z(L2))
ops.node(7001, L2, 0.0, z(L2))
ops.node(7002, L2,  L4, z(L2))
ops.fix(7000,1,0,1,0,1,0)
ops.fix(7002,1,1,1,0,1,0)
ops.equalDOF(7001,3001+N,3)   # @todo switch back to only dof 3
ops.element(element_type, 7000, 7000, 7001, transf_tag['beam'], beamint_tag['beam'])
ops.element(element_type, 7001, 7001, 7002, transf_tag['beam'], beamint_tag['beam'])


###########################################################
# Define Ponding Load Cells
###########################################################

PondingLoadManager = PondingLoadManager3d()

for i in range(4):
    for j in range(N+2):
        id = '%i_%03i' % (i,j)
        vertexI = ('node',(i+1)*1000+j)
        vertexJ = ('node',(i+1)*1000+j+1)
        vertexK = ('node',(i+2)*1000+j+1)
        vertexL = ('node',(i+2)*1000+j)
        PondingLoadManager.add_cell(id,vertexI,vertexJ,vertexK,vertexL,gamma,na,nb)


###########################################################
# Define Dead Load
###########################################################

ops.timeSeries("Constant", dead_load_ts_tag)
ops.pattern("Plain", dead_load_pattern_tag, dead_load_ts_tag)

# Define uniform dead load on joists
for i in range(N+2):

    n = 2000+i
    wD = 0.5*(1-xD_center)*(2*L4+2*L5)*qD
    nodes = ops.eleNodes(n)
    coordi = ops.nodeCoord(nodes[0])
    coordj = ops.nodeCoord(nodes[1])
    Lx = coordj[0]-coordi[0]
    Lz = coordj[2]-coordi[2]
    L = sqrt(Lx**2 + Lz**2)
    Wy = -wD*(Lx/L)*(Lx/L)
    Wx = -wD*(Lx/L)*(Lz/L)
    ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

    n = 3000+i
    wD = xD_center*(2*L4+2*L5)*qD
    nodes = ops.eleNodes(n)
    coordi = ops.nodeCoord(nodes[0])
    coordj = ops.nodeCoord(nodes[1])
    Lx = coordj[0]-coordi[0]
    Lz = coordj[2]-coordi[2]
    L = sqrt(Lx**2 + Lz**2)
    Wy = -wD*(Lx/L)*(Lx/L)
    Wx = -wD*(Lx/L)*(Lz/L)
    ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

    n = 4000+i
    wD = 0.5*(1-xD_center)*(2*L4+2*L5)*qD
    nodes = ops.eleNodes(n)
    coordi = ops.nodeCoord(nodes[0])
    coordj = ops.nodeCoord(nodes[1])
    Lx = coordj[0]-coordi[0]
    Lz = coordj[2]-coordi[2]
    L = sqrt(Lx**2 + Lz**2)
    Wy = -wD*(Lx/L)*(Lx/L)
    Wx = -wD*(Lx/L)*(Lz/L)
    ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)


###########################################################
# Run Analysis
###########################################################

# Set up analysis
ops.numberer("RCM")
ops.constraints("Plain")
ops.system("ProfileSPD")
ops.test("NormUnbalance", 1.0e-6, 25, test_flag)
ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0)
ops.analysis("Static")
ops.analyze(1)
ops.reactions()



#opsplt.plot_model()  # command from Get_Rendering module
#opsv.plot_model()  # command from ops_vis module



# Initialize results
water_volume = np.empty((num_steps+1,1))
water_level  = np.empty((num_steps+1,1))
water_volume[:] = np.NaN
water_level[:]  = np.NaN

#print(ops.nodeReaction(6000))
#print(ops.nodeReaction(6001))
#print(ops.nodeReaction(6002))

#print(ops.nodeReaction(7000))
#print(ops.nodeReaction(7001))
#print(ops.nodeReaction(7002))

#print(ops.nodeReaction(2001))
#print(ops.nodeReaction(3001))
#print(ops.nodeReaction(4001))

# Find lowest point
zo = PondingLoadManager.find_lowest_point()
print(f'{zo = } in.')

# Store Reuslts
water_volume[0] = 0.
water_level[0] = zo

# Run Ponding Analysis
if analysis_type == 'target_level':

    ops.timeSeries("Constant", ponding_load_ts_tag)
    for iStep in range(1,num_steps+1):

        # Update ponding load cells
        if include_ponding_effect:
            PondingLoadManager.update()

        # Compute load vector
        izw = zo + (iStep/num_steps)*(target_zw[specimen]-zo)
        (iV,idVdz) = PondingLoadManager.get_volume(izw)
        PondingLoadManager.compute_current_load_vector(izw)

        # Apply difference to model
        ops.pattern("Plain", ponding_load_pattern_tag_start+iStep, ponding_load_ts_tag)
        PondingLoadManager.apply_load_increment()
        PondingLoadManager.commit_current_load_vector()

        # Run analysis
        ret = ops.analyze(1)
        if ret < 0:
            break
            
        ops.reactions()

        # Store Reuslts
        water_volume[iStep] = iV
        water_level[iStep]  = izw

elif analysis_type == 'target_volume':

    izw = zo+0.1
    ops.timeSeries("Constant", ponding_load_ts_tag)
    for iStep in range(1,num_steps+1):

        # Update ponding load cells
        if include_ponding_effect:
            PondingLoadManager.update()
            
        # Compute load vector
        nsteps_vol = 30
        vol_tol = 0.1
        volume_found = False
        if iStep <= num_steps1:
            target_iV = (iStep/num_steps1)*target_V1[specimen]
        else:
            target_iV = target_V1[specimen] + ((iStep-num_steps1)/num_steps2)*(target_V2[specimen]-target_V1[specimen])
        print('Step: %i, Target W: %.3f' % (iStep,target_iV*gamma))
        for i in range(nsteps_vol):
            (iV,idVdz) = PondingLoadManager.get_volume(izw)
            izw = izw - (iV-target_iV)/idVdz
            if abs(target_iV-iV) <= vol_tol:
                volume_found = True
                break 
        if volume_found == False:
            break
        PondingLoadManager.compute_current_load_vector(izw)

        # Store Reuslts
        water_volume[iStep] = iV
        water_level[iStep]  = izw

        # Apply difference to model
        ops.pattern("Plain", ponding_load_pattern_tag_start+iStep, ponding_load_ts_tag)
        PondingLoadManager.apply_load_increment()
        PondingLoadManager.commit_current_load_vector()

        # Run analysis
        ret = ops.analyze(1)
        if ret < 0:
            break
            #ops.integrator("LoadControl", 0.1)
            #ops.test("NormUnbalance", 1.0e-6, 1, 5)
            #ret = ops.analyze(10)
        
            #if ret == 0:
            #    ops.integrator("LoadControl", 0.1)
            #    ops.test("NormUnbalance", 1.0e-6, 25, test_flag)
            #else:
            #    break
        ops.reactions()

else:
    raise Exception('Unknown analysis_type: %s' % analysis_type)

###########################################################
# Plot Deformed Shape
###########################################################
scale_factor = 1.0

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

NodeCoord = dict()
NodeDisp = dict()
NodeTags = ops.getNodeTags()
for i in NodeTags:
    NodeCoord[i] = ops.nodeCoord(i)
    NodeDisp[i] = ops.nodeDisp(i)
    #ax.scatter(coord[0], coord[1], coord[2], marker='o')

EleTags = ops.getEleTags()
for i in EleTags:
    nodes = ops.eleNodes(i)
    coordi = NodeCoord[nodes[0]]
    coordj = NodeCoord[nodes[1]]
    dispi = NodeDisp[nodes[0]]
    dispj = NodeDisp[nodes[1]]
    xplt = [coordi[0]+scale_factor*dispi[0],coordj[0]+scale_factor*dispj[0]]
    yplt = [coordi[1]+scale_factor*dispi[1],coordj[1]+scale_factor*dispj[1]]
    zplt = [coordi[2]+scale_factor*dispi[2],coordj[2]+scale_factor*dispj[2]]
    ax.plot(xplt,yplt,zplt)

plt.show()


###########################################################
# Save Results to Text File
###########################################################
out_file = open('OSU_data\PyPonding_%s%s.csv' % (specimen,file_suffix),'w') 
out_file.write('total_load_kips,water_level_inches\n')
for i in range(num_steps+1):
    if np.isnan(water_volume[i]):
        break
    out_file.write('%.6f,%.6f\n' % (water_volume[i]*gamma,water_level[i]))
out_file.close()


###########################################################
# Make Plot
###########################################################

# Data from OpenSeesPy and PyPonding 
plt.plot(water_volume*gamma, water_level, 'o-')

# Data from OSU Experiment
exp_data = pd.read_csv('OSU_data\Figure_68_%s.csv' % specimen)
plt.plot(exp_data['total_load_kips'].to_numpy(), exp_data['water_level_inches'].to_numpy())

plt.xlabel('Total Load (kips)')
plt.ylabel('Water Level (in)')
plt.show()





# # Plot Ponding Load Cells
# if plot_load_cells:
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')
#     for i in PondingLoadManager.cells:
#         coordI = PondingLoadManager.cells[i].vertexI.coord()
#         coordJ = PondingLoadManager.cells[i].vertexJ.coord()
#         coordK = PondingLoadManager.cells[i].vertexK.coord()
#         coordL = PondingLoadManager.cells[i].vertexL.coord()
#         #if i=='BC_001x_001y':
#         #    print(coordI)
#         #    print(coordJ)
#         #    print(coordK)
#         #    print(coordL)
#         X = np.array([[coordJ[0], coordK[0]], [coordI[0], coordL[0]]])
#         Y = np.array([[coordJ[1], coordK[1]], [coordI[1], coordL[1]]])
#         Z = np.array([[coordJ[2], coordK[2]], [coordI[2], coordL[2]]])
#         surf = ax.plot_surface(X, Y, Z, linewidth=1)
#     plt.show()
