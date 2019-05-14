#import openseespy.opensees as ops
import sys
sys.path.append('/home/mhscott/OpenSees/SRC/interpreter')
import opensees as ops

from math import pi,cos,cosh,ceil
from scipy.stats import norm

sys.path.append('/home/mhscott/PyPonding')
from PyPonding.PondingLoadCell import PondingLoadCell2d
import numpy as np
import matplotlib.pyplot as plt

class PondingLoadCell2d_OPS(PondingLoadCell2d):
    def __init__(self,id,nodeI,nodeJ,gamma,tw):
        self.id = id
        self.nodeI = nodeI
        self.nodeJ = nodeJ
        self.gamma = gamma
        self.tw = tw
        
        # Retreive Node Coordinates
        self.xI = ops.nodeCoord(self.nodeI,1)
        self.yI = ops.nodeCoord(self.nodeI,2)
        self.xJ = ops.nodeCoord(self.nodeJ,1)
        self.yJ = ops.nodeCoord(self.nodeJ,2)
        
    def update(self):
        # Code currently only updates y postion of nodes - @todo maybe update x position also
        # self.dxI = ops.nodeDisp(self.nodeI,1)
        self.dyI = ops.nodeDisp(self.nodeI,2)
        # self.dxJ = ops.nodeDisp(self.nodeJ,1)
        self.dyJ = ops.nodeDisp(self.nodeJ,2)
        

class wf:
    def __init__(self,d,tw,bf,tf,Fy,E):
        self.d  = d
        self.tw = tw
        self.bf = bf
        self.tf = tf
        self.Fy = Fy
        self.E  = E
        self.material_type = 'ElasticPP'
        self.num_fiber = 20
        
    def dw(self):
        dw = self.d-2*self.tf
        return dw
    
    def A(self):
        A = 2*self.bf*self.tf + (self.d-2*self.tf)*self.tw
        return A
    
    def Iz(self):
        Iz = (1/12)*self.bf*self.d**3 - (1/12)*(self.bf-self.tw)*self.dw()**3
        return Iz
        
    def define_fiber_section(self,secTag,matTag):
        if self.material_type == 'Elastic':
            ops.uniaxialMaterial('Elastic', matTag, self.E)
        elif self.material_type == 'ElasticPP':
            # ops.uniaxialMaterial('ElasticPP', matTag, self.E, self.Fy/self.E)
            ops.uniaxialMaterial('Steel01', matTag, self.Fy, self.E, 0.001)
            #ops.uniaxialMaterial('Hardening', matTag, self.E, self.Fy, 0.0, 0.001*E)
        else:
            raise Exception('Input Error - unknown material type (%s)' % self.material_type)
        ops.section('Fiber', secTag)
        numSubdivY = ceil(self.tf*(self.num_fiber/self.d))
        ops.patch('rect', matTag, numSubdivY, 1, self.dw()/2, -self.bf/2, self.d/2, self.bf/2)
        numSubdivY = ceil(self.dw()*(self.num_fiber/self.d))
        ops.patch('rect', matTag, numSubdivY, 1, -self.dw()/2, -self.tw/2, self.dw()/2, self.tw/2)
        numSubdivY = ceil(self.tf*(self.num_fiber/self.d))
        ops.patch('rect', matTag, numSubdivY, 1, -self.d/2, -self.bf/2, -self.dw()/2, self.bf/2)
        return
  

nsteps = 100
data_volume = np.zeros((nsteps+1,2))
data_height = np.zeros((nsteps+1,2))
end_step = [nsteps, nsteps]

# Number of Monte Carlo trials
Ntrials = 10

material_types = ['Elastic','ElasticPP']
  
# input parameters
for iAnalysis in range(1):

    E = 29000.0
    Fy = 50.0
    wf_section = wf(5.99,0.230,5.99,0.260,Fy,E) # W6x15
    wf_section.material_type = material_types[iAnalysis]
    # print(wf_section.A()) # 4.43 from the Steel Manual
    # print(wf_section.Iz()) # 29.1 from the Steel Manual

    L       = 480.0
    A       = wf_section.A()
    Iz      = wf_section.Iz()
    gamma   = 62.4/1000/12**3
    tw      = 60.0
    zi      = 0.0
    zj      = 10.0

    max_volume = 10*L*tw

    nsteps_vol = 30
    nele = 20
    vol_tol = max_volume/nsteps/100
    mid_node = int(nele/2)

    # set modelbuilder
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # create nodes
    for i in range(nele+1):
        ops.node(i,L*i/(nele),zi+(zj-zi)*i/(nele))

    # set boundary condition
    ops.fix(   0, 1, 1, 0)
    ops.fix(nele, 0, 1, 0)

    # define coordinate transformation
    ops.geomTransf('Linear',1)

    # define cross section
    wf_section.define_fiber_section(1,1)
    #ops.beamIntegration('Lobatto', 1, 1, 3)
    ops.beamIntegration('Legendre', 1, 1, 2)    

    ops.randomVariable(1,'lognormal','-mean',Fy,'-stdv',0.1*Fy)
    ops.randomVariable(2,'lognormal','-mean', E,'-stdv',0.02*E)    

    ops.parameter(1)
    ops.parameter(2)
    
    # define elements
    for i in range(0,nele):
        # ops.element("elasticBeamColumn",i,i,i+1,A,E,Iz,1)
        ops.element("forceBeamColumn",i,i,i+1,1,1)
        #ops.element("dispBeamColumn",i,i,i+1,1,1)
        ops.addToParameter(1,'element',i,'fy')
        ops.addToParameter(2,'element',i,'E')


        
    # ------------------------------
    # Start of analysis generation
    # ------------------------------

    # create SOE
    ops.system("BandSPD")

    # create DOF number
    ops.numberer("RCM")

    # create constraint handler
    ops.constraints("Plain")

    # create integrator
    ops.integrator("LoadControl", 1.0)

    # create algorithm
    ops.algorithm("Newton")

    ops.probabilityTransformation('Nataf')
    
    # create analysis object
    ops.analysis("Static")

    # Number of failed MC trials
    Nfailed = 0
    
    Nrv = len(ops.getRVTags())
    u = np.zeros(Nrv)
    
    for j in range(Ntrials):

        ops.reset()

        # Transform random variables from standard normal to actual space
        # 1. Create random values in standard normal space
        jj = 0
        for rv in ops.getRVTags():
            u[jj] = norm.ppf(np.random.rand())
            jj = jj+1

        # 2. Transform to real space
        x = ops.transformUtoX(*u)
        print(j,x)
        # 3. Update parameters with random realizations
        jj = 0
        for rv in ops.getRVTags():
            ops.updateParameter(rv,x[jj])
            jj = jj+1


            
        # define ponding load cells    
        PondingLoadCells = dict()
        for i in range(0,nele):
            PondingLoadCells[i] = PondingLoadCell2d_OPS(id,i,i+1,gamma,tw)
            
        # ------------------------------
        # Finally perform the analysis
        # ------------------------------

        # Create dict of each node that can have ponding load applied and initilize load to zero
        EmptyPondingLoad = dict()
        for iCell in PondingLoadCells:
            if not PondingLoadCells[iCell].nodeI in EmptyPondingLoad:
                EmptyPondingLoad[PondingLoadCells[iCell].nodeI] = 0.0    
            if not PondingLoadCells[iCell].nodeJ in EmptyPondingLoad:
                EmptyPondingLoad[PondingLoadCells[iCell].nodeJ] = 0.0


            
        # Perform analysis, ramping up volume      
        zw = 0.1
        CurrentPondingLoad = EmptyPondingLoad.copy()
        for iStep in range(0,nsteps):

            target_volume = (iStep+1)/nsteps*max_volume

            # Update ponding load cells
            for iCell in PondingLoadCells:
                PondingLoadCells[iCell].update()

            # Estimate water height
            for i in range(nsteps_vol):
                V = 0
                dVdz = 0
                for iCell in PondingLoadCells:
                    (iV,idVdz) = PondingLoadCells[iCell].get_volume(zw)
                    V += iV
                    dVdz += idVdz
                    zw = zw - (V-target_volume)/dVdz
                if abs(target_volume-V) <= vol_tol:
                    break 

            # Compute load vector
            UpdatedPondingLoad = EmptyPondingLoad.copy()
            for iCell in PondingLoadCells:    
                f = PondingLoadCells[iCell].get_load_vector(zw)
                UpdatedPondingLoad[PondingLoadCells[iCell].nodeI] += f.item(0)
                UpdatedPondingLoad[PondingLoadCells[iCell].nodeJ] += f.item(1)

            # Apply difference to model
            #
            # Remove time series and pattern so there's not
            # duplicate tag errors in MC loop
            ops.remove('timeSeries',iStep)
            ops.remove('pattern',iStep)
            ops.timeSeries("Linear", iStep)
            ops.pattern("Plain", iStep, iStep)
            for iNode in UpdatedPondingLoad:        
                fy = UpdatedPondingLoad[iNode] - CurrentPondingLoad[iNode]
                ops.load(iNode, 0.0, fy, 0.0)
                CurrentPondingLoad = UpdatedPondingLoad

            # Run analysis
            ops.analyze(1)
            ops.loadConst('-time',0.0)

            # Store Data
            data_volume[iStep+1,iAnalysis] = target_volume
            data_height[iStep+1,iAnalysis] = zw

            # Stop analysis if water level too low
            if zw <= -1:
                end_step[iAnalysis] = iStep+1
                break        
        
    # Wipe Analysis
    ops.wipe()




    
#wi = gamma*zw*tw
#deltai = -5*wi*L**4/(384*E*Iz)
#C = gamma*tw*L**4/(pi**4*E*Iz)
#delta  = deltai/(5*pi**4*C/192/(1/(cos(pi/2*C**0.25))+1/(cosh(pi/2*C**0.25))-2))
#    
#uy = nodeDisp(mid_node,2)
#print('Ponding Amplification: %.3f' % (delta/deltai))    
#print('Closed-form: %.5f' % delta)    
#print('OpenSees:    %.5f' % uy)
#print('Percent Diff %.2f%%' % (100*(uy-delta)/delta))

# Show plot
plt.plot(data_volume[:end_step[0]+1,0], data_height[:end_step[0]+1,0])
plt.plot(data_volume[:end_step[1]+1,1], data_height[:end_step[1]+1,1])
plt.xlabel('Water Volume (in^3)')
plt.ylabel('Water Height (in)')
plt.show()

