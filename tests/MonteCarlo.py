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
    def __init__(self,d,tw,bf,tf,Fy,E,Hk):
        self.d  = d
        self.tw = tw
        self.bf = bf
        self.tf = tf
        self.Fy = Fy
        self.E  = E
        self.Hk = Hk
        self.material_type = 'ElasticPP'
        self.num_fiber = 20
        
    def dw(self):
        dw = self.d-2*self.tf
        return dw
    
    def A(self):
        A = 2*self.bf*self.tf + (self.d-2*self.tf)*self.tw
        return A
    
    def Iz(self):
        Iz = (1.0/12)*self.bf*self.d**3 - (1.0/12)*(self.bf-self.tw)*self.dw()**3
        return Iz
        
    def define_fiber_section(self,secTag,matTag):
        if self.material_type == 'Elastic':
            ops.uniaxialMaterial('Elastic', matTag, self.E)
        elif self.material_type == 'ElasticPP':
            # ops.uniaxialMaterial('ElasticPP', matTag, self.E, self.Fy/self.E)
            b = self.Hk/(self.E+self.Hk)
            #ops.uniaxialMaterial('Steel01', matTag, self.Fy, self.E, b)
            ops.uniaxialMaterial('Hardening', matTag, self.E, self.Fy, 0.0, self.Hk)
        else:
            raise Exception('Input Error - unknown material type (%s)' % self.material_type)
        #ops.section('Fiber', secTag)
        #numSubdivY = ceil(self.tf*(self.num_fiber/self.d))
        #ops.patch('rect', matTag, numSubdivY, 1, self.dw()/2, -self.bf/2, self.d/2, self.bf/2)
        #numSubdivY = ceil(self.dw()*(self.num_fiber/self.d))
        #ops.patch('rect', matTag, numSubdivY, 1, -self.dw()/2, -self.tw/2, self.dw()/2, self.tw/2)
        #numSubdivY = ceil(self.tf*(self.num_fiber/self.d))
        #ops.patch('rect', matTag, numSubdivY, 1, -self.d/2, -self.bf/2, -self.dw()/2, self.bf/2)
        
        Nfw = ceil(self.dw()*(self.num_fiber/self.d))
        Nff = ceil(self.tf*(self.num_fiber/self.d))
        ops.section('WFSection2d',secTag,matTag,self.d,self.tw,self.bf,self.tf,Nfw,Nff)
        return
  

nsteps = 100
nsteps = 250

material_type = 'Elastic'
material_type = 'ElasticPP'

inch = 1.0
kip = 1.0
minute = 1.0

lb = kip/1000.0
ft = 12.0*inch
hr = 60*minute

gal = 0.133681*ft**3

ksi = kip/inch**2
psf = lb/ft**2
pcf = psf/ft

Fy = 50.0*ksi
E = 29000.0*ksi
Hk = 29.0*ksi

# W6x15
d = 5.99*inch
tweb = 0.230*inch
bf = 5.99*inch
tf = 0.260*inch

# W14x22
d = 13.74*inch
bf = 5.0*inch
tf = 0.335*inch
tweb = 0.23*inch

#wf_section = wf(5.99,0.230,5.99,0.260,Fy,E) # W6x15
wf_section = wf(d,tweb,bf,tf,Fy,E,Hk) # W6x15
wf_section.material_type = material_type
# print(wf_section.A()) # 4.43 from the Steel Manual
# print(wf_section.Iz()) # 29.1 from the Steel Manual

L       = 45*ft # span length
A       = wf_section.A()
Iz      = wf_section.Iz()
gamma   = 62.4*pcf
tw      = 5*ft # tributary width
zi      = 0.0*inch
# Max water height below hFail is failure
zj      = 0*inch;     hFail = 4.9*inch
#zj      = 11.25*inch; hFail = 10.5*inch

# Hourly rainfall rate (in/hr)
rate = 3.75*inch/hr




qD = 20.0*psf # dead load
#qD = 0.01*psf
wD = qD*tw # Dead load per length

max_volume = (15*inch)*L*tw

nsteps_vol = 30
nele = 20
vol_tol = max_volume/nsteps/100.0
mid_node = int(nele/2)

# set modelbuilder
ops.wipe()
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
Np = 4
ops.beamIntegration('Lobatto', 1, 1, Np)
#ops.beamIntegration('Legendre', 1, 1, 2)    

EIplot = np.zeros((nele*Np,nsteps+1))

# Time series for loads
ops.timeSeries("Constant", 1)

# Dead load
ops.pattern('Plain',-1,1)

    
ops.randomVariable(1,'lognormal','-mean',Fy,'-stdv',0.1*Fy)
ops.randomVariable(2,'lognormal','-mean', E,'-stdv',0.02*E)
ops.randomVariable(3,'normal','-mean',d,'-stdv',0.02*d)
ops.randomVariable(4,'normal','-mean',tweb,'-stdv',0.02*tweb)
ops.randomVariable(5,'normal','-mean',bf,'-stdv',0.02*bf)
ops.randomVariable(6,'normal','-mean',tf,'-stdv',0.02*tf)
ops.randomVariable(7,'normal','-mean',-wD,'-stdv',-0.1*wD)
ops.randomVariable(8,'lognormal','-mean',Hk,'-stdv',0.05*Hk)

ops.probabilityTransformation('Nataf')

ops.parameter(1)
ops.parameter(2)
ops.parameter(3)
ops.parameter(4)
ops.parameter(5)
ops.parameter(6)
ops.parameter(7)
ops.parameter(8)

gammaRVTag = 9
ops.randomVariable(gammaRVTag,'lognormal','-mean',gamma,'-stdv',0.1*gamma)
ops.parameter(gammaRVTag)
ops.updateParameter(gammaRVTag,gamma)

rateRVTag = 10
ops.randomVariable(rateRVTag,'lognormal','-mean',rate,'-stdv',0.2*rate)
ops.parameter(rateRVTag)
ops.updateParameter(rateRVTag,rate)

# define elements
for i in range(nele):
    # ops.element("elasticBeamColumn",i,i,i+1,A,E,Iz,1)
    ops.element("forceBeamColumn",i,i,i+1,1,1)
    #ops.element("dispBeamColumn",i,i,i+1,1,1)
    ops.addToParameter(1,'element',i,'fy')
    ops.addToParameter(2,'element',i,'E')
    ops.addToParameter(3,'element',i,'d')
    ops.addToParameter(4,'element',i,'tw')
    ops.addToParameter(5,'element',i,'bf')
    ops.addToParameter(6,'element',i,'tf')
    ops.addToParameter(8,'element',i,'Hkin')    

    ops.eleLoad('-ele',i,'-type','beamUniform',-wD)
    ops.addToParameter(7,'loadPattern',-1,'elementLoad',i,'wy')

legendLabel = {
    9: 'gamma',
    10: 'rate',
    1: 'Fy',
    2: 'E',
    3: 'd',
    4: 'tw',
    5: 'bf',
    6: 'tf',
    8: 'Hkin',
    7: 'wD'
}

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
    
ops.probabilityTransformation('Nataf')
    
# create analysis object
ops.analysis("Static")


Nparam = len(ops.getParamTags())
Nrv = len(ops.getRVTags())
u = np.zeros(Nrv)

# Number of Monte Carlo trials
Ntrials = 1000
Ntrials = 100

# Number of failed MC trials
Nfailed = 0

duPlot = np.zeros((nsteps+1,Nparam))
meanPlot = np.zeros((nsteps+1,2))

plt.figure(1)
plt.subplot(2,1,1)

for j in range(Ntrials+1):

    ops.reset()

    # Transform random variables from standard normal to actual space
    # 1. Create random values in standard normal space
    jj = 0
    for rv in ops.getRVTags():
        u[jj] = norm.ppf(np.random.rand())
        if j == Ntrials:
            u[jj] = 0
        jj = jj+1

    # 2. Transform to real space
    x = ops.transformUtoX(*u)
    #print(j,x)
    print(j)
    # 3. Update parameters with random realizations
    jj = 0
    for rv in ops.getRVTags():
        ops.updateParameter(rv,x[jj])
        jj = jj+1

    # Dead load analysis
    ops.analyze(1)

    # define ponding load cells
    # Inside the MC loop so we can make gamma random
    #
    gamma = x[gammaRVTag-1]
    PondingLoadCells = dict()
    for i in range(nele):
        PondingLoadCells[i] = PondingLoadCell2d_OPS(id,i,i+1,gamma,tw)

    # ------------------------------
    # Finally perform the analysis
    # ------------------------------

    data_volume = np.zeros(nsteps+1)
    data_height = np.zeros(nsteps+1)
    end_step = nsteps

    # Create dict of each node that can have ponding load applied and initilize load to zero
    EmptyPondingLoad = dict()
    for iCell in PondingLoadCells:
        if not PondingLoadCells[iCell].nodeI in EmptyPondingLoad:
            EmptyPondingLoad[PondingLoadCells[iCell].nodeI] = 0.0    
        if not PondingLoadCells[iCell].nodeJ in EmptyPondingLoad:
            EmptyPondingLoad[PondingLoadCells[iCell].nodeJ] = 0.0
        

    if j == Ntrials:
        ops.sensitivityAlgorithm('-computeAtEachStep')
            
    # Perform analysis, ramping up volume      
    zw = 0.1
    CurrentPondingLoad = EmptyPondingLoad.copy()
    for iStep in range(nsteps):

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
        ops.pattern("Plain", iStep, 1)
        for iNode in UpdatedPondingLoad:
            fy = UpdatedPondingLoad[iNode] - CurrentPondingLoad[iNode]
            ops.load(iNode, 0.0, fy, 0.0)
        CurrentPondingLoad = UpdatedPondingLoad

        # Run analysis
        ok = ops.analyze(1)

        # Store Data
        data_volume[iStep+1] = target_volume
        data_height[iStep+1] = zw

        if j == Ntrials:
            meanPlot[iStep+1,0] = target_volume
            meanPlot[iStep+1,1] = zw
            for rv in ops.getRVTags():
                cov = ops.getStdv(rv)/ops.getMean(rv)
                # Negative sign because disp is downward
                duPlot[iStep+1,rv-1] = -ops.sensNodeDisp(mid_node,2,rv)*abs(ops.getStdv(rv))
                    
            ix = 0
            for e in range(nele):
                [I,J] = ops.eleNodes(e)
                XI = ops.nodeCoord(I,'X')
                XJ = ops.nodeCoord(J,'X')
                Lele = XJ-XI
                for ip in range(Np):
                    x = ops.sectionLocation(e,ip+1)
                    EIplot[ix,0] = XI + x
                    ks = ops.sectionStiffness(e,ip+1)
                    EIplot[ix,iStep+1] = ks[3]
                    ix = ix+1
                  
        # Stop analysis if water level too low or analysis failed
        if zw <= -4*inch or ok < 0:
            end_step = iStep+1
            break        

    if j < Ntrials and max(data_height) <= hFail:
        Nfailed = Nfailed+1

    # Remove patterns so there's not
    # duplicate tag errors in MC loop
    for iStep in range(nsteps):
        ops.remove('loadPattern',iStep)


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

    #print(zw,end_step)
    # Show plot
    plt.plot(data_volume[:end_step+1], data_height[:end_step+1],'y',linewidth=0.5)
    #plt.plot(data_volume, data_height)

pf = 1.0*Nfailed/Ntrials
print('Probability of failure at less than %.1f inch water height is %.3f\n' % (hFail,pf))

plt.title('W14x22, L=%.1f ft, zj=%.2f in -- Pf=%.3f' % (L/ft,zj,pf))
plt.plot([0,max_volume],[hFail,hFail],'k:')
plt.plot(meanPlot[:end_step+1,0],meanPlot[:end_step+1,1],'k',linewidth=2)
#plt.xlabel('Water Volume (in^3)')
plt.ylabel('Water Height (in)')
plt.ylim([-10*inch,15*inch])
plt.xlim([0,max_volume])
plt.grid()
#plt.show()

plt.subplot(2,1,2)
#plt.title('W14x22, L=%.1f ft, zj=%.2f in -- Pf=%.3f' % (L/ft,zj,pf))
for rv in ops.getRVTags():
    plt.plot(data_volume[:end_step+1],duPlot[:end_step+1,rv-1],label=legendLabel[rv])
plt.xlabel('Water Volume (in^3)')
plt.ylabel('$\Delta u$ (in)')
plt.ylim([-1*inch,0])
plt.xlim([0,max_volume])
plt.legend(loc='lower left')
plt.grid()
plt.show()


plt.figure()
for i in range(nsteps):
    plt.plot(EIplot[:,0],EIplot[:,i+1],'k')
plt.grid()
plt.xlabel('x (in)')
plt.ylabel('EI (kip-in^2)')
plt.title('Flexural Stiffness along Beam at all %d Steps' % nsteps)
plt.show()
