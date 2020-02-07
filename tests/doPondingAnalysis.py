#import openseespy.opensees as ops
import sys
sys.path.append('/home/mhscott/OpenSees/SRC/interpreter')
import opensees as ops

from math import pi,cos,cosh,ceil,log
from scipy.stats import norm

sys.path.append('/home/mhscott/PyPonding')
from PyPonding.PondingLoadCell import PondingLoadCell2d
import numpy as np

from wide_flange import wf,wf_shapes

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

        
def doPondingAnalysis(shape_name,L,slope,qD,label):

    inch = 1.0
    kip = 1.0
    minute = 1.0
    
    lb = kip/1000.0
    ft = 12.0*inch
    hr = 60*minute
    sec = minute/60.0

    ksi = kip/inch**2
    psf = lb/ft**2
    pcf = psf/ft
    
    g = 386.4*inch/sec**2
    gamma   = 62.4*pcf

    gal = 0.133681*ft**3


    material_type = 'Elastic'
    material_type = 'ElasticPP'
    material_type = 'Hardening'
    
    Fy = 50.0*ksi # yield stress
    E = 29000.0*ksi # elastic modulus
    Hk = 29.0*ksi # kinematic hardening modulus
    Fr = 0.2*Fy # residual stress
    Fr = 0.001*Fy
    
    # Hourly rainfall rate (in/hr)
    locations = ['Denver','New York','New Orleans']
    rate = {}
    perc95 = {}
    rate['Denver'] = 1.26*inch/(0.25*hr)
    perc95['Denver'] = 1.72*inch/(0.25*hr)
    rate['New York'] = 1.67*inch/(0.25*hr)
    perc95['New York'] = 2.33*inch/(0.25*hr)
    rate['New Orleans'] = 2.33*inch/(0.25*hr)
    perc95['New Orleans'] = 3.11*inch/(0.25*hr)
    
    # "Drag coefficient" in scupper flow equation
    cd = 0.6
    
    # Scupper width
    ws = 6*inch
    ws = 12*inch

    # Scupper spacing
    Ss = 40*ft
    
    # Tributary width of beams
    tw = 10*ft

    # Number of analysis steps
    nsteps = 100
    nsteps = 250
    nsteps = 3000
    nsteps = 1000
    
    # Number of Monte Carlo trials
    Ntrials = 1000
    #Ntrials = 200
    Ntrials = 99
    Ntrials = 1

    # From function input
    wfs = wf_shapes[shape_name]
    d = wfs['d']*inch
    tweb = wfs['tw']*inch
    bf = wfs['bf']*inch
    tf = wfs['tf']*inch
    dw = d-2*tf
    
    wf_section = wf(d,tweb,bf,tf,Fy,E,Hk)
    wf_section.material_type = material_type

    A  = wf_section.A()
    Iz = wf_section.Iz()

    # From function input
    zi      = 0.0*inch
    zj      = slope*L

    # From function input
    wD = qD*tw # Dead load per length



    Atrib = L*tw
    wf_section.gamma    = gamma
    wf_section.TW    = tw
    wf_section.wD    = wD
    wf_section.L = L
    wf_section.zi = zi
    wf_section.zj = zj
    wf_section.frc = Fr

    # Nominal value for static head
    As = Ss*L
    dhnom = {}
    for city in locations:
        q = rate[city]*As
        dhnom[city] = (1.5*q/(cd*ws*(2*g)**0.5))**(2.0/3)
    
    methods = ['AISC Appendix 2','DAMP','Proposed for ASCE 7','Neglect Ponding']
    ds = {}
    zw_lim = {}
    # Set design method and compute zw_lim (ds+dh)
    for method in methods:
        zw_lim[method] = wf_section.maximum_permitted_zw(method)
        #ds[method] = zw_lim[method] - dh # Now computed in ReadMonteCarlo.py

    max_volume = (300*inch)*Atrib

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

    
    ops.randomVariable(1,'lognormal','-mean',Fy,'-stdv',0.1*Fy); ops.parameter(1)
    ops.randomVariable(2,'lognormal','-mean', E,'-stdv',0.02*E); ops.parameter(2)
    ops.randomVariable(3,'normal','-mean',-wD,'-stdv',0.1*wD); ops.parameter(3)
    ops.randomVariable(4,'lognormal','-mean',Fr,'-stdv',0.15*Fr); ops.parameter(4)

    ops.parameter(8) # For the web
    Nregions = wf_section.num_regions # Number of fibers in each flange
    for i in range(Nregions):
        ops.parameter(9+i)

    
    #ops.randomVariable(3,'normal','-mean',d,'-stdv',0.02*d)
    #ops.randomVariable(4,'normal','-mean',tweb,'-stdv',0.02*tweb)
    #ops.randomVariable(5,'normal','-mean',bf,'-stdv',0.02*bf)
    #ops.randomVariable(6,'normal','-mean',tf,'-stdv',0.02*tf)
    #ops.randomVariable(8,'lognormal','-mean',Hk,'-stdv',0.05*Hk)
    
    ops.probabilityTransformation('Nataf')

    rateRVTag = 5
    for city in locations:
        #ops.randomVariable(rateRVTag,'lognormal','-mean',rate,'-stdv',0.2*rate)
        eulergamma = 0.57721566490153286061
        scale = (rate[city]-perc95[city])/(log(-log(0.95))+eulergamma)
        loc = rate[city] - scale*eulergamma
        #ops.randomVariable(rateRVTag,'type1LargestValue','-parameters',loc,scale)
        ops.randomVariable(rateRVTag,'type1LargestValue','-mean',rate[city],'-stdv',(scale*6**0.5)/pi)
        ops.parameter(rateRVTag)
        ops.updateParameter(rateRVTag,rate[city])
        rateRVTag += 1
    
    
    # define elements
    for i in range(nele):
        # ops.element("elasticBeamColumn",i,i,i+1,A,E,Iz,1)
        #ops.element("forceBeamColumn",i,i,i+1,1,1)
        ops.element("dispBeamColumn",i,i,i+1,1,1)
        ops.addToParameter(1,'element',i,'fy')
        ops.addToParameter(2,'element',i,'E')
        ops.eleLoad('-ele',i,'-type','beamUniform',-wD)
        ops.addToParameter(3,'loadPattern',-1,'elementLoad',i,'wy')
        ops.addToParameter(6,'element',i,'material',1,'F0') # Web
        for j in range(Nregions):
            ops.addToParameter(7+j,'element',i,'material',1+2*(j+1),'F0') # Flange fibers

        #ops.addToParameter(3,'element',i,'d')
        #ops.addToParameter(4,'element',i,'tw')
        #ops.addToParameter(5,'element',i,'bf')
        #ops.addToParameter(6,'element',i,'tf')
        #ops.addToParameter(8,'element',i,'Hkin')

        



    legendLabel = {
        #9: 'gamma',
        9: 'rate',
        10: 'cd',
        1: 'Fy',
        2: 'E',
        3: 'd',
        4: 'tw',
        5: 'bf',
        6: 'tf',
        #8: 'Hkin',
        8: 'Fr',
        7: 'wD'
    }

    # ------------------------------
    # Start of analysis generation
    # ------------------------------

    # create SOE
    ops.system("SparseGeneral")
    
    # create DOF number
    ops.numberer("RCM")
    
    # create constraint handler
    ops.constraints("Plain")
    
    # create integrator
    ops.integrator("LoadControl", 0.0)
    
    ops.test('NormUnbalance',1.0e-6,20,0)
    
    # create algorithm
    ops.algorithm("Newton")
    
    ops.probabilityTransformation('Nataf')
    
    # create analysis object
    ops.analysis("Static")
    

    Nparam = len(ops.getParamTags())
    Nrv = len(ops.getRVTags())
    u = np.zeros(Nrv)

    duPlot = np.zeros((nsteps+1,Nparam))
    meanPlot = np.zeros((nsteps+1,2))

    #plt.figure(1)
    #plt.subplot(2,1,1)

    output = open(f'trials{label}.csv','w')
    output.write(f'{shape_name}, L = {L} in, slope = {slope} in/in, dead load = {qD} kip/in^2\n')
    for method in methods:
        output.write(f'zwlim {method},')
    for city in locations:
        output.write(f'dhnom {city},')
    for city in locations:
        output.write(f'dh {city},')        
    output.write('dmax\n')

    for j in range(Ntrials+1):
        print(j,label)
        
        ops.reset()

        # Transform random variables from standard normal to actual space
        # 1. Create random values in standard normal space
        jj = 0
        for rv in ops.getRVTags():
            u[jj] = norm.ppf(np.random.rand())
            if j == Ntrials: 
                u[jj] = 0 # mean realizations
            jj = jj+1

        # 2. Transform to real space
        x = ops.transformUtoX(*u)

        # 3. Update parameters with random realizations
        jj = 0
        for rv in ops.getRVTags():
            if rv != 4:
                ops.updateParameter(rv,x[jj])
            else:
                Frj = x[jj]
                Frt = -Frj*(bf*tf)/(bf*tf+tw*dw)
                ops.updateParameter(8,Frt)
                for j in range(Nregions):
                    Fri = Frj + (j+0.5)/Nregions*(Frt-Frj)
                    ops.updateParameter(9+j,Fri)
            jj = jj+1

        # Dead load analysis
        ops.analyze(1)

        # define ponding load cells
        # Inside the MC loop so we can make gamma random
        #
        #gamma = x[gammaRVTag-1]
        PondingLoadCells = dict()
        for i in range(nele):
            PondingLoadCells[i] = PondingLoadCell2d_OPS(id,i,i+1,gamma,tw)

        for method in methods:
            output.write(f'{zw_lim[method]},')
        for city in locations:
            output.write(f'{dhnom[city]},')
        rateRVTag = 5
        for city in locations:
            q = x[rateRVTag-1]*As
            dh = (1.5*q/(cd*ws*(2*g)**0.5))**(2.0/3)
            output.write(f'{dh},')
            rateRVTag += 1
    
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
            meanHV = open(f'meanHV{label}.csv','w')
            #ops.sensitivityAlgorithm('-computeAtEachStep')
            
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
            if ok < 0:
                print(target_volume,max_volume,zw)
                        
            # Store Data
            data_volume[iStep+1] = target_volume
            data_height[iStep+1] = zw
            if j == Ntrials:
                meanPlot[iStep+1,0] = target_volume
                meanPlot[iStep+1,1] = zw
                meanHV.write(f'{target_volume},{zw}\n')
                for rv in ops.getRVTags():
                    cov = ops.getStdv(rv)/ops.getMean(rv)
                    # Negative sign because disp is downward
                    #duPlot[iStep+1,rv-1] = -ops.sensNodeDisp(mid_node,2,rv)*abs(ops.getStdv(rv))


            # Stop analysis if water level too low or analysis failed
            if zw <= -4*inch or ok < 0:
                end_step = iStep+1
                break        

        output.write(f'{max(data_height)}')
        output.write('\n')

        if j == Ntrials:
            meanHV.close()
                
        # Remove patterns so there's not
        # duplicate tag errors in MC loop
        for iStep in range(nsteps):
            ops.remove('loadPattern',iStep)



    output.close()

    # Remove random variables bc ops.wipe() doesn't do this yet
    ops.wipeReliability()
