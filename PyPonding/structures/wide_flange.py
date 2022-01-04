import numpy as np
from math import pi, ceil
from PyPonding import PondingLoadCell2d_OPS
from PyPonding.structures import steel_beam

try:
    import opensees as ops
except ImportError:
    try:
        import openseespy.opensees as ops
    except ImportError:
        import warnings
        warnings.warn(
            'OpenSeesPy not found. OpenSees functionality is not available.')


class wf:
    geomTransfType = 'Corotational'
    extra_results = False

    def __init__(self, d, tw, bf, tf, Fy, E, Hk):
        self.d = d             # Depth
        self.tw = tw            # Thickness of the web
        self.bf = bf            # Width of the flange
        self.tf = tf            # Thickness of the flnage
        self.Fy = Fy            # Yield stree
        self.E = E             # Elastic modulus
        self.Hk = Hk            # Kinematic hardening modulus
        # Maximum compressive residual stress (Lehigh pattern)
        self.frc = 0.0
        self.num_regions = 10   # Number of regions for the discretization of residual stress

        self.L = float('nan')  # Beam length (span)
        self.TW = float('nan')  # Tributary width
        self.gamma = 0.0

        self.zi = 0.0
        self.zj = 0.0
        self.c = 0.0           # Camber

        self.wD = 0.0           # Dead load (force per unit length)

        self.material_type = 'Hardening'
        self.num_elements = 20
        self.num_fiber = 20
        self.num_steps = 100

        self.nsteps_vol = 30
        self.max_volume = float('nan')
        self.vol_tol = float('nan')
        self.percent_drop = 0.05

    def dw(self):
        dw = self.d-2*self.tf
        return dw

    def A(self):
        A = 2*self.bf*self.tf + (self.d-2*self.tf)*self.tw
        return A

    def Iz(self):
        Iz = (1.0/12)*self.bf*self.d**3 - (1.0/12) * \
              (self.bf-self.tw)*self.dw()**3
        return Iz

    def Sz(self):
        Sz = self.Iz()/(self.d/2)
        return Sz

    def Zz(self):
        Zz = 2*((self.tf*self.bf)*(self.d/2-self.tf/2) +
                (self.dw()/2*self.tw)*(self.dw()/4));
        return Zz

    def Mp(self):
        return self.Fy*self.Zz()

    def C(self):
        return (self.gamma*self.TW*self.L**4)/(pi**4*self.E*self.Iz())

    def define_fiber_section(self, secTag, matTag):
        Nfw = ceil(self.dw()*(self.num_fiber/self.d))
        Nff = ceil(self.tf*(self.num_fiber/self.d))

        if self.frc == 0 or self.material_type == 'Elastic':
            if self.material_type == 'Elastic':
                ops.uniaxialMaterial('Elastic', matTag, self.E)
            elif self.material_type == 'ElasticPP':
                ops.uniaxialMaterial('ElasticPP', matTag,
                                     self.E, self.Fy/self.E)
            elif self.material_type == 'Steel01':
                b = self.Hk/(self.E+self.Hk)
                ops.uniaxialMaterial('Steel01', matTag, self.Fy, self.E, b)
            elif self.material_type == 'Hardening':
                ops.uniaxialMaterial('Hardening', matTag,
                                     self.E, self.Fy, 0.0, self.Hk)
            else:
                raise Exception(
                    'Input Error - unknown material type (%s)' % self.material_type)

            ops.section('WFSection2d', secTag, matTag, self.d,
                        self.tw, self.bf, self.tf, Nfw, Nff)
        else:
            ops.section('Fiber', secTag)

            frt = -self.frc*(self.bf*self.tf) / \
                             (self.bf*self.tf+self.tw*self.dw())

            # Define web fibers
            if self.material_type == 'ElasticPP':
                ops.uniaxialMaterial(
                    'ElasticPP', matTag, self.E, self.Fy/self.E, -self.Fy/self.E, frt/self.E)
            elif self.material_type == 'Steel01':
                b = self.Hk/(self.E+self.Hk)
                ops.uniaxialMaterial('Steel01', matTag+1, self.Fy, self.E, b)
                ops.uniaxialMaterial('InitStressMaterial',
                                     matTag, matTag+1, frt)
            elif self.material_type == 'Hardening':
                ops.uniaxialMaterial('Hardening', matTag+1,
                                     self.E, self.Fy, 0.0, self.Hk)
                ops.uniaxialMaterial('InitStressMaterial',
                                     matTag, matTag+1, frt)
            else:
                raise Exception(
                    'Input Error - unknown material type (%s)' % self.material_type)
            ops.patch('rect', matTag, Nfw, 1, -self.dw() /
                      2, -self.tw/2, self.dw()/2, self.tw/2)

            # Define flange fibers
            region_width = self.bf/self.num_regions
            for i in range(self.num_regions):
                fri = self.frc + ((i+0.5)/self.num_regions)*(frt-self.frc)

                matTagi = matTag+2*(i+1)
                if self.material_type == 'ElasticPP':
                    ops.uniaxialMaterial(
                        'ElasticPP', matTagi, self.E, self.Fy/self.E, -self.Fy/self.E, fri/self.E)
                elif self.material_type == 'Steel01':
                    b = self.Hk/(self.E+self.Hk)
                    ops.uniaxialMaterial(
                        'Steel01', matTagi+1, self.Fy, self.E, b)
                    ops.uniaxialMaterial(
                        'InitStressMaterial', matTagi, matTagi+1, fri)
                elif self.material_type == 'Hardening':
                    ops.uniaxialMaterial(
                        'Hardening', matTagi+1, self.E, self.Fy, 0.0, self.Hk)
                    ops.uniaxialMaterial(
                        'InitStressMaterial', matTagi, matTagi+1, fri)
                else:
                    raise Exception(
                        'Input Error - unknown material type (%s)' % self.material_type)

                ops.patch('rect', matTagi, Nff, 1, self.dw()/2, -
                          region_width/2, self.d/2, region_width/2)
                ops.patch('rect', matTagi, Nff, 1, -self.d/2, -
                          region_width/2, -self.dw()/2, region_width/2)
        return

    def perform_OpenSees_analysis(self):
        mid_node = int(self.num_elements/2)

        # set modelbuilder
        ops.model('basic', '-ndm', 2, '-ndf', 3)

        # create nodes
        for i in range(self.num_elements+1):
            x_over_L = float(i)/(self.num_elements)
            ix = self.L*x_over_L
            iy = self.zi+(self.zj-self.zi)*x_over_L
            if self.c != 0:
                r = (self.c**2 + (self.L/2)**2)/(2*self.c)
                iy += sqrt(r**2 - (self.L*x_over_L - self.L/2)**2) + (self.c-r)
            ops.node(i, ix, iy)

        # set boundary condition
        ops.fix(0, 1, 1, 0)
        ops.fix(self.num_elements, 0, 1, 0)

        # define coordinate transformation
        ops.geomTransf(self.geomTransfType, 1)

        # define cross section
        self.define_fiber_section(1, 1)
        ops.beamIntegration('Lobatto', 1, 1, 3)

        # Time series for loads
        ops.timeSeries("Constant", 1)

        # define elements
        for i in range(0, self.num_elements):
            # ops.element("forceBeamColumn",i,i,i+1,1,1)
            ops.element("dispBeamColumn", i, i, i+1, 1, 1)

        # Dead load
        ops.pattern('Plain', -1, 1)
        for i in range(0, self.num_elements):
            ops.eleLoad('-ele', i, '-type', 'beamUniform', -self.wD)
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

        ops.test('NormDispIncr', 1e-6, 100, 1)

        # create analysis object
        ops.analysis("Static")

        # Run dead load analysis
        ops.analyze(1)
        zo = float('inf')
        for i in range(self.num_elements+1):
            iz = ops.nodeCoord(i, 2) + ops.nodeDisp(i, 2)
            if iz < zo:
                zo = iz
            print(f"node {i}: {ops.nodeDisp(i, 2)}")

        exit()
        
        # Initilize data        
        data_volume = np.zeros((self.num_steps+1,1))
        data_height = np.full((self.num_steps+1,1),zo)
        end_step = self.num_steps+1        
        
        
        # ------------------------------
        # Finally perform the analysis
        # ------------------------------

        # define ponding load cells    
        PondingLoadCells = dict()
        for i in range(0,self.num_elements):
            PondingLoadCells[i] = PondingLoadCell2d_OPS(i,i,i+1,self.gamma,self.TW)
        
        # Create dict of each node that can have ponding load applied and initilize load to zero
        EmptyPondingLoad = dict()
        for iCell in PondingLoadCells:
            if not PondingLoadCells[iCell].nodeI in EmptyPondingLoad:
                EmptyPondingLoad[PondingLoadCells[iCell].nodeI] = 0.0    
            if not PondingLoadCells[iCell].nodeJ in EmptyPondingLoad:
                EmptyPondingLoad[PondingLoadCells[iCell].nodeJ] = 0.0

        # Perform analysis, ramping up volume      
        zw = zo+0.1
        CurrentPondingLoad = EmptyPondingLoad.copy()
        for iStep in range(0,self.num_steps):
            
            target_volume = (iStep+1)/self.num_steps*self.max_volume
            
            # Update ponding load cells
            for iCell in PondingLoadCells:
                 PondingLoadCells[iCell].update()
            
            # Estimate water height
            for i in range(self.nsteps_vol):
                V = 0
                dVdz = 0
                for iCell in PondingLoadCells:
                    (iV,idVdz) = PondingLoadCells[iCell].get_volume(zw)
                    V += iV
                    dVdz += idVdz
                zw = zw - (V-target_volume)/dVdz
                if abs(target_volume-V) <= self.vol_tol:
                    break 
            
            # Compute load vector
            UpdatedPondingLoad = EmptyPondingLoad.copy()
            for iCell in PondingLoadCells:    
                f = PondingLoadCells[iCell].get_load_vector(zw)
                UpdatedPondingLoad[PondingLoadCells[iCell].nodeI] += f.item(0)
                UpdatedPondingLoad[PondingLoadCells[iCell].nodeJ] += f.item(1)
                
            # Apply difference to model
            ops.pattern("Plain", iStep, 1)
            for iNode in UpdatedPondingLoad:        
                fy = UpdatedPondingLoad[iNode] - CurrentPondingLoad[iNode]
                ops.load(iNode, 0.0, fy, 0.0)
            CurrentPondingLoad = UpdatedPondingLoad

            # Run analysis
            ops.analyze(1)
            
            # Store Data
            data_volume[iStep+1] = target_volume
            data_height[iStep+1] = zw

            # Stop analysis if water level too low
            if (zw-zo) <= (1-self.percent_drop)*(np.amax(data_height)-zo):
                end_step = iStep+1
                break
            
        # Extra Results
        if self.extra_results:
            x = np.linspace(0.0, self.L, num=(self.num_elements+1))
            y = np.zeros((self.num_elements+1,1))
            M = np.zeros((self.num_elements+1,1))
            for i in range(self.num_elements+1):
                y[i] = ops.nodeDisp(i,2)
                if i != self.num_elements:
                    M[i] = -ops.eleForce(i,3)
            extra_results = (x,y,M)
            
        # Wipe Analysis
        ops.wipe()
    
        data_volume = data_volume[:end_step]
        data_height = data_height[:end_step]
    
    
        if self.extra_results:
            return (data_volume,data_height,extra_results)
        else:
            return (data_volume,data_height)

    def steel_beam_object(self):
        elastic_beam = steel_beam();
        elastic_beam.L   = self.L
        elastic_beam.tw  = self.TW
        elastic_beam.zi  = self.zi
        elastic_beam.zj  = self.zj
        elastic_beam.c   = self.c

        elastic_beam.E   = self.E
        elastic_beam.A   = self.A()
        elastic_beam.I   = self.Iz()

        elastic_beam.Mc  = self.Fy*self.Zz()
        # elastic_beam.Vc  = 0.6*self.Fy*self.d*self.tw   @todo add this?

        elastic_beam.alpha   = 1.0
        elastic_beam.LF_D    = 1.0
        elastic_beam.wd      = self.wD/self.TW
        elastic_beam.LF_P    = 1.0
        elastic_beam.gamma   = self.gamma
        elastic_beam.LF_S1   = 0.0
        elastic_beam.LF_S2   = 0.0
        elastic_beam.gammas  = 0.0
        elastic_beam.hs      = 0.0
        elastic_beam.BuildModel();
        return elastic_beam
     

    def maximum_permitted_zw(self,method):
        if method == 'AISC Appendix 2':
            Cs = (self.gamma*self.TW*self.L**4)/(pi**4*self.E*self.Iz())
            # Us_limit = Cs/(1-Cs)
            # fo_limit = 0.8*self.Fy/(Us_limit+1)
            # Mo_limit = fo_limit*self.S
            Mo_limit = 0.8*self.Fy*self.Sz()*(1-Cs)
            
            elastic_beam = self.steel_beam_object()
            elastic_beam.nele   = 40
            elastic_beam.include_ponding_effect = False
            elastic_beam.Mc     = Mo_limit
            zw1 = elastic_beam.Run_To_Strength_Limit()

            elastic_beam.LF_D   = 1.2
            elastic_beam.LF_P   = 1.6
            elastic_beam.Mc     = 0.9*self.Fy*self.Zz()
            zw2 = elastic_beam.Run_To_Strength_Limit()
            
            zw = min(zw1,zw2)
            
        elif method == 'DAMP':
            elastic_beam = self.steel_beam_object()
            elastic_beam.nele   = 40
            elastic_beam.E      = 0.8*self.E
            elastic_beam.LF_D   = 1.2
            elastic_beam.LF_P   = 1.2
            elastic_beam.Mc     = 0.9*self.Fy*self.Zz()
            zw = elastic_beam.Run_To_Strength_Limit()

        elif method == 'DAMP without 0.8':
            elastic_beam = self.steel_beam_object()
            elastic_beam.nele   = 40
            elastic_beam.LF_D   = 1.2
            elastic_beam.LF_P   = 1.2
            elastic_beam.Mc     = 0.9*self.Fy*self.Zz()
            zw = elastic_beam.Run_To_Strength_Limit()
            
        elif method == 'DAMP 1.4':
            elastic_beam = self.steel_beam_object()
            elastic_beam.nele   = 40
            elastic_beam.E      = 0.8*self.E
            elastic_beam.LF_D   = 1.2
            elastic_beam.LF_P   = 1.4
            elastic_beam.Mc     = 0.9*self.Fy*self.Zz()
            zw = elastic_beam.Run_To_Strength_Limit()

        elif method == 'DAMP 1.6':
            elastic_beam = self.steel_beam_object()
            elastic_beam.nele   = 40
            elastic_beam.E      = 0.8*self.E
            elastic_beam.LF_D   = 1.2
            elastic_beam.LF_P   = 1.6
            elastic_beam.Mc     = 0.9*self.Fy*self.Zz()
            zw = elastic_beam.Run_To_Strength_Limit()            
            
        elif method == 'Proposed for ASCE 7' or method == 'Modified Rain Load':
            elastic_beam = self.steel_beam_object()
            elastic_beam.nele   = 40
            elastic_beam.LF_D   = 1.2
            elastic_beam.LF_P   = 1.6
            elastic_beam.Mc     = 0.9*self.Fy*self.Zz()
            elastic_beam.modified_rain_load = True
            zw = elastic_beam.Run_To_Strength_Limit()

        elif method == 'Modified Rain Load with 0.8':
            elastic_beam = self.steel_beam_object()
            elastic_beam.nele   = 40
            elastic_beam.E      = 0.8*self.E
            elastic_beam.LF_D   = 1.2
            elastic_beam.LF_P   = 1.6
            elastic_beam.Mc     = 0.9*self.Fy*self.Zz()
            elastic_beam.modified_rain_load = True
            zw = elastic_beam.Run_To_Strength_Limit()
            
        elif method == 'Neglect Ponding':
            elastic_beam = self.steel_beam_object()
            elastic_beam.nele   = 40
            elastic_beam.LF_D   = 1.2
            elastic_beam.LF_P   = 1.6
            elastic_beam.Mc     = 0.9*self.Fy*self.Zz()
            elastic_beam.include_ponding_effect = False
            zw = elastic_beam.Run_To_Strength_Limit()
            
        else:
            print('Unknown method: %s' % method)
            zw = float('nan')
            
        return zw
