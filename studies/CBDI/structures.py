import numpy as np
import matplotlib.pyplot as plt
import time
from math import atan2,pi,pow,sqrt,ceil
from mpl_toolkits.mplot3d import Axes3D
from PyPonding import PondingLoadManager2d,PondingLoadManager3d
from PyPonding import opensees as ops

def camber(xi,L,c):
    if c == 0:
        z = 0.
    else:
        r = (c**2 + (L/2)**2)/(2*c)
        z = sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
    return z

class AnalysisResults:
    print_each_analysis_time_increment = True
    total_analysis_time = 0.
    
    def __init__(self):
        pass
        
    def add_to_analysis_time(self,tic,toc):
        self.total_analysis_time += toc - tic
        if self.print_each_analysis_time_increment:
            print(f"Adding {toc - tic:0.4f} seconds to total analysis run time.")

    def print_total_analysis_time(self):
        print(f"Total time for ops.analyze: {self.total_analysis_time:0.4f} seconds")

class ExampleStructure:

    # Common Properties
    include_ponding_effect = True
    print_each_analysis_time_increment = False
    test_flag = 0        # Print flag for OpenSees test object
    
    def __init__(self):
        pass

class ExampleBeam(ExampleStructure):

    # Geometric Properties
    yi  = 0.
    yj  = 0.
    c   = 0.
    
    # Support Spring Stiffnessnes
    yi_fixed = True
    yj_fixed = True
    xj_fixed = False
    kyi = 0.
    kyj = 0.
    kxj = 0.

    # Material Properties 
    material_type = 'Hardening'
    hardening_ratio = 0.001      # Ratio of hardening stiffness to elastic stiffness
    residual_stress_ratio = -0.3 # Ratio of peak compressive residual stress to yield stress (should be negative)
    num_regions_RS = 10          # Number of regions for the discretization of residual stress

    # Cross-Sectional Properties
    _A  = None  # Area override
    _Iz = None  # Moment of inertia override
    
    # Loads
    gamma   = 62.4/1000/12**3

    # Analysis Options
    num_elements = 20
    num_divisions_per_element = 1
    
    transf_type = 'Linear'
    num_fiber = 20          # Nominal number of fibers along depth of section
    percent_drop = 20       # Percent drop in simple step volume analysis to halt analysis
    tol_volume = 0.1        # Tolerance for volume iterations
    max_iter_volume = 30    # Maximum number of volume iterations
    tol_load_level = 0.0001 # Tolerance for force in 'IterativeLevel' analyses
    max_iter_level = 30     # Maximum number of iterations for 'IterativeLevel' analyses
    nIP = 3                 # Number of integration points override
    element_type = 'dispBeamColumn'
    beamint_type = 'Legendre'
    
    # Tags    
    transf_tag  = 1
    section_tag = 1
    beamint_tag = 1
    dead_load_pattern_tag = 1
    dead_load_ts_tag = 1
    ponding_load_pattern_tag_start = 2
    ponding_load_ts_tag = 2  
    
    def __init__(self,d,tw,bf,tf,Fy,E,L,S,qD):
        self.d = d              # Section Depth
        self.tw = tw            # Thickness of the web
        self.bf = bf            # Width of the flange
        self.tf = tf            # Thickness of the flnage
        self.Fy = Fy            # Yield stress
        self.E = E              # Elastic modulus
        self.L = L              # Member length
        self.S = S              # Member spacing (tributary width)
        self.qD = qD            # Uniform dead load (forcer per unit area, downward positive)
 
    def lowest_point(self):
        yo = float('inf')
        for i in ops.getNodeTags():
            iy = ops.nodeCoord(i, 2) + ops.nodeDisp(i, 2)
            if iy < yo:
                yo = iy
        return yo
 
    @property
    def num_divisions(self):
        return self.num_divisions_per_element*self.num_elements

    @property
    def dw(self):
        dw = self.d-2*self.tf
        return dw
    
    @property
    def A(self):
        if self._A is None:
            return 2*self.bf*self.tf + (self.d-2*self.tf)*self.tw
        else: 
            return self._A
    
    @A.setter
    def A(self, value):
        self._A = value
       
    @property
    def Iz(self):
        if self._Iz is None:
            return (1.0/12)*self.bf*self.d**3 - (1.0/12)*(self.bf-self.tw)*self.dw**3
        else:
            return self._Iz

    @Iz.setter
    def Iz(self, value):
        self._Iz = value

    def define_fiber_section(self,secTag,matTag):
        Nfw = ceil(self.dw*(self.num_fiber/self.d))
        Nff = ceil(self.tf*(self.num_fiber/self.d))
        
        Hk = self.hardening_ratio*self.E
        
        if self.residual_stress_ratio == 0 or self.material_type == 'Elastic':
            if self.material_type == 'Elastic':
                ops.uniaxialMaterial('Elastic', matTag, self.E)
            elif self.material_type == 'ElasticPP':
                ops.uniaxialMaterial('ElasticPP', matTag, self.E, self.Fy/self.E)
            elif self.material_type == 'Steel01':
                b = Hk/(self.E+Hk)
                ops.uniaxialMaterial('Steel01', matTag, self.Fy, self.E, b)
            elif self.material_type == 'Hardening':
                ops.uniaxialMaterial('Hardening', matTag, self.E, self.Fy, 0.0, Hk)
            else:
                raise Exception('Input Error - unknown material type (%s)' % self.material_type)
            
            ops.section('WFSection2d',secTag,matTag,self.d,self.tw,self.bf,self.tf,Nfw,Nff)
        else: 
            ops.section('Fiber', secTag)    
        
            frc = self.residual_stress_ratio*self.Fy
            frt = -frc*(self.bf*self.tf)/(self.bf*self.tf+self.tw*self.dw)
            
            # Define web fibers
            if self.material_type == 'ElasticPP':
                ops.uniaxialMaterial('ElasticPP', matTag, self.E, self.Fy/self.E, -self.Fy/self.E, frt/self.E)
            elif self.material_type == 'Steel01':
                b = Hk/(self.E+Hk)
                ops.uniaxialMaterial('Steel01', matTag+1, self.Fy, self.E, b)
                ops.uniaxialMaterial('InitStressMaterial', matTag , matTag+1, frt)
            elif self.material_type == 'Hardening':
                ops.uniaxialMaterial('Hardening', matTag+1, self.E, self.Fy, 0.0, Hk)
                ops.uniaxialMaterial('InitStressMaterial', matTag , matTag+1, frt)
            else:
                raise Exception('Input Error - unknown material type (%s)' % self.material_type)
            ops.patch('rect', matTag, Nfw, 1, -self.dw/2, -self.tw/2, self.dw/2, self.tw/2)
                      
            # Define flange fibers
            region_width = self.bf/self.num_regions_RS
            for i in range(self.num_regions_RS):
                fri = frc + ((i+0.5)/self.num_regions_RS)*(frt-frc)
            
                matTagi = matTag+2*(i+1)
                if self.material_type == 'ElasticPP':
                    ops.uniaxialMaterial('ElasticPP', matTagi, self.E, self.Fy/self.E, -self.Fy/self.E, fri/self.E)
                elif self.material_type == 'Steel01':
                    b = Hk/(self.E+Hk)
                    ops.uniaxialMaterial('Steel01', matTagi+1, self.Fy, self.E, b)
                    ops.uniaxialMaterial('InitStressMaterial', matTagi , matTagi+1, fri)
                elif self.material_type == 'Hardening':
                    ops.uniaxialMaterial('Hardening', matTagi+1, self.E, self.Fy, 0.0, Hk)
                    ops.uniaxialMaterial('InitStressMaterial', matTagi , matTagi+1, fri)
                else:
                    raise Exception('Input Error - unknown material type (%s)' % self.material_type)
            
                ops.patch('rect', matTagi, Nff, 1, self.dw/2, -region_width/2,   self.d/2, region_width/2)
                ops.patch('rect', matTagi, Nff, 1, -self.d/2, -region_width/2, -self.dw/2, region_width/2)
        return

    def RunAnalysis(self,analysis_type,target_zw=None,target_Vw=None,num_steps=None):

        if analysis_type.lower() == 'simplesteplevel':
            # Path analysis, ramping up level using a simple step incremental procedure
            if num_steps == None:
                raise Exception('num_steps required for simple step level analysis')
            if target_zw == None:
                raise Exception('target_zw required for simple step level analysis')
     
        elif analysis_type.lower() == 'simplestepvolume':
            # Path analysis, ramping up volume using a simple step incremental procedure
            if num_steps == None:
                raise Exception('num_steps required for simple step volume analysis')
            if target_Vw == None:
                raise Exception('target_Vw required for simple step volume analysis')        
             
        elif analysis_type.lower() == 'iterativelevel':
            # Lumped analysis, going directly to zw and iterating
            if target_zw == None:
                raise Exception('target_zw required for simple step level analysis')
            
        else:
            raise Exception('Unknown analysis type: %s' % analysis_type)

        # Create OpenSees model
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)

        # Define nodes
        for i in range(self.num_elements+1):
            xi = (i/self.num_elements)
            x = xi*self.L
            y = self.yi + xi*(self.yj-self.yi) + camber(xi,self.L,self.c)
            ops.node(i,x,y)
        
        # Define boundary conditions
        if self.yi_fixed:
            ops.fix(0,1,1,0)
        else:
            ops.fix(0,1,0,0)
            ops.node(1001,0.0,self.yi)
            ops.fix(1001,1,1,1)
            ops.uniaxialMaterial('Elastic', 1001, self.kyi)
            ops.element('zeroLength', 1001, 1001, 0, '-mat', 1001, '-dir', 2)

        if self.yj_fixed and self.xj_fixed:
            ops.fix(self.num_elements,1,1,0)
        elif self.yj_fixed:
            ops.fix(self.num_elements,0,1,0)
        elif self.xj_fixed:
            ops.fix(self.num_elements,1,0,0)
        
        if (not self.yj_fixed) and (self.kyj != 0.):
            ops.node(1002,self.L,self.yj)
            ops.fix(1002,1,1,1)
            ops.uniaxialMaterial('Elastic', 1002, self.kyj)
            ops.element('zeroLength', 1002, 1002, self.num_elements, '-mat', 1002, '-dir', 2)
        
        if (not self.xj_fixed) and (self.kxj != 0.):
            ops.node(1003,self.L,self.yj)
            ops.fix(1003,1,1,1)
            ops.uniaxialMaterial('Elastic', 1003, self.kxj)
            ops.element('zeroLength', 1003, 1003, self.num_elements, '-mat', 1003, '-dir', 1)
            
        # Define elements
        ops.geomTransf(self.transf_type, self.transf_tag)
        if self.material_type.lower() == 'elasticsection':
            ops.section('Elastic', self.section_tag, self.E, self.A, self.Iz)
        else:
            self.define_fiber_section(self.section_tag,1)
        if self.beamint_type == 'OpenNewtonCotes':
            xi_list = [None]*self.nIP
            for i in range(self.nIP):
                xi_list[i] = (2*i+1)/(2*self.nIP)
            sec_tag_list = [self.section_tag]*self.nIP
            ops.beamIntegration('FixedLocation', self.beamint_tag, self.nIP, *sec_tag_list, *xi_list)
        elif self.beamint_type == 'MidDistance':
            xi_list = [None]*self.nIP
            for i in range(self.nIP):
                xi_list[i] = (2*i+1)/(2*self.nIP)
            sec_tag_list = [self.section_tag]*self.nIP
            ops.beamIntegration('MidDistance', self.beamint_tag, self.nIP, *sec_tag_list, *xi_list)
        else:
            ops.beamIntegration(self.beamint_type, self.beamint_tag, self.section_tag, self.nIP)
        for i in range(self.num_elements):
            ops.element(self.element_type, i, i, i+1, self.transf_tag, self.beamint_tag)

        # Define Ponding Load Cells
        PondingLoadManager = PondingLoadManager2d()
        for i in range(self.num_divisions):
            y_offsetI = 0.
            y_offsetJ = 0.

            # I end
            if i % self.num_divisions_per_element == 0:
                n = i//self.num_divisions_per_element
                endI = ('node',n)
            else:
                e = i//self.num_divisions_per_element
                x = (i%self.num_divisions_per_element)/self.num_divisions_per_element
                endI = ('element',e,x)
                if self.num_elements == 1 or self.c == 0:
                    y_offsetI = camber(i/self.num_divisions,self.L,self.c)
                else:
                    raise Exception('Camber not yet implemnted for multiple elements')
            
            # J end
            if (i+1) % self.num_divisions_per_element == 0:
                n = (i+1)//self.num_divisions_per_element
                endJ = ('node',n)
            else:
                e = (i+1)//self.num_divisions_per_element
                x = ((i+1)%self.num_divisions_per_element)/self.num_divisions_per_element
                endJ = ('element',e,x)
                if self.num_elements == 1 or self.c == 0:
                    y_offsetJ = camber((i+1)/self.num_divisions,self.L,self.c)
                else:
                    raise Exception('Camber not yet implemnted for multiple elements')
            
            # Add load cell
            PondingLoadManager.add_cell(i,endI,endJ,self.gamma,self.S)
            PondingLoadManager.cells[i].endI.y_offset = y_offsetI 
            PondingLoadManager.cells[i].endJ.y_offset = y_offsetJ
            PondingLoadManager.cells[i].update_coord()
                
        # Define Dead Load
        ops.timeSeries("Constant", self.dead_load_ts_tag)
        ops.pattern("Plain", self.dead_load_pattern_tag, self.dead_load_ts_tag)

        # Define uniform dead load on joists
        for i in range(self.num_elements):
            nodes = ops.eleNodes(i)
            coordi = ops.nodeCoord(nodes[0])
            coordj = ops.nodeCoord(nodes[1])
            Lx = coordj[0]-coordi[0]
            Ly = coordj[1]-coordi[1]
            L = sqrt(Lx**2 + Ly**2)
            Wy = -self.qD*self.S*(Lx/L)*(Lx/L)
            Wx = -self.qD*self.S*(Lx/L)*(Ly/L)
            ops.eleLoad('-ele', i, '-type', '-beamUniform', Wy, Wx)

        # Initilize data
        results = AnalysisResults()
        results.print_each_analysis_time_increment = self.print_each_analysis_time_increment
        results.analysis_type = analysis_type
        
        # Set up analysis
        ops.numberer("RCM")
        ops.constraints("Plain")
        ops.system("BandSPD")
        ops.test("NormUnbalance", 1.0e-6, 25, self.test_flag)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        tic = time.perf_counter()
        ops.analyze(1)
        toc = time.perf_counter()
        results.add_to_analysis_time(tic,toc)
        ops.reactions()

        # Run Ponding Analysis
        ops.timeSeries("Constant", self.ponding_load_ts_tag)
        
        if analysis_type.lower() == 'simplesteplevel':

            # Initialize results
            results.target_zw = target_zw
            results.water_volume = np.zeros((num_steps+1,1))
            results.water_level = np.zeros((num_steps+1,1))
            results.Rxi = np.zeros((num_steps+1,1))
            results.Ryi = np.zeros((num_steps+1,1))
            results.Ryj = np.zeros((num_steps+1,1))  

            # Find lowest point
            yo = self.lowest_point()

            # Store dead load results
            results.water_volume[0] = 0.
            results.water_level[0] = yo
            results.Rxi[0] = self.Rxi()
            results.Ryi[0] = self.Ryi()
            results.Ryj[0] = self.Ryj()
        
            for iStep in range(1,num_steps+1):

                # Update ponding load cells
                if self.include_ponding_effect:
                    PondingLoadManager.update()

                # Compute load vector
                izw = yo + (iStep/num_steps)*(target_zw-yo)
                (iV,idVdz) = PondingLoadManager.get_volume(izw)
                PondingLoadManager.compute_current_load_vector(izw)

                # Apply difference to model
                ops.pattern("Plain", self.ponding_load_pattern_tag_start+iStep, self.ponding_load_ts_tag)
                PondingLoadManager.apply_load_increment()
                PondingLoadManager.commit_current_load_vector()

                # Run analysis
                tic = time.perf_counter()
                ops.analyze(1)
                toc = time.perf_counter()
                results.add_to_analysis_time(tic,toc)
                ops.reactions()

                # Store Reuslts
                results.water_volume[iStep] = iV
                results.water_level[iStep] = izw
                results.Rxi[iStep] = self.Rxi()
                results.Ryi[iStep] = self.Ryi()
                results.Ryj[iStep] = self.Ryj()
                
        elif analysis_type.lower() == 'simplestepvolume':
        
            # Initialize results
            results.target_Vw = target_Vw
            results.water_volume = np.zeros((num_steps+1,1))
            results.water_level = np.zeros((num_steps+1,1))
            results.Rxi = np.zeros((num_steps+1,1))
            results.Ryi = np.zeros((num_steps+1,1))
            results.Ryj = np.zeros((num_steps+1,1))  

            # Find lowest point
            yo = self.lowest_point()
            
            # Store dead load results
            results.water_volume[0] = 0.
            results.water_level[0] = yo
            results.Rxi[0] = self.Rxi()
            results.Ryi[0] = self.Ryi()
            results.Ryj[0] = self.Ryj()
        
            for iStep in range(1,num_steps+1):

                # Update ponding load cells
                if self.include_ponding_effect:
                    PondingLoadManager.update()

                # Estimate water height
                step_Vw = (iStep+1)/num_steps*target_Vw
                if iStep == 1:
                    izw = yo+0.1 # Initial guess
                for i in range(self.max_iter_volume):
                    (iV,idVdz) = PondingLoadManager.get_volume(izw)
                    izw = izw - (iV-step_Vw)/idVdz
                    if abs(iV-step_Vw) <= self.tol_volume:
                        break 

                # Compute load vector
                PondingLoadManager.compute_current_load_vector(izw)

                # Apply difference to model
                ops.pattern("Plain", self.ponding_load_pattern_tag_start+iStep, self.ponding_load_ts_tag)
                PondingLoadManager.apply_load_increment()
                PondingLoadManager.commit_current_load_vector()

                # Run analysis
                tic = time.perf_counter()
                ops.analyze(1)
                toc = time.perf_counter()
                results.add_to_analysis_time(tic,toc)
                ops.reactions()

                # Store Reuslts
                results.water_volume[iStep] = step_Vw
                results.water_level[iStep] = izw
                results.Rxi[iStep] = self.Rxi()
                results.Ryi[iStep] = self.Ryi()
                results.Ryj[iStep] = self.Ryj()
            
                # Stop analysis if water level too low
                if (izw-yo) <= (1-0.01*self.percent_drop)*(np.amax(results.water_level)-yo):
                    # Truncate remainging results and break out of stepping loop
                    results.water_volume = results.water_volume[:(iStep+1)]
                    results.water_level = results.water_level[:(iStep+1)]
                    results.Rxi = results.Rxi[:(iStep+1)]
                    results.Ryi = results.Ryi[:(iStep+1)]
                    results.Ryj = results.Ryj[:(iStep+1)]
                    break            
            
        elif analysis_type.lower() == 'iterativelevel':

            for iStep in range(1,self.max_iter_level+1):

                # Update ponding load cells
                if self.include_ponding_effect:
                    PondingLoadManager.update()

                # Compute load vector
                PondingLoadManager.compute_current_load_vector(target_zw)

                # Check for convergence
                if PondingLoadManager.sub_abs_diff_load_increment() < self.tol_load_level:
                    print('Converged')
                    break
                    
                # Print data on iteration
                print('Iteration: %3.i, Total Water Load: %0.3f' % (iStep,PondingLoadManager.total_current_load()))
                    
                # Apply difference to model
                ops.pattern("Plain", self.ponding_load_pattern_tag_start+iStep, self.ponding_load_ts_tag)
                PondingLoadManager.apply_load_increment()
                PondingLoadManager.commit_current_load_vector()

                # Run analysis
                tic = time.perf_counter()
                ops.analyze(1)
                toc = time.perf_counter()
                results.add_to_analysis_time(tic,toc)
                ops.reactions()

            # Store Reuslts
            (V,dVdz) = PondingLoadManager.get_volume(target_zw)
            results.water_volume = V
            results.water_level  = target_zw            
            results.Rxi = self.Rxi()
            results.Ryi = self.Ryi()
            results.Ryj = self.Ryj()
            
        else:
            raise Exception('Unknown analysis type: %s' % analysis_type)

        results.print_total_analysis_time()

        return results

    def Rxi(self):
        R = ops.nodeReaction(0, 1)
        return R        
        
    def Ryi(self):
        if self.yi_fixed:
            R = ops.nodeReaction(0, 2)
        else:
            R = ops.nodeReaction(1001, 2)
        return R

    def Ryj(self):
        if self.yj_fixed:
            R = ops.nodeReaction(self.num_elements, 2)
        else:
            R = ops.nodeReaction(1002, 2)
        return R
 
class ExampleRoof(ExampleStructure):

    # Bay Widths
    L_AB = 480
    L_BC = 480
    L_CD = 480
    L_12 = 480
    L_23 = 480

    # Roof Elevations at Columns
    z_A1 =   0
    z_A2 =  -8
    z_A3 =   0
    z_B1 =   0
    z_B2 = -10
    z_B3 =   0
    z_C1 =   0
    z_C2 = -10
    z_C3 =   0
    z_D1 =   0
    z_D2 =  -8
    z_D3 =   0

    # Joist and Joist Girder Section Properties
    E = 29000

    A_J   = 100.
    Iz_J  = 215.
    Iy_J  = 100.
    GJ_J  = 100.

    A_JG  = 100.
    Iz_JG = 2029.
    Iy_JG = 100.
    GJ_JG = 100.

    # Other Properties
    c_J     = 1.
    c_JG    = 1.
    nspaces = 8

    # Loading Parameters
    wdJG    = 50./1000/12 # Joist Girder self-weight (forcer per unit length, downward positive)
    qD      = 10./1000/12**2 # Uniform dead load (forcer per unit area, downward positive)
    gamma   = 62.4/1000/12**3

    # Analysis Options
    use_CBDI = False    
    ndiv_J  = 10
    na      = 4
    nb      = 4
    tol_volume = 0.1        # Tolerance for volume iterations
    max_iter_volume = 30    # Maximum number of volume iterations
    tol_load_level = 0.0001 # Tolerance for force in 'IterativeLevel' analyses
    max_iter_level = 30     # Maximum number of iterations for 'IterativeLevel' analyses
    _nIP_J = None           # Number of integration points override
    nIP_JG = 4              # Number of integration points for Joist Girder
    _element_type_J_override = None
    element_type_JG = 'dispBeamColumn'
    
    # Analysis Tags
    girder_transf_tag   = 1
    girder_section_tag  = 1
    girder_beamint_tag  = 1

    joist_transf_tag    = 2
    joist_section_tag   = 2
    joist_beamint_tag   = 2

    dead_load_pattern_tag = 1
    dead_load_ts_tag      = 1

    ponding_load_pattern_tag_start = 2
    ponding_load_ts_tag = 2

    # Run Options
    plot_load_cells = False

    def __init__(self):
        pass

    def lowest_point(self):
        zo = float('inf')
        for i in ops.getNodeTags():
            iz = ops.nodeCoord(i, 3) + ops.nodeDisp(i, 3)
            if iz < zo:
                zo = iz
        return zo

    @property
    def nele_J(self):
        if self.use_CBDI:
            return 1
        else:
            return self.ndiv_J

    @property
    def nIP_J(self):
        if self._nIP_J is None:
            if self.use_CBDI:
                return 8
            else:
                return 4
        else:
            return self._nIP_J

    @nIP_J.setter
    def nIP_J(self, value):
        self._nIP_J = value
        
    @property
    def element_type_J(self):
        if self._element_type_J_override is None:
            if self.use_CBDI:
                return 'forceBeamColumn'
            else:
                return 'dispBeamColumn'
        else:
            return self._element_type_J_override

    @element_type_J.setter
    def element_type_J(self, value):
        self._element_type_J_override = value        
        
    def RunAnalysis(self,analysis_type,target_zw=None,target_Vw=None,num_steps=None):

        if analysis_type.lower() == 'simplesteplevel':
            # Path analysis, ramping up level using a simple step incremental procedure
            if num_steps == None:
                raise Exception('num_steps required for simple step level analysis')
            if target_zw == None:
                raise Exception('target_zw required for simple step level analysis')
     
        elif analysis_type.lower() == 'simplestepvolume':
            # Path analysis, ramping up volume using a simple step incremental procedure
            if num_steps == None:
                raise Exception('num_steps required for simple step volume analysis')
            if target_Vw == None:
                raise Exception('target_Vw required for simple step volume analysis')        
             
        elif analysis_type.lower() == 'iterativelevel':
            # Lumped analysis, going directly to zw and iterating
            if target_zw == None:
                raise Exception('target_zw required for simple step level analysis')
            
        else:
            raise Exception('Unknown analysis type: %s' % analysis_type)

        # Create OpenSees model
        ops.wipe()
        ops.model('basic', '-ndm', 3, '-ndf', 6)

        ###########################################################
        # Define Joist Girders
        ###########################################################
        ops.geomTransf('Linear', self.girder_transf_tag, 0.0, -1.0, 0.0)
        ops.section('Elastic', self.girder_section_tag, self.E, self.A_JG, self.Iz_JG, self.Iy_JG, 1, self.GJ_JG)
        ops.beamIntegration('Lobatto', self.girder_beamint_tag, self.girder_section_tag, self.nIP_JG)
        
        z_position = dict()
        
        # Define Joist Girder B12
        for i in range(self.nspaces+1):
            n = 110101+i

            xi = (i/self.nspaces)
            L = self.L_12
            x = xi*L
            y = self.L_AB
            z = self.z_B1 + xi*(self.z_B2-self.z_B1) + camber(xi,L,self.c_JG)
            z_position[n] = z

            ops.node(n,x,y,z)
            if i == 0:
                ops.fix(n,1,1,1,1,0,1)
            elif i == self.nspaces:
                ops.fix(n,0,1,1,1,0,1)
            else:
                ops.fix(n,0,1,0,1,0,1)

        for i in range(self.nspaces):
            ops.element(self.element_type_JG, 110101+i, 110101+i, 110102+i, self.girder_transf_tag, self.girder_beamint_tag)

        # Define Joist Girder B23
        for i in range(self.nspaces+1):
            n = 110201+i

            xi = (i/self.nspaces)
            L = self.L_23
            x = self.L_12 + xi*L
            y = self.L_AB
            z = self.z_B2 + xi*(self.z_B3-self.z_B2) + camber(xi,L,self.c_JG)
            z_position[n] = z

            ops.node(n,x,y,z)
            if i == 0:
                ops.fix(n,1,1,1,1,0,1)
            elif i == self.nspaces:
                ops.fix(n,0,1,1,1,0,1)
            else:
                ops.fix(n,0,1,0,1,0,1)

        for i in range(self.nspaces):
            ops.element(self.element_type_JG, 110201+i, 110201+i, 110202+i, self.girder_transf_tag, self.girder_beamint_tag)

        # Define Joist Girder C12
        for i in range(self.nspaces+1):
            n = 120101+i

            xi = (i/self.nspaces)
            L = self.L_12
            x = xi*L
            y = self.L_AB + self.L_BC
            z = self.z_C1 + xi*(self.z_C2-self.z_C1) + camber(xi,L,self.c_JG)
            z_position[n] = z

            ops.node(n,x,y,z)
            if i == 0:
                ops.fix(n,1,1,1,1,0,1)
            elif i == self.nspaces:
                ops.fix(n,0,1,1,1,0,1)
            else:
                ops.fix(n,0,1,0,1,0,1)

        for i in range(self.nspaces):
            ops.element(self.element_type_JG, 120101+i, 120101+i, 120102+i, self.girder_transf_tag, self.girder_beamint_tag)

        # Define Joist Girder C23
        for i in range(self.nspaces+1):
            n = 120201+i

            xi = (i/self.nspaces)
            L = self.L_23
            x = self.L_12 + xi*L
            y = self.L_AB + self.L_BC
            z = self.z_C2 + xi*(self.z_C3-self.z_C2) + camber(xi,L,self.c_JG)
            z_position[n] = z

            ops.node(n,x,y,z)
            if i == 0:
                ops.fix(n,1,1,1,1,0,1)
            elif i == self.nspaces:
                ops.fix(n,0,1,1,1,0,1)
            else:
                ops.fix(n,0,1,0,1,0,1)

        for i in range(self.nspaces):
            ops.element(self.element_type_JG, 120201+i, 120201+i, 120202+i, self.girder_transf_tag, self.girder_beamint_tag)


        ###########################################################
        # Define Joists
        ##########################################################
        ops.geomTransf('Linear', self.joist_transf_tag, 1.0, 0.0, 0.0)
        ops.section('Elastic', self.joist_section_tag, self.E, self.A_J, self.Iz_J, self.Iy_J, 1, self.GJ_J)
        ops.beamIntegration('Lobatto', self.joist_beamint_tag, self.joist_section_tag, self.nIP_J)
        
        # Define joists between grid lines A and B
        for i in range(1,2*self.nspaces):

            if i < self.nspaces:
                nj = 110101+i
                zi = self.z_A1 + (i/self.nspaces)*(self.z_A2-self.z_A1)
                zj = z_position[nj]
                x = (i/self.nspaces)*self.L_12
            elif i == self.nspaces:
                zi = self.z_A2
                zj = self.z_B2
                x = self.L_12
            else:
                nj = 110201+(i-self.nspaces)
                zi = self.z_A2 + ((i-self.nspaces)/self.nspaces)*(self.z_A3-self.z_A2)
                zj = z_position[nj]
                x = self.L_12 + ((i-self.nspaces)/self.nspaces)*self.L_23


            for j in range(self.nele_J+1):
                n = 210001+100*i+j

                xi = j/self.nele_J
                L = self.L_AB
                y = xi*L
                z = zi + xi*(zj-zi) + camber(xi,L,self.c_J)

                ops.node(n,x,y,z)
                if j == 0:
                    ops.fix(n,1,1,1,0,1,1)
                elif j == self.nele_J:
                    if i == self.nspaces:
                        ops.fix(n,1,0,1,0,1,1)
                    else:
                        ops.fix(n,1,0,0,0,1,1)
                        ops.equalDOF(nj,n,3)
                else:
                    ops.fix(n,1,0,0,0,1,1)

            for j in range(self.nele_J):
                ops.element(self.element_type_J, 210001+100*i+j, 210001+100*i+j, 210002+100*i+j, self.joist_transf_tag, self.joist_beamint_tag)


        # Define joists between grid lines B and C
        for i in range(1,2*self.nspaces):

            if i < self.nspaces:
                ni = 110101+i
                nj = 120101+i
                zi = z_position[ni]
                zj = z_position[nj]
                x = (i/self.nspaces)*self.L_12
            elif i == self.nspaces:
                zi = self.z_B2
                zj = self.z_C2
                x = self.L_12
            else:
                ni = 110201+(i-self.nspaces)
                nj = 120201+(i-self.nspaces)
                zi = z_position[ni]
                zj = z_position[nj]
                x = self.L_12 + ((i-self.nspaces)/self.nspaces)*self.L_23


            for j in range(self.nele_J+1):
                n = 220001+100*i+j

                xi = j/self.nele_J
                L = self.L_BC
                y = self.L_AB + xi*L
                z = zi + xi*(zj-zi) + camber(xi,L,self.c_J)

                ops.node(n,x,y,z)
                if j == 0:
                    if i == self.nspaces:
                        ops.fix(n,1,1,1,0,1,1)
                    else:
                        ops.fix(n,1,1,0,0,1,1)
                        ops.equalDOF(ni,n,3)
                elif j == self.nele_J:
                    if i == self.nspaces:
                        ops.fix(n,1,0,1,0,1,1)
                    else:
                        ops.fix(n,1,0,0,0,1,1)
                        ops.equalDOF(nj,n,3)
                else:
                    ops.fix(n,1,0,0,0,1,1)

            for j in range(self.nele_J):
                ops.element(self.element_type_J, 220001+100*i+j, 220001+100*i+j, 220002+100*i+j, self.joist_transf_tag, self.joist_beamint_tag)


        # Define joists between grid lines C and D
        for i in range(1,2*self.nspaces):

            if i < self.nspaces:
                ni = 120101+i
                zi = z_position[ni]
                zj = self.z_D1 + (i/self.nspaces)*(self.z_D2-self.z_D1)
                x = (i/self.nspaces)*self.L_12
            elif i == self.nspaces:
                zi = self.z_C2
                zj = self.z_D2
                x = self.L_12
            else:
                ni = 120201+(i-self.nspaces)
                zi = z_position[ni]
                zj = self.z_D2 + ((i-self.nspaces)/self.nspaces)*(self.z_D3-self.z_D2)
                x = self.L_12 + ((i-self.nspaces)/self.nspaces)*self.L_23


            for j in range(self.nele_J+1):
                n = 230001+100*i+j

                xi = j/self.nele_J
                L = self.L_CD
                y = self.L_AB + self.L_BC + xi*L
                z = zi + xi*(zj-zi) + camber(xi,L,self.c_J)

                ops.node(n,x,y,z)
                if j == 0:
                    if i == self.nspaces:
                        ops.fix(n,1,1,1,0,1,1)
                    else:
                        ops.fix(n,1,1,0,0,1,1)
                        ops.equalDOF(ni,n,3)
                elif j == self.nele_J:
                    ops.fix(n,1,0,1,0,1,1)
                else:
                    ops.fix(n,1,0,0,0,1,1)

            for j in range(self.nele_J):
                ops.element(self.element_type_J, 230001+100*i+j, 230001+100*i+j, 230002+100*i+j, self.joist_transf_tag, self.joist_beamint_tag)


        ###########################################################
        # Define Ponding Load Cells
        ###########################################################

        PondingLoadManager = PondingLoadManager3d()

        # Define ponding load cells between A and B
        for i in range(1,2*self.nspaces+1):
            for j in range(1,self.ndiv_J+1):
                z_offsetI = 0.
                z_offsetJ = 0.
                z_offsetK = 0.
                z_offsetL = 0.

                # Vertex I
                if i == 1:
                    x = 0.
                    y = ((j-1)/self.ndiv_J)*self.L_AB
                    z = self.z_A1 + ((j-1)/self.ndiv_J)*(self.z_B1-self.z_A1)
                    vertexI = ('fixed',x,y,z)
                else:
                    if j == 1:
                        y = 0.
                        if i <= self.nspaces:
                            x = ((i-1)/self.nspaces)*self.L_12
                            z = self.z_A1 + ((i-1)/self.nspaces)*(self.z_A2-self.z_A1)
                        else:
                            x = self.L_12 + ((i-1-self.nspaces)/self.nspaces)*self.L_23
                            z = self.z_A2 + ((i-1-self.nspaces)/self.nspaces)*(self.z_A3-self.z_A2)
                        vertexI = ('fixed',x,y,z)
                    else:
                        if self.use_CBDI:
                            n = 210001 + 100*(i-1)
                            xi = (j-1)/self.ndiv_J
                            z_offsetI = camber(xi,self.L_AB,self.c_J)
                            vertexI = ('element',n,xi)
                        else:
                            n = 210001 + 100*(i-1) + (j-1)
                            vertexI = ('node',n)

                # Vertex J
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = ((j-1)/self.ndiv_J)*self.L_AB
                    z = self.z_A3 + ((j-1)/self.ndiv_J)*(self.z_B3-self.z_A3)
                    vertexJ = ('fixed',x,y,z)
                else:
                    if j == 1:
                        y = 0.
                        if i <= self.nspaces:
                            x = (i/self.nspaces)*self.L_12
                            z = self.z_A1 + (i/self.nspaces)*(self.z_A2-self.z_A1)
                        else:
                            x = self.L_12 + ((i-self.nspaces)/self.nspaces)*self.L_23
                            z = self.z_A2 + ((i-self.nspaces)/self.nspaces)*(self.z_A3-self.z_A2)
                        vertexJ = ('fixed',x,y,z)
                    else:
                        if self.use_CBDI:
                            n = 210001 + 100*i
                            xi = (j-1)/self.ndiv_J
                            z_offsetJ = camber(xi,self.L_AB,self.c_J)
                            vertexJ = ('element',n,xi)
                        else:
                            n = 210001 + 100*i + (j-1)
                            vertexJ = ('node',n)

                # Vertex K
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = (j/self.ndiv_J)*self.L_AB
                    z = self.z_A3 + (j/self.ndiv_J)*(self.z_B3-self.z_A3)
                    vertexK = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        if i <= self.nspaces:
                            n = 110101 + i
                        elif i == self.nspaces:
                            n = 210001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 110201 + (i-self.nspaces)
                        vertexK = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 210001 + 100*i
                            xi = j/self.ndiv_J
                            z_offsetK = camber(xi,self.L_AB,self.c_J)
                            vertexK = ('element',n,xi)
                        else:
                            n = 210001 + 100*i + j
                            vertexK = ('node',n)

                # Vertex L
                if i == 1:
                    x = 0.
                    y = (j/self.ndiv_J)*self.L_AB
                    z = self.z_A1 + (j/self.ndiv_J)*(self.z_B1-self.z_A1)
                    vertexL = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        if i <= self.nspaces:
                            n = 110101 + (i-1)
                        elif i == self.nspaces+1:
                            n = 210001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 110201 + (i-1-self.nspaces)
                        vertexL = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 210001 + 100*(i-1)
                            xi = j/self.ndiv_J
                            z_offsetL = camber(xi,self.L_AB,self.c_J)
                            vertexL = ('element',n,xi)
                        else:
                            n = 210001 + 100*(i-1) + j
                            vertexL = ('node',n)

                # add load cell
                id = 'AB_%03ix_%03iy' % (i,j)
                PondingLoadManager.add_cell(id,vertexI,vertexJ,vertexK,vertexL,self.gamma,self.na,self.nb)
                PondingLoadManager.cells[id].vertexI.z_offset = z_offsetI
                PondingLoadManager.cells[id].vertexJ.z_offset = z_offsetJ
                PondingLoadManager.cells[id].vertexK.z_offset = z_offsetK
                PondingLoadManager.cells[id].vertexL.z_offset = z_offsetL
                PondingLoadManager.cells[id].update_coord()


        # Define ponding load cells between B and C
        for i in range(1,2*self.nspaces+1):
            for j in range(1,self.ndiv_J+1):
                z_offsetI = 0.
                z_offsetJ = 0.
                z_offsetK = 0.
                z_offsetL = 0.

                # Vertex I
                if i == 1:
                    x = 0.
                    y = self.L_AB + ((j-1)/self.ndiv_J)*self.L_BC
                    z = self.z_B1 + ((j-1)/self.ndiv_J)*(self.z_C1-self.z_B1)
                    vertexI = ('fixed',x,y,z)
                else:
                    if j == 1:
                        if i <= self.nspaces:
                            n = 110101 + (i-1)
                        elif i == self.nspaces+1:
                            n = 210001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 110201 + (i-1-self.nspaces)
                        vertexI = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 220001 + 100*(i-1)
                            xi = (j-1)/self.ndiv_J
                            z_offsetI = camber(xi,self.L_BC,self.c_J)
                            vertexI = ('element',n,xi)
                        else:
                            n = 220001 + 100*(i-1) + (j-1)
                            vertexI = ('node',n)

                # Vertex J
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = self.L_AB + ((j-1)/self.ndiv_J)*self.L_BC
                    z = self.z_B3 + ((j-1)/self.ndiv_J)*(self.z_C3-self.z_B3)
                    vertexJ = ('fixed',x,y,z)
                else:
                    if j == 1:
                        if i <= self.nspaces:
                            n = 110101 + i
                        elif i == self.nspaces:
                            n = 210001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 110201 + (i-self.nspaces)
                        vertexJ = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 220001 + 100*i
                            xi = (j-1)/self.ndiv_J
                            z_offsetJ = camber(xi,self.L_BC,self.c_J)
                            vertexJ = ('element',n,xi)
                        else:
                            n = 220001 + 100*i + (j-1)
                            vertexJ = ('node',n)

                # Vertex K
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = self.L_AB + (j/self.ndiv_J)*self.L_BC
                    z = self.z_B3 + (j/self.ndiv_J)*(self.z_C3-self.z_B3)
                    vertexK = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        if i <= self.nspaces:
                            n = 120101 + i
                        elif i == self.nspaces:
                            n = 220001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 120201 + (i-self.nspaces)
                        vertexK = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 220001 + 100*i
                            xi = j/self.ndiv_J
                            z_offsetK = camber(xi,self.L_BC,self.c_J)
                            vertexK = ('element',n,xi)
                        else:
                            n = 220001 + 100*i + j
                            vertexK = ('node',n)

                # Vertex L
                if i == 1:
                    x = 0.
                    y = self.L_AB + (j/self.ndiv_J)*self.L_BC
                    z = self.z_B1 + (j/self.ndiv_J)*(self.z_C1-self.z_B1)
                    vertexL = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        if i <= self.nspaces:
                            n = 120101 + (i-1)
                        elif i == self.nspaces+1:
                            n = 220001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 120201 + (i-1-self.nspaces)
                        vertexL = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 220001 + 100*(i-1)
                            xi = j/self.ndiv_J
                            z_offsetL = camber(xi,self.L_BC,self.c_J)
                            vertexL = ('element',n,xi)
                        else:
                            n = 220001 + 100*(i-1) + j
                            vertexL = ('node',n)

                # add load cell
                id = 'BC_%03ix_%03iy' % (i,j)
                PondingLoadManager.add_cell(id,vertexI,vertexJ,vertexK,vertexL,self.gamma,self.na,self.nb)
                PondingLoadManager.cells[id].vertexI.z_offset = z_offsetI
                PondingLoadManager.cells[id].vertexJ.z_offset = z_offsetJ
                PondingLoadManager.cells[id].vertexK.z_offset = z_offsetK
                PondingLoadManager.cells[id].vertexL.z_offset = z_offsetL
                PondingLoadManager.cells[id].update_coord()


        # Define ponding load cells between C and D
        for i in range(1,2*self.nspaces+1):
            for j in range(1,self.ndiv_J+1):
                z_offsetI = 0.
                z_offsetJ = 0.
                z_offsetK = 0.
                z_offsetL = 0.

                # Vertex I
                if i == 1:
                    x = 0.
                    y = self.L_AB + self.L_BC + ((j-1)/self.ndiv_J)*self.L_CD
                    z = self.z_C1 + ((j-1)/self.ndiv_J)*(self.z_D1-self.z_C1)
                    vertexI = ('fixed',x,y,z)
                else:
                    if j == 1:
                        if i <= self.nspaces:
                            n = 120101 + (i-1)
                        elif i == self.nspaces+1:
                            n = 220001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 120201 + (i-1-self.nspaces)
                        vertexI = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 230001 + 100*(i-1)
                            xi = (j-1)/self.ndiv_J
                            z_offsetI = camber(xi,self.L_CD,self.c_J)
                            vertexI = ('element',n,xi)
                        else:
                            n = 230001 + 100*(i-1) + (j-1)
                            vertexI = ('node',n)

                # Vertex J
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = self.L_AB + self.L_BC + ((j-1)/self.ndiv_J)*self.L_CD
                    z = self.z_C3 + ((j-1)/self.ndiv_J)*(self.z_D3-self.z_C3)
                    vertexJ = ('fixed',x,y,z)
                else:
                    if j == 1:
                        if i <= self.nspaces:
                            n = 120101 + i
                        elif i == self.nspaces:
                            n = 220001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 120201 + (i-self.nspaces)
                        vertexJ = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 230001 + 100*i
                            xi = (j-1)/self.ndiv_J
                            z_offsetJ = camber(xi,self.L_CD,self.c_J)
                            vertexJ = ('element',n,xi)
                        else:
                            n = 230001 + 100*i + (j-1)
                            vertexJ = ('node',n)

                # Vertex K
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = self.L_AB + self.L_BC + (j/self.ndiv_J)*self.L_CD
                    z = self.z_C3 + (j/self.ndiv_J)*(self.z_D3-self.z_C3)
                    vertexK = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        y = self.L_AB + self.L_BC + self.L_CD
                        if i <= self.nspaces:
                            x = (i/self.nspaces)*self.L_12
                            z = self.z_D1 + (i/self.nspaces)*(self.z_D2-self.z_D1)
                        else:
                            x = self.L_12 + ((i-self.nspaces)/self.nspaces)*self.L_23
                            z = self.z_D2 + ((i-self.nspaces)/self.nspaces)*(self.z_D3-self.z_D2)
                        vertexK = ('fixed',x,y,z)
                    else:
                        if self.use_CBDI:
                            n = 230001 + 100*i
                            xi = j/self.ndiv_J
                            z_offsetK = camber(xi,self.L_CD,self.c_J)
                            vertexK = ('element',n,xi)
                        else:
                            n = 230001 + 100*i + j
                            vertexK = ('node',n)

                # Vertex L
                if i == 1:
                    x = 0.
                    y = self.L_AB + self.L_BC + (j/self.ndiv_J)*self.L_CD
                    z = self.z_C1 + (j/self.ndiv_J)*(self.z_D1-self.z_C1)
                    vertexL = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        y = self.L_AB + self.L_BC + self.L_CD
                        if i <= self.nspaces:
                            x = ((i-1)/self.nspaces)*self.L_12
                            z = self.z_D1 + ((i-1)/self.nspaces)*(self.z_D2-self.z_D1)
                        else:
                            x = self.L_12 + ((i-1-self.nspaces)/self.nspaces)*self.L_23
                            z = self.z_D2 + ((i-1-self.nspaces)/self.nspaces)*(self.z_D3-self.z_D2)
                        vertexL = ('fixed',x,y,z)
                    else:
                        if self.use_CBDI:
                            n = 230001 + 100*(i-1)
                            xi = j/self.ndiv_J
                            z_offsetL = camber(xi,self.L_CD,self.c_J)
                            vertexL = ('element',n,xi)
                        else:
                            n = 230001 + 100*(i-1) + j
                            vertexL = ('node',n)

                # add load cell
                id = 'CD_%03ix_%03iy' % (i,j)
                PondingLoadManager.add_cell(id,vertexI,vertexJ,vertexK,vertexL,self.gamma,self.na,self.nb)
                PondingLoadManager.cells[id].vertexI.z_offset = z_offsetI
                PondingLoadManager.cells[id].vertexJ.z_offset = z_offsetJ
                PondingLoadManager.cells[id].vertexK.z_offset = z_offsetK
                PondingLoadManager.cells[id].vertexL.z_offset = z_offsetL
                PondingLoadManager.cells[id].update_coord()


        ###########################################################
        # Define Dead Load
        ###########################################################

        ops.timeSeries("Constant", self.dead_load_ts_tag)
        ops.pattern("Plain", self.dead_load_pattern_tag, self.dead_load_ts_tag)

        # Define uniform dead load on joists
        for i in range(1,2*self.nspaces):

            if i < self.nspaces:
                wD = self.qD*self.L_12/self.nspaces
            elif i == self.nspaces:
                wD = self.qD*(0.5*self.L_12/self.nspaces+0.5*self.L_23/self.nspaces)
            else:
                wD = self.qD*self.L_23/self.nspaces

            for j in range(self.nele_J):
                n = 210001+100*i+j
                nodes = ops.eleNodes(n)
                coordi = ops.nodeCoord(nodes[0])
                coordj = ops.nodeCoord(nodes[1])
                Ly = coordj[1]-coordi[1]
                Lz = coordj[2]-coordi[2]
                L = sqrt(Ly**2 + Lz**2)
                Wy = -wD*(Ly/L)*(Ly/L)
                Wx = -wD*(Ly/L)*(Lz/L)
                ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

                n = 220001+100*i+j
                nodes = ops.eleNodes(n)
                coordi = ops.nodeCoord(nodes[0])
                coordj = ops.nodeCoord(nodes[1])
                Ly = coordj[1]-coordi[1]
                Lz = coordj[2]-coordi[2]
                L = sqrt(Ly**2 + Lz**2)
                Wy = -wD*(Ly/L)*(Ly/L)
                Wx = -wD*(Ly/L)*(Lz/L)
                ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

                n = 230001+100*i+j
                nodes = ops.eleNodes(n)
                coordi = ops.nodeCoord(nodes[0])
                coordj = ops.nodeCoord(nodes[1])
                Ly = coordj[1]-coordi[1]
                Lz = coordj[2]-coordi[2]
                L = sqrt(Ly**2 + Lz**2)
                Wy = -wD*(Ly/L)*(Ly/L)
                Wx = -wD*(Ly/L)*(Lz/L)
                ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

        # Define Self Weight on Joist Girder
        for i in range(self.nspaces):
            n = 110101+i
            nodes = ops.eleNodes(n)
            coordi = ops.nodeCoord(nodes[0])
            coordj = ops.nodeCoord(nodes[1])
            Lx = coordj[0]-coordi[0]
            Lz = coordj[2]-coordi[2]
            L = sqrt(Lx**2 + Lz**2)
            Wy = -self.wdJG*(Lx/L)*(Lx/L)
            Wx = -self.wdJG*(Lx/L)*(Lz/L)
            ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

            n = 110201+i
            nodes = ops.eleNodes(n)
            coordi = ops.nodeCoord(nodes[0])
            coordj = ops.nodeCoord(nodes[1])
            Lx = coordj[0]-coordi[0]
            Lz = coordj[2]-coordi[2]
            L = sqrt(Lx**2 + Lz**2)
            Wy = -self.wdJG*(Lx/L)*(Lx/L)
            Wx = -self.wdJG*(Lx/L)*(Lz/L)
            ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

            n = 120101+i
            nodes = ops.eleNodes(n)
            coordi = ops.nodeCoord(nodes[0])
            coordj = ops.nodeCoord(nodes[1])
            Lx = coordj[0]-coordi[0]
            Lz = coordj[2]-coordi[2]
            L = sqrt(Lx**2 + Lz**2)
            Wy = -self.wdJG*(Lx/L)*(Lx/L)
            Wx = -self.wdJG*(Lx/L)*(Lz/L)
            ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

            n = 120201+i
            nodes = ops.eleNodes(n)
            coordi = ops.nodeCoord(nodes[0])
            coordj = ops.nodeCoord(nodes[1])
            Lx = coordj[0]-coordi[0]
            Lz = coordj[2]-coordi[2]
            L = sqrt(Lx**2 + Lz**2)
            Wy = -self.wdJG*(Lx/L)*(Lx/L)
            Wx = -self.wdJG*(Lx/L)*(Lz/L)
            ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

        ###########################################################
        # Run Analysis
        ###########################################################

        # Initilize data
        results = AnalysisResults()
        results.print_each_analysis_time_increment = self.print_each_analysis_time_increment
        results.analysis_type = analysis_type
        
        # Set up analysis
        ops.numberer("RCM")
        ops.constraints("Plain")
        ops.system("ProfileSPD")
        ops.test("NormUnbalance", 1.0e-6, 25, self.test_flag)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        tic = time.perf_counter()
        ops.analyze(1)
        toc = time.perf_counter()
        results.add_to_analysis_time(tic,toc)
        ops.reactions()

        #print(ops.systemSize())

        if analysis_type.lower() == 'simplesteplevel':

            # Initialize results
            results.water_volume = np.zeros((num_steps+1,1))
            results.water_level  = np.zeros((num_steps+1,1))
            results.col_react_B2 = np.zeros((num_steps+1,1))

            # Find lowest point
            zo = self.lowest_point()
        
            # Store Reuslts
            results.water_volume[0] = 0.
            results.water_level[0] = zo
            results.col_react_B2[0] = self.ColumnReaction('B2')

            # Run Ponding Analysis
            ops.timeSeries("Constant", self.ponding_load_ts_tag)
            for iStep in range(1,num_steps+1):

                # Update ponding load cells
                if self.include_ponding_effect:
                    PondingLoadManager.update()

                # Compute load vector
                izw = zo + (iStep/num_steps)*(target_zw-zo)
                (iV,idVdz) = PondingLoadManager.get_volume(izw)
                PondingLoadManager.compute_current_load_vector(izw)

                # Apply difference to model
                ops.pattern("Plain", self.ponding_load_pattern_tag_start+iStep, self.ponding_load_ts_tag)
                PondingLoadManager.apply_load_increment()
                PondingLoadManager.commit_current_load_vector()

                # Run analysis
                tic = time.perf_counter()
                ops.analyze(1)
                toc = time.perf_counter()
                results.add_to_analysis_time(tic,toc)
                ops.reactions()

                # Store Reuslts
                results.water_volume[iStep] = iV
                results.water_level[iStep]  = izw
                results.col_react_B2[iStep] = self.ColumnReaction('B2')

        elif analysis_type.lower() == 'simplestepvolume':
            raise Exception('Simple step volume analysis not yet implemented')
            
        elif analysis_type.lower() == 'iterativelevel':

            # Run Ponding Analysis
            ops.timeSeries("Constant", self.ponding_load_ts_tag)
            for iStep in range(1,self.max_iter_level+1):

                # Update ponding load cells
                if self.include_ponding_effect:
                    PondingLoadManager.update()

                # Compute load vector
                PondingLoadManager.compute_current_load_vector(target_zw)

                # Check for convergence
                if PondingLoadManager.sub_abs_diff_load_increment() < self.tol_load_level:
                    print('Converged')
                    break
                    
                # Print data on iteration
                print('Iteration: %3.i, Total Water Load: %0.3f' % (iStep,PondingLoadManager.total_current_load()))
                    
                # Apply difference to model
                ops.pattern("Plain", self.ponding_load_pattern_tag_start+iStep, self.ponding_load_ts_tag)
                PondingLoadManager.apply_load_increment()
                PondingLoadManager.commit_current_load_vector()

                # Run analysis
                tic = time.perf_counter()
                ops.analyze(1)
                toc = time.perf_counter()
                results.add_to_analysis_time(tic,toc)
                ops.reactions()

            # Store Reuslts
            (iV,idVdz) = PondingLoadManager.get_volume(target_zw)
            results.water_volume = iV
            results.water_level  = target_zw
            results.col_react_B2 = self.ColumnReaction('B2')
            
        else:
            raise Exception('Unknown analysis type: %s' % analysis_type)

        results.print_total_analysis_time()

        # Plot Ponding Load Cells
        if self.plot_load_cells:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            for i in PondingLoadManager.cells:
                coordI = PondingLoadManager.cells[i].vertexI.coord()
                coordJ = PondingLoadManager.cells[i].vertexJ.coord()
                coordK = PondingLoadManager.cells[i].vertexK.coord()
                coordL = PondingLoadManager.cells[i].vertexL.coord()
                #if i=='BC_001x_001y':
                #    print(coordI)
                #    print(coordJ)
                #    print(coordK)
                #    print(coordL)
                X = np.array([[coordJ[0], coordK[0]], [coordI[0], coordL[0]]])
                Y = np.array([[coordJ[1], coordK[1]], [coordI[1], coordL[1]]])
                Z = np.array([[coordJ[2], coordK[2]], [coordI[2], coordL[2]]])
                surf = ax.plot_surface(X, Y, Z, linewidth=1)
            plt.show()
        
        return results

    def PlotModel(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        NodeCoord = dict()
        NodeTags = ops.getNodeTags()
        for i in NodeTags:
            coord = ops.nodeCoord(i)
            NodeCoord[i] = coord
            #ax.scatter(coord[0], coord[1], coord[2], marker='o')

        EleTags = ops.getEleTags()
        for i in EleTags:
            nodes = ops.eleNodes(i)
            coordi = NodeCoord[nodes[0]]
            coordj = NodeCoord[nodes[1]]
            ax.plot([coordi[0],coordj[0]],[coordi[1],coordj[1]],[coordi[2],coordj[2]])

        plt.show()


    def ColumnReaction(self,column):
        if column == 'B2':
            #print(ops.nodeReaction(110101+self.nspaces))
            #print(ops.nodeReaction(110201))
            #print(ops.nodeReaction(210001+100*self.nspaces+self.nele_J))
            #print(ops.nodeReaction(220001+100*self.nspaces))

            R = ops.nodeReaction(110101+self.nspaces, 3) + \
                ops.nodeReaction(110201, 3) + \
                ops.nodeReaction(210001+100*self.nspaces+self.nele_J, 3) + \
                ops.nodeReaction(220001+100*self.nspaces, 3)
        elif column == 'C2':
            R = ops.nodeReaction(120101+self.nspaces, 3) + \
                ops.nodeReaction(120201, 3) + \
                ops.nodeReaction(220001+100*self.nspaces+self.nele_J, 3) + \
                ops.nodeReaction(230001+100*self.nspaces, 3)
        else:
            raise Exception('Unknown column')

        return R

