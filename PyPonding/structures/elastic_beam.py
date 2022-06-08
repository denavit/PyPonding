import numpy as np
from math import pi
from .. import PondingLoadManager2d
from .. import opensees as ops
from libdenavit.OpenSees import AnalysisResults
from libdenavit import camber

class ElasticBeam2d:

    def __init__(self,L,S,E,I,gamma,A=1000.0,zi=0.0,zj=0.0,c=0.0,qD=0.0):
        self.L = L          # Beam Length (span)
        self.S = S          # Tributary width
        self.E = E          # Elastic modulus
        self.I = I          # Moment of inertia
        self.A = A          # Cross-sectional area
        self.gamma = gamma  # Fluid density
        self.zi = zi        # Elevation of i-end (left)
        self.zj = zj        # Elevation of j-end (right)
        self.c = c          # Camber
        self.qD = qD        # Dead load (force per unit area)
                
        # General Analysis Options
        self.include_ponding_effect = True
        self.num_elements = 20
        self.maximum_number_of_iterations = 20
        self.maximum_number_of_iterations_level_from_volume = 20
        self.load_tolerance = 1e-3
        self.volume_tolerance = 1e-3
        
        # OpenSees Analysis Options
        self.OPS_geom_transf_type = 'Linear'
        self.OPS_element_type = 'dispBeamColumn'
        self.OPS_ponding_load_ts_tag = 100
        self.OPS_ponding_load_pattern_tag_start = 100

    def C(self):
        return (self.gamma*self.S*self.L**4)/(pi**4*self.E*self.I)

    def run_analysis_OPS(self,analysis_type,target_zw=None,target_Vw=None,num_steps=None):

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
                raise Exception('target_zw required for iterative level analysis')
             
        elif analysis_type.lower() == 'iterativevolume':
            # Lumped analysis, going directly to Vw and iterating
            if target_Vw == None:
                raise Exception('target_Vw required for iterative volume analysis')

        else:
            raise Exception('Unknown analysis type: %s' % analysis_type)

        # Create OpenSees model
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)

        # Define nodes
        for i in range(self.num_elements+1):
            x_over_L = float(i)/(self.num_elements)
            ix = self.L*x_over_L
            iy = self.zi+(self.zj-self.zi)*x_over_L + camber(ix,self.L,self.c)
            ops.node(i,ix,iy)

        ops.fix(0,1,1,0)
        ops.fix(self.num_elements,0,1,0)

        # Define elements
        ops.geomTransf(self.OPS_geom_transf_type,1)
        ops.section('Elastic',1,self.E,self.A,self.I)
        ops.beamIntegration('Lobatto', 1, 1, 3)
        for i in range(0,self.num_elements):
            ops.element(self.OPS_element_type,i,i,i+1,1,1)

        # Define dead load
        ops.timeSeries("Constant", 1)
        ops.pattern('Plain', -1, 1)
        for i in range(0,self.num_elements):
            ops.eleLoad('-ele',i,'-type','beamUniform',-self.qD*self.S)
            
        # Define analysis
        ops.system("UmfPack")
        ops.numberer("RCM")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.0)
        ops.algorithm("Newton")
        ops.analysis("Static")

        # Run dead load analysis
        ops.analyze(1)
        
        # Determine low point along beam
        zo = float('inf')
        for i in range(self.num_elements+1):
            if self.include_ponding_effect:
                iz = ops.nodeCoord(i, 2) + ops.nodeDisp(i, 2)
            else:
                iz = ops.nodeCoord(i, 2)
            if iz < zo:
                zo = iz
        
        
        # Initilize data
        results = AnalysisResults()
        results.analysis_type = analysis_type
        
        # Define ponding load cells
        PondingLoadManager = PondingLoadManager2d()
        for i in range(0,self.num_elements):
            PondingLoadManager.add_cell(i,i,i+1,self.gamma,self.S)
        
        # Run Ponding Analysis
        ops.timeSeries("Constant", self.OPS_ponding_load_ts_tag)
        
        if analysis_type.lower() == 'simplesteplevel':

            # Initialize results
            results.target_zw = target_zw
            results.water_volume = np.zeros((num_steps+1,1))
            results.water_level = np.zeros((num_steps+1,1))

            # Store dead load results
            results.water_volume[0] = 0.
            results.water_level[0] = zo
        
            for iStep in range(1,num_steps+1):

                # Update ponding load cells
                if self.include_ponding_effect:
                    PondingLoadManager.update()

                # Compute load vector
                izw = zo + (iStep/num_steps)*(target_zw-zo)
                (iV,idVdz) = PondingLoadManager.get_volume(izw) # @todo - should the volume be computed before or after the analysis step?
                PondingLoadManager.compute_current_load_vector(izw)

                # Apply difference to model
                ops.pattern("Plain", self.ponding_load_pattern_tag_start+iStep, self.OPS_ponding_load_ts_tag)
                PondingLoadManager.apply_load_increment()
                PondingLoadManager.commit_current_load_vector()

                # Run analysis
                ops.analyze(1)
                ops.reactions()

                # Store Reuslts
                results.water_volume[iStep] = iV
                results.water_level[iStep] = izw
                
        elif analysis_type.lower() == 'simplestepvolume':
        
            # Initialize results
            results.target_Vw = target_Vw
            results.water_volume = np.zeros((num_steps+1,1))
            results.water_level = np.zeros((num_steps+1,1))
            
            # Store dead load results
            results.water_volume[0] = 0.
            results.water_level[0] = zo
        
            for iStep in range(1,num_steps+1):

                # Update ponding load cells
                if self.include_ponding_effect:
                    PondingLoadManager.update()

                # Estimate water height
                step_Vw = (iStep+1)/num_steps*target_Vw
                if iStep == 1:
                    izw = zo+0.1 # Initial guess
                for i in range(self.maximum_number_of_iterations_level_from_volume):
                    (iV,idVdz) = PondingLoadManager.get_volume(izw)
                    izw = izw - (iV-step_Vw)/idVdz
                    if abs(iV-step_Vw) <= self.volume_tolerance:
                        break

                # Compute load vector
                PondingLoadManager.compute_current_load_vector(izw)

                # Apply difference to model
                ops.pattern("Plain", self.OPS_ponding_load_pattern_tag_start+iStep, self.OPS_ponding_load_ts_tag)
                PondingLoadManager.apply_load_increment()
                PondingLoadManager.commit_current_load_vector()

                # Run analysis
                ops.analyze(1)
                ops.reactions()

                # Store Reuslts
                results.water_volume[iStep] = step_Vw
                results.water_level[iStep] = izw 
            
        elif analysis_type.lower() == 'iterativelevel':

            iStep = 1
            while True:
                if iStep > self.maximum_number_of_iterations:
                    raise Exception(f'Analysis did not converge in {self.maximum_number_of_iterations} iterations')

                # Update ponding load cells
                if self.include_ponding_effect:
                    PondingLoadManager.update()

                # Compute load vector
                PondingLoadManager.compute_current_load_vector(target_zw)

                # Check for convergence
                if PondingLoadManager.sub_abs_diff_load_increment() < self.load_tolerance:
                    print('Converged')
                    break
                    
                # Print data on iteration
                print('Iteration: %3.i, Total Water Load: %0.3f' % (iStep,PondingLoadManager.total_current_load()))
                    
                # Apply difference to model
                ops.pattern("Plain", self.OPS_ponding_load_pattern_tag_start+iStep, self.OPS_ponding_load_ts_tag)
                PondingLoadManager.apply_load_increment()
                PondingLoadManager.commit_current_load_vector()

                # Run analysis
                ops.analyze(1)
                ops.reactions()
            
                # Increment step counter
                iStep += 1
                
            # Store Results
            (V,dVdz) = PondingLoadManager.get_volume(target_zw)
            results.water_volume = V
            results.water_level  = target_zw  
            self.add_results_along_length_OPS(results)
            
        elif analysis_type.lower() == 'iterativevolume':

            iStep = 1
            while True:
                if iStep > self.maximum_number_of_iterations:
                    raise Exception(f'Analysis did not converge in {self.maximum_number_of_iterations} iterations')

                # Update ponding load cells
                if self.include_ponding_effect:
                    PondingLoadManager.update()

                # Estimate water height
                if iStep == 1:
                    izw = zo+0.1 # Initial guess
                for i in range(self.maximum_number_of_iterations_level_from_volume):
                    (iV,idVdz) = PondingLoadManager.get_volume(izw)
                    izw = izw - (iV-target_Vw)/idVdz
                    if abs(iV-target_Vw) <= self.volume_tolerance:
                        break

                # Compute load vector
                PondingLoadManager.compute_current_load_vector(izw)

                # Check for convergence
                if PondingLoadManager.sub_abs_diff_load_increment() < self.load_tolerance:
                    print('Converged')
                    break
                    
                # Print data on iteration
                print('Iteration: %3.i, Total Water Load: %0.3f' % (iStep,PondingLoadManager.total_current_load()))
                    
                # Apply difference to model
                ops.pattern("Plain", self.OPS_ponding_load_pattern_tag_start+iStep, self.OPS_ponding_load_ts_tag)
                PondingLoadManager.apply_load_increment()
                PondingLoadManager.commit_current_load_vector()

                # Run analysis
                ops.analyze(1)
                ops.reactions()
            
                # Increment step counter
                iStep += 1
                
            # Store Results
            (V,dVdz) = PondingLoadManager.get_volume(izw)
            results.water_volume = V
            results.water_level  = izw  
            self.add_results_along_length_OPS(results)
            
        else:
            raise Exception('Unknown analysis type: %s' % analysis_type)
        
        return results
        
    def add_results_along_length_OPS(self,results):
        results.position_along_length       = np.zeros((2*self.num_elements,1))
        results.deflection_along_length     = np.zeros((2*self.num_elements,1))
        results.axial_load_along_length     = np.zeros((2*self.num_elements,1))
        results.shear_along_length          = np.zeros((2*self.num_elements,1))
        results.bending_moment_along_length = np.zeros((2*self.num_elements,1))
        
        for i in range(self.num_elements+1):
            if i == 0:
                results.position_along_length[2*i] = (float(i)/self.num_elements)*self.L
                results.deflection_along_length[2*i] = ops.nodeDisp(i,2)
            elif i < self.num_elements:
                results.position_along_length[2*i-1] = (float(i)/self.num_elements)*self.L
                results.deflection_along_length[2*i-1] = ops.nodeDisp(i,2)
                results.position_along_length[2*i] = results.position_along_length[2*i-1]
                results.deflection_along_length[2*i] = results.deflection_along_length[2*i-1]
            else:
                results.position_along_length[2*i-1] = (float(i)/self.num_elements)*self.L
                results.deflection_along_length[2*i-1] = ops.nodeDisp(i,2)
        
        for i in range(self.num_elements):
            ele_forces = ops.eleForce(i)
            results.axial_load_along_length[2*i]     = ele_forces[0]
            results.shear_along_length[2*i]          = ele_forces[1]
            results.bending_moment_along_length[2*i] = -ele_forces[2]   
            results.axial_load_along_length[2*i+1]     = -ele_forces[3]
            results.shear_along_length[2*i+1]          = -ele_forces[4]
            results.bending_moment_along_length[2*i+1] = ele_forces[5]   
        
        return