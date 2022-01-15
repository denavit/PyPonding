import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, ceil, pi
from PyPonding import FE, PondingLoadCell
from PyPonding import opensees as ops
from libdenavit import camber
from libdenavit.OpenSees import AnalysisResults

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

class IdealizedBay:
    """
    Class to perform ponding analyses on an idealized rectangular bay. 
    
    
    """
    
    def __init__(self, **attrs):
    
        # Geometry
        self.primary_member_span    = attrs['primary_member_span']
        self.secondary_member_span  = attrs['secondary_member_span']
        self.number_of_joist_spaces = attrs['number_of_joist_spaces']

        # Loads and load factors
        self.dead_load_uniform         = attrs['dead_load_uniform'] 
        self.dead_load_primary_member  = attrs['dead_load_primary_member'] # Self-weight of joist girder
        self.water_density             = attrs['water_density'] 
        self.snow_density              = attrs['snow_density'] 
        self.snow_height               = attrs['snow_height'] 
        self.alpha                     = attrs['alpha']
        self.load_factor_dead          = attrs['load_factor_dead']    # Dead
        self.load_factor_ponding       = attrs['load_factor_ponding'] # Impounded Water
        self.load_factor_snow1         = attrs['load_factor_snow1']   # Snow in Ponding Load Cell
        self.load_factor_snow2         = attrs['load_factor_snow2']   # Snow as Simple Load

        # Top of roof elevation
        self.z_TL = attrs['z_TL']
        self.z_TR = attrs['z_TR']
        self.z_BL = attrs['z_BL']
        self.z_BR = attrs['z_BR']

        # Camber
        self.primary_member_camber_T = attrs['primary_member_camber_T']
        self.primary_member_camber_B = attrs['primary_member_camber_B']
        self.secondary_member_camber = attrs['secondary_member_camber']

        # Edge conditions ('mirrored', 'rigid', or '')
        self.edge_condition_L = attrs['edge_condition_L']
        self.edge_condition_R = attrs['edge_condition_R']
        self.edge_condition_T = attrs['edge_condition_T']
        self.edge_condition_B = attrs['edge_condition_B']

        # Material and section properties
        self.E  = attrs['E']
        self.Ap = attrs['Ap']
        self.As = attrs['As']
        self.Ip = attrs['Ip']
        self.Is = attrs['Is']

        # Analysis options
        self.analsis_engine         = attrs.get('analsis_engine', 'FE')
        self.include_ponding_effect = attrs.get('include_ponding_effect', True)
        self.num_ele_secondary      = attrs.get('num_ele_secondary', 20)
        self.num_subcell_X          = attrs.get('num_subcell_X', 4) # Number of ponding sub-cells along joist direction
        self.num_subcell_Y          = attrs.get('num_subcell_Y', 4) # Number of ponding sub-cells along joist girder direction       
        self.MAX_ITER               = attrs.get('MAX_ITER', 50) # Maximum number of ponding analysis iterations
        self.tol                    = attrs.get('tol', 0.00001) # Ponding analysis tolerance
        
    def Run_Analysis(self,water_level):
           
        # Determine profile based on top of roof elevations and camber
        roof_profile = np.zeros((self.number_of_joist_spaces+1,self.num_ele_secondary+1))
        for i in range(self.number_of_joist_spaces+1):
            x = i/self.number_of_joist_spaces
            L = self.primary_member_span
            
            zB = self.z_BL + x*(self.z_BR - self.z_BL) + camber(x*L,L,self.primary_member_camber_B)
            zT = self.z_TL + x*(self.z_TR - self.z_TL) + camber(x*L,L,self.primary_member_camber_T)
            
            for j in range(self.num_ele_secondary+1):
                x = j/self.num_ele_secondary
                L = self.secondary_member_span
                c = self.secondary_member_camber
                
                if i == 0 and self.edge_condition_L == 'rigid':
                    c = 0 # Don't add camber to the left secondary member if it is rigid
                if i == self.number_of_joist_spaces and self.edge_condition_R == 'rigid':    
                    c = 0 # Don't add camber to the right secondary member if it is rigid
                    
                roof_profile[i,j] = zT + x*(zB-zT) + camber(x*L,L,self.secondary_member_camber)
                
        #print(np.transpose(roof_profile))


        # Define ponding load cells
        ponding_load_cells = [[0] * self.num_ele_secondary for i in range(self.number_of_joist_spaces)]
        for i in range(self.number_of_joist_spaces):
            for j in range(self.num_ele_secondary):
                
                iCell = PondingLoadCell.PondingLoadCell3d()
                iCell.id      = '%i,%i' % (i,j)
                
                iCell.xI      = i/self.number_of_joist_spaces*self.primary_member_span
                iCell.yI      = -(j+1)/self.num_ele_secondary*self.secondary_member_span
                iCell.zI      = roof_profile[i,j+1]
                iCell.xJ      = (i+1)/self.number_of_joist_spaces*self.primary_member_span
                iCell.yJ      = -(j+1)/self.num_ele_secondary*self.secondary_member_span
                iCell.zJ      = roof_profile[i+1,j+1]
                iCell.xK      = (i+1)/self.number_of_joist_spaces*self.primary_member_span
                iCell.yK      = -j/self.num_ele_secondary*self.secondary_member_span
                iCell.zK      = roof_profile[i+1,j]
                iCell.xL      = i/self.number_of_joist_spaces*self.primary_member_span
                iCell.yL      = -j/self.num_ele_secondary*self.secondary_member_span
                iCell.zL      = roof_profile[i,j]

                iCell.gamma   = self.alpha*self.load_factor_ponding*self.water_density
                iCell.na      = self.num_subcell_X
                iCell.nb      = self.num_subcell_Y

                iCell.gammas  = self.alpha*self.load_factor_snow1*self.snow_density
                iCell.hs      = self.snow_height

                ponding_load_cells[i][j] = iCell
        
        
        # Run Ponding Analysis
        ponding_load      = np.zeros((self.number_of_joist_spaces+1,self.num_ele_secondary+1))
        ponding_load_last = np.zeros((self.number_of_joist_spaces+1,self.num_ele_secondary+1))
        for iteration in range(self.MAX_ITER):
            
            # Reset ponding load arrays
            for i in range(self.number_of_joist_spaces+1):
                for j in range(self.num_ele_secondary+1):
                    ponding_load_last[i,j] = ponding_load[i,j]
                    ponding_load[i,j] = 0
                    
            # Compute new ponding load
            for i in range(self.number_of_joist_spaces):
                for j in range(self.num_ele_secondary):
                    if iteration == 0:
                        ponding_load_cells[i][j].dzI = 0.0
                        ponding_load_cells[i][j].dzJ = 0.0
                        ponding_load_cells[i][j].dzK = 0.0
                        ponding_load_cells[i][j].dzL = 0.0
                    else:
                        ponding_load_cells[i][j].dzI = analysis_results.bay_total_deflection[i,j+1]
                        ponding_load_cells[i][j].dzJ = analysis_results.bay_total_deflection[i+1,j+1]
                        ponding_load_cells[i][j].dzK = analysis_results.bay_total_deflection[i+1,j]
                        ponding_load_cells[i][j].dzL = analysis_results.bay_total_deflection[i,j]
                        
                    f = ponding_load_cells[i][j].get_load_vector(water_level)
                    ponding_load[i,j+1]   += f[0]
                    ponding_load[i+1,j+1] += f[1]
                    ponding_load[i+1,j]   += f[2]
                    ponding_load[i,j]     += f[3]

            # Perform analysis on secondary members
            if self.analsis_engine.lower() == 'opensees':
                analysis_results = self.run_static_analysis_OPS(ponding_load)
            elif self.analsis_engine.lower() == 'fe':
                analysis_results = self.run_static_analysis_FE(ponding_load)
            else:
                raise Exception(f'Unknown analysis engine {self.analsis_engine}')

            # Check for convergence            
            if not self.include_ponding_effect:
                break
            
            sum_of_force = -self.alpha*self.load_factor_snow1*self.primary_member_span*self.secondary_member_span*self.snow_height*self.snow_density
            sum_of_diff  = 0
            for i in range(self.number_of_joist_spaces+1):
                for j in range(self.num_ele_secondary+1):            
                    sum_of_force += abs(ponding_load[i,j])
                    sum_of_diff  += abs(ponding_load_last[i,j]-ponding_load[i,j])

            print('Iteration %02i, Total Force: %.5f' % (iteration,sum_of_force))

            if sum_of_force == 0:
                if sum_of_diff <= self.tol and iteration > 0:
                    # Converged
                    break
            else:
                if sum_of_diff/sum_of_force <= self.tol:
                    # Converged
                    break

            if iteration == self.MAX_ITER-1:
                print('The maximum number iterations has been reached without convergence')
              
        return analysis_results
        
    def run_static_analysis_FE(self,ponding_load):
    
        # Initilize Results
        deflection                      = np.zeros((self.number_of_joist_spaces+1,self.num_ele_secondary+1))
        secondary_member_deflection     = np.zeros((self.number_of_joist_spaces+1,self.num_ele_secondary+1))
        primary_member_deflection_T     = np.zeros((self.number_of_joist_spaces+1,1))
        primary_member_deflection_B     = np.zeros((self.number_of_joist_spaces+1,1))
        primary_member_ponding_load_T   = np.zeros((self.number_of_joist_spaces+1,1))
        primary_member_ponding_load_B   = np.zeros((self.number_of_joist_spaces+1,1))
        secondary_member_results        = dict()    
    
        # Build Primary Member Model
        primary_member_model = FE.Model('Primary Member')
        
        # Define Nodes
        for i in range(self.number_of_joist_spaces+1):
            n = 'n%02i' % i
            x = (i/self.number_of_joist_spaces)*self.primary_member_span
            primary_member_model.AddNode(n,(x,0.0),('UX','UY','RZ'))

            if i == 0:
                primary_member_model.Nodes[n].dofs['UX'].constrained = True
                primary_member_model.Nodes[n].dofs['UY'].constrained = True
            if i == self.number_of_joist_spaces:
                primary_member_model.Nodes[n].dofs['UY'].constrained = True
                
            # Only apply self-weight or girder, other loads will come from the reaction of the secondary members
            P_dead = self.dead_load_primary_member*(self.primary_member_span/self.number_of_joist_spaces)
            if i == 0 or i == self.number_of_joist_spaces:
                primary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead/2
            else:
                primary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead
                
        # Define Elements
        for i in range(self.number_of_joist_spaces):
            ni = 'n%02i' % i
            nj = 'n%02i' % (i+1)
            primary_member_model.AddElement('e%02i'%i,'ElasticBeam2d',(ni,nj),self.E,self.Ip,self.Ap)
        
        # Build Secondary Member Model
        secondary_member_model = FE.Model('Secondary Member')
        
        # Define Nodes
        for i in range(self.num_ele_secondary+1):
            n = 'n%02i' % i
            x = (i/self.num_ele_secondary)*self.secondary_member_span
            secondary_member_model.AddNode(n,(x,0.0),('UX','UY','RZ'))

            if i == 0:
                secondary_member_model.Nodes[n].dofs['UX'].constrained = True
                secondary_member_model.Nodes[n].dofs['UY'].constrained = True
            if i == self.num_ele_secondary:
                secondary_member_model.Nodes[n].dofs['UY'].constrained = True
              
            P_dead = self.dead_load_uniform*(self.primary_member_span/self.number_of_joist_spaces)*(self.secondary_member_span/self.num_ele_secondary)
            P_snow = self.snow_density*self.snow_height*(self.primary_member_span/self.number_of_joist_spaces)*(self.secondary_member_span/self.num_ele_secondary)
            # @todo - check this if it is mirrored
            # @todo - check if first and last nodes should have half the load
            if i == 0 or i == self.num_ele_secondary:
                secondary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead/2
                secondary_member_model.Nodes[n].dofs['UY'].loads['SNOW'] = -P_snow/2
            else:
                secondary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead
                secondary_member_model.Nodes[n].dofs['UY'].loads['SNOW'] = -P_snow
            
        # Define Elements
        for i in range(self.num_ele_secondary):
            ni = 'n%02i' % i
            nj = 'n%02i' % (i+1)
            secondary_member_model.AddElement('e%02i'%i,'ElasticBeam2d',(ni,nj),self.E,self.Is,self.As)   

        # Perform analyses on secondary members
        for i in range(self.number_of_joist_spaces+1):
            if i == 0 and self.edge_condition_L == 'rigid':
                break
            if i == self.number_of_joist_spaces and self.edge_condition_R == 'rigid':
                break
                
            # Apply water loads to the secondary members
            for j in range(self.num_ele_secondary + 1):
                if (i == 0 and self.edge_condition_L == 'mirrored') or (i == self.number_of_joist_spaces and self.edge_condition_R == 'mirrored'):
                    secondary_member_model.Nodes['n%02i' % j].dofs['UY'].loads['PONDING'] = 2*ponding_load[i,j]
                else:
                    secondary_member_model.Nodes['n%02i' % j].dofs['UY'].loads['PONDING'] = ponding_load[i,j]
        
            # Run secondary member analyses
            res = FE.LinearAnalysis(secondary_member_model)
            res.run({'DEAD':self.alpha*self.load_factor_dead,'SNOW':self.alpha*self.load_factor_snow2,'PONDING':1.0})
            
            # Get reactions (load on the primary members)
            primary_member_ponding_load_T[i] = -secondary_member_model.Nodes['n00'].dofs['UY'].react(res)
            primary_member_ponding_load_B[i] = -secondary_member_model.Nodes['n%02i' % self.num_ele_secondary].dofs['UY'].react(res)
            
            # Get member defelctions
            for j in range(self.num_ele_secondary + 1):
                secondary_member_deflection[i,j] = secondary_member_model.Nodes['n%02i' % j].dofs['UY'].disp(res)
        
            # Save results
            secondary_member_results[i] = res
            
        # Perform analysis on top primary member
        if self.edge_condition_T != 'rigid':
            # Apply water loads to the primary members
            for i in range(self.number_of_joist_spaces+1):
                if self.edge_condition_T == 'mirrored':
                    primary_member_model.Nodes['n%02i' % i].dofs['UY'].loads['PONDING'] = 2*primary_member_ponding_load_T[i]
                else:
                    primary_member_model.Nodes['n%02i' % i].dofs['UY'].loads['PONDING'] = primary_member_ponding_load_T[i]
        
            # Run analyses
            res = FE.LinearAnalysis(primary_member_model)
            res.run({'DEAD':self.alpha*self.load_factor_dead,'PONDING':1.0})
            
            # Get member defelctions
            for i in range(self.number_of_joist_spaces+1):
                primary_member_deflection_T[i] = primary_member_model.Nodes['n%02i' % i].dofs['UY'].disp(res)
                
            # Save results
            top_primary_member_results = res    
                
        # Perform analysis on bottom primary member
        if self.edge_condition_B != 'rigid':
            # Apply water loads to the primary members
            for i in range(self.number_of_joist_spaces+1):
                if self.edge_condition_B == 'mirrored':
                    primary_member_model.Nodes['n%02i' % i].dofs['UY'].loads['PONDING'] = 2*primary_member_ponding_load_B[i]
                else:
                    primary_member_model.Nodes['n%02i' % i].dofs['UY'].loads['PONDING'] = primary_member_ponding_load_B[i]
        
            # Run analyses
            res = FE.LinearAnalysis(primary_member_model)
            res.run({'DEAD':self.alpha*self.load_factor_dead,'PONDING':1.0})
            
            # Get member defelctions
            for i in range(self.number_of_joist_spaces+1):
                primary_member_deflection_B[i] = primary_member_model.Nodes['n%02i' % i].dofs['UY'].disp(res)   
        
            # Save results
            bot_primary_member_results = res  
        
        # Compute roof deflection
        for i in range(self.number_of_joist_spaces+1):
            for j in range(self.num_ele_secondary+1):            
                x = j/self.num_ele_secondary
                primary_member_deflection = primary_member_deflection_T[i] + x*(primary_member_deflection_B[i] - primary_member_deflection_T[i])
                deflection[i,j] = primary_member_deflection + secondary_member_deflection[i,j]

        # Save results (shears, moments, deflections, reactions, forces)
        results = AnalysisResults()
        results.bay_total_deflection = deflection                    
        results.top_primary_member_deflection     = primary_member_deflection_T   
        results.top_primary_member_deflection     = primary_member_deflection_B   
        results.top_primary_member_ponding_load   = primary_member_ponding_load_T 
        results.top_primary_member_ponding_load   = primary_member_ponding_load_B 

        x = np.zeros((2*self.num_ele_secondary,1))
        for i in range(2*self.num_ele_secondary):
            x[i] = self.secondary_member_span*ceil(i/2)/self.num_ele_secondary
        results.secondary_members_position = x
        
        results.secondary_members_moment = np.zeros((self.number_of_joist_spaces+1,2*self.num_ele_secondary))
        results.secondary_members_shear = np.zeros((self.number_of_joist_spaces+1,2*self.num_ele_secondary))
        results.secondary_members_top_reaction = np.zeros(self.number_of_joist_spaces+1)
        results.secondary_members_bot_reaction = np.zeros(self.number_of_joist_spaces+1)
        for i in range(self.number_of_joist_spaces+1):
            if i == 0 and self.edge_condition_L == 'rigid':
                break
            if i == self.number_of_joist_spaces and self.edge_condition_R == 'rigid':
                break            
            for j in range(self.num_ele_secondary):
                ele_force = secondary_member_model.Elements['e%02i'%j].force(secondary_member_results[i])/self.alpha
                results.secondary_members_moment[i,2*j+0] = -ele_force.item(2)
                results.secondary_members_moment[i,2*j+1] =  ele_force.item(5)
                results.secondary_members_shear[i,2*j+0] =  ele_force.item(1)
                results.secondary_members_shear[i,2*j+1] = -ele_force.item(4)
            results.secondary_members_top_reaction[i] = secondary_member_model.Nodes['n00'].dofs['UY'].react(secondary_member_results[i])/self.alpha
            results.secondary_members_bot_reaction[i] = secondary_member_model.Nodes['n%02i' % self.num_ele_secondary].dofs['UY'].react(secondary_member_results[i])/self.alpha          
               
        x = np.zeros((2*self.number_of_joist_spaces,1))
        for i in range(2*self.number_of_joist_spaces):
            x[i] = self.primary_member_span*ceil(i/2)/self.number_of_joist_spaces
        results.primary_members_position = x

        results.top_primary_member_moment = np.zeros(2*self.number_of_joist_spaces)
        results.top_primary_member_shear = np.zeros(2*self.number_of_joist_spaces)
        results.bot_primary_member_moment = np.zeros(2*self.number_of_joist_spaces)
        results.bot_primary_member_shear = np.zeros(2*self.number_of_joist_spaces)
        for i in range(self.number_of_joist_spaces):
            if self.edge_condition_T != 'rigid':
                ele_force = primary_member_model.Elements['e%02i'%i].force(top_primary_member_results)/self.alpha
                results.top_primary_member_moment[2*i+0] = -ele_force.item(2)
                results.top_primary_member_moment[2*i+1] =  ele_force.item(5)
                results.top_primary_member_shear[2*i+0] =  ele_force.item(1)
                results.top_primary_member_shear[2*i+1] = -ele_force.item(4)
        
            if self.edge_condition_B != 'rigid':
                ele_force = primary_member_model.Elements['e%02i'%i].force(bot_primary_member_results)/self.alpha
                results.bot_primary_member_moment[2*i+0] = -ele_force.item(2)
                results.bot_primary_member_moment[2*i+1] =  ele_force.item(5)
                results.bot_primary_member_shear[2*i+0] =  ele_force.item(1)
                results.bot_primary_member_shear[2*i+1] = -ele_force.item(4)
        
        return results

    def run_static_analysis_OPS(self,ponding_load):
    
        # Initilize Results
        deflection                      = np.zeros((self.number_of_joist_spaces+1,self.num_ele_secondary+1))
        secondary_member_deflection     = np.zeros((self.number_of_joist_spaces+1,self.num_ele_secondary+1))
        secondary_M                     = np.zeros((self.number_of_joist_spaces+1,2*self.num_ele_secondary))
        secondary_V                     = np.zeros((self.number_of_joist_spaces+1,2*self.num_ele_secondary))
        secondary_R_top                 = np.zeros(self.number_of_joist_spaces+1)
        secondary_R_bot                 = np.zeros(self.number_of_joist_spaces+1)  
        
        primary_member_deflection_T     = np.zeros(self.number_of_joist_spaces+1)
        primary_member_deflection_B     = np.zeros(self.number_of_joist_spaces+1)
        top_primary_M                   = np.zeros(2*self.number_of_joist_spaces)
        top_primary_V                   = np.zeros(2*self.number_of_joist_spaces)
        bot_primary_M                   = np.zeros(2*self.number_of_joist_spaces)
        bot_primary_V                   = np.zeros(2*self.number_of_joist_spaces)
  
        # Perform analyses on secondary members
        for i in range(self.number_of_joist_spaces+1):
            if i == 0 and self.edge_condition_L == 'rigid':
                break
            if i == self.number_of_joist_spaces and self.edge_condition_R == 'rigid':
                break
                
            # Create OpenSees model
            ops.wipe()
            ops.model('basic', '-ndm', 2, '-ndf', 3)

            # Define nodes
            for j in range(self.num_ele_secondary+1):
                ops.node(j,self.secondary_member_span*j/self.num_ele_secondary,0.0)

            ops.fix(0,1,1,0)
            ops.fix(self.num_ele_secondary,0,1,0)

            # Define elements
            ops.geomTransf('Linear',1)
            for j in range(0,self.num_ele_secondary):
                ops.element('elasticBeamColumn',j,j,j+1,self.As,self.E,self.Is,1)

            # Define load
            uniform_load = -self.alpha* \
                (self.load_factor_dead*self.dead_load_uniform + self.load_factor_snow2*self.snow_density*self.snow_height)* \
                (self.primary_member_span/self.number_of_joist_spaces)*(self.secondary_member_span/self.num_ele_secondary)
                
            ops.timeSeries("Constant", 1)
            ops.pattern('Plain', 1, 1)
            for j in range(self.num_ele_secondary+1):
                if i == 0:
                    # Left secondary member
                    if j == 0:
                        # Top left corner
                        if self.edge_condition_L == 'mirrored':
                            P = 0.5*uniform_load + 2*ponding_load[i,j]
                        else:
                            P = 0.25*uniform_load + ponding_load[i,j]
                    elif j < self.num_ele_secondary:
                        if self.edge_condition_L == 'mirrored':
                            P = uniform_load + 2*ponding_load[i,j]
                        else:
                            P = 0.5*uniform_load + ponding_load[i,j]
                    else:
                        # Bottom left corner
                        if self.edge_condition_L == 'mirrored':
                            P = 0.5*uniform_load + 2*ponding_load[i,j]
                        else:
                            P = 0.25*uniform_load + ponding_load[i,j]
                
                elif i < self.number_of_joist_spaces:
                    if j == 0:
                        # Top edge
                        P = 0.5*uniform_load + ponding_load[i,j]
                    elif j < self.num_ele_secondary:
                        P = uniform_load + ponding_load[i,j]
                    else:
                        # Bottom edge
                        P = 0.5*uniform_load + ponding_load[i,j]
                
                else:
                    # Right secondary member
                    if j == 0:
                        # Top right corner
                        if self.edge_condition_R == 'mirrored':
                            P = 0.5*uniform_load + 2*ponding_load[i,j]
                        else:
                            P = 0.25*uniform_load + ponding_load[i,j]
                    elif j < self.num_ele_secondary:
                        if self.edge_condition_R == 'mirrored':
                            P = uniform_load + 2*ponding_load[i,j]
                        else:
                            P = 0.5*uniform_load + ponding_load[i,j]
                    else:
                        # Bottom right corner
                        if self.edge_condition_R == 'mirrored':
                            P = 0.5*uniform_load + 2*ponding_load[i,j]
                        else:
                            P = 0.25*uniform_load + ponding_load[i,j]                    
                    
                ops.load(j,0.0,P,0.0)
                
            # Define analysis
            ops.system("UmfPack")
            ops.numberer("RCM")
            ops.constraints("Plain")
            ops.integrator("LoadControl", 0.0)
            ops.algorithm("Newton")
            ops.analysis("Static")

            # Run dead load analysis
            ops.analyze(1)
            ops.reactions()
            
            # Get reactions (load on the primary members)
            secondary_R_top[i] = ops.nodeReaction(0,2)
            secondary_R_bot[i] = ops.nodeReaction(self.num_ele_secondary,2)
            
            # Get member defelctions
            for j in range(self.num_ele_secondary + 1):
                secondary_member_deflection[i,j] = ops.nodeDisp(j, 2)
            
            # Get member forces
            for j in range(self.num_ele_secondary):
                ele_forces = ops.eleForce(j)
                secondary_V[i,2*j]   =  ele_forces[1]
                secondary_M[i,2*j]   = -ele_forces[2] 
                secondary_V[i,2*j+1] = -ele_forces[4]
                secondary_M[i,2*j+1] =  ele_forces[5] 
            
        # Perform analysis on top primary member
        if self.edge_condition_T != 'rigid':
        
            # Create OpenSees model
            ops.wipe()
            ops.model('basic', '-ndm', 2, '-ndf', 3)

            # Define nodes
            for i in range(self.number_of_joist_spaces+1):
                ops.node(i,self.primary_member_span*i/self.number_of_joist_spaces,0.0)

            ops.fix(0,1,1,0)
            ops.fix(self.number_of_joist_spaces,0,1,0)

            # Define elements
            ops.geomTransf('Linear',1)
            for i in range(0,self.number_of_joist_spaces):
                ops.element('elasticBeamColumn',i,i,i+1,self.Ap,self.E,self.Ip,1)

            # Define load
            ops.timeSeries("Constant", 1)
            ops.pattern('Plain', 1, 1)
            self_weight = -self.alpha*self.load_factor_dead*self.dead_load_primary_member* \
                (self.primary_member_span/self.number_of_joist_spaces)
            for i in range(self.number_of_joist_spaces+1):
                if self.edge_condition_T == 'mirrored':
                    P = -2*secondary_R_top[i]
                else:
                    P = secondary_R_top[i]
                if i == 0 or i == self.number_of_joist_spaces:
                    P += 0.5*self_weight
                else:
                    P += self_weight
                ops.load(i,0.0,P,0.0)
                
            # Define analysis
            ops.system("UmfPack")
            ops.numberer("RCM")
            ops.constraints("Plain")
            ops.integrator("LoadControl", 0.0)
            ops.algorithm("Newton")
            ops.analysis("Static")

            # Run dead load analysis
            ops.analyze(1)
            ops.reactions()
            
            # Get member defelctions
            for i in range(self.number_of_joist_spaces+1):
                primary_member_deflection_T[i] = ops.nodeDisp(i, 2)
                
            # Get member forces
            for i in range(self.number_of_joist_spaces):
                ele_forces = ops.eleForce(i)
                top_primary_V[2*i]   =  ele_forces[1]
                top_primary_M[2*i]   = -ele_forces[2] 
                top_primary_V[2*i+1] = -ele_forces[4]
                top_primary_M[2*i+1] =  ele_forces[5]
                
        # Perform analysis on bottom primary member
        if self.edge_condition_B != 'rigid':
        
            # Create OpenSees model
            ops.wipe()
            ops.model('basic', '-ndm', 2, '-ndf', 3)

            # Define nodes
            for i in range(self.number_of_joist_spaces+1):
                ops.node(i,self.primary_member_span*i/self.number_of_joist_spaces,0.0)

            ops.fix(0,1,1,0)
            ops.fix(self.number_of_joist_spaces,0,1,0)

            # Define elements
            ops.geomTransf('Linear',1)
            for i in range(0,self.number_of_joist_spaces):
                ops.element('elasticBeamColumn',i,i,i+1,self.Ap,self.E,self.Ip,1)

            # Define load
            ops.timeSeries("Constant", 1)
            ops.pattern('Plain', 1, 1)
            self_weight = -self.alpha*self.load_factor_dead*self.dead_load_primary_member* \
                (self.primary_member_span/self.number_of_joist_spaces)
            for i in range(self.number_of_joist_spaces+1):
                if self.edge_condition_B == 'mirrored':
                    P = -2*secondary_R_bot[i]
                else:
                    P = -secondary_R_bot[i]
                if i == 0 or i == self.number_of_joist_spaces:
                    P += 0.5*self_weight
                else:
                    P += self_weight
                ops.load(i,0.0,P,0.0)
                
            # Define analysis
            ops.system("UmfPack")
            ops.numberer("RCM")
            ops.constraints("Plain")
            ops.integrator("LoadControl", 0.0)
            ops.algorithm("Newton")
            ops.analysis("Static")

            # Run dead load analysis
            ops.analyze(1)
            ops.reactions()
            
            # Get member defelctions
            for i in range(self.number_of_joist_spaces+1):
                primary_member_deflection_B[i] = ops.nodeDisp(i, 2)
        
            # Get member forces
            for i in range(self.number_of_joist_spaces):
                ele_forces = ops.eleForce(i)
                bot_primary_V[2*i]   =  ele_forces[1]
                bot_primary_M[2*i]   = -ele_forces[2] 
                bot_primary_V[2*i+1] = -ele_forces[4]
                bot_primary_M[2*i+1] =  ele_forces[5]         
        
        # Compute roof deflection
        for i in range(self.number_of_joist_spaces+1):
            for j in range(self.num_ele_secondary+1):            
                x = j/self.num_ele_secondary
                primary_member_deflection = primary_member_deflection_T[i] + x*(primary_member_deflection_B[i] - primary_member_deflection_T[i])
                deflection[i,j] = primary_member_deflection + secondary_member_deflection[i,j]

        # Save results (shears, moments, deflections, reactions, forces)
        results = AnalysisResults()
        results.bay_total_deflection = deflection                    
        results.top_primary_member_deflection = primary_member_deflection_T   
        results.bot_primary_member_deflection = primary_member_deflection_B   

        x = np.zeros((2*self.num_ele_secondary,1))
        for i in range(2*self.num_ele_secondary):
            x[i] = self.secondary_member_span*ceil(i/2)/self.num_ele_secondary
        results.secondary_members_position = x
        
        results.secondary_members_moment = secondary_M
        results.secondary_members_shear = secondary_V
        results.secondary_members_top_reaction = secondary_R_top
        results.secondary_members_bot_reaction = secondary_R_bot
               
        x = np.zeros((2*self.number_of_joist_spaces,1))
        for i in range(2*self.number_of_joist_spaces):
            x[i] = self.primary_member_span*ceil(i/2)/self.number_of_joist_spaces
        results.primary_members_position = x

        results.top_primary_member_moment = top_primary_M
        results.top_primary_member_shear = top_primary_V
        results.bot_primary_member_moment = bot_primary_M
        results.bot_primary_member_shear = bot_primary_V

        return results

    def amplification_factor(self):
        gamma = self.alpha*self.load_factor_ponding*self.water_density
        Cp = gamma*self.secondary_member_span*self.primary_member_span**4/(pi**4*self.E*self.Ip)
        Cs = gamma*(self.primary_member_span/self.number_of_joist_spaces)*self.secondary_member_span**4/(pi**4*self.E*self.Is)
        Bp = 1/(1-1.15*Cp-Cs)
        return Bp
          
def run_example():

    # Define units
    inch = 1.0
    kip = 1.0
    lb  = kip/1000.0
    ft  = 12.0*inch
    in_per_ft = inch/ft
    ksi = kip/inch**2
    plf = lb/ft
    psf = lb/ft**2
    pcf = psf/ft
    kipft = kip*ft

    input = {
        'primary_member_span': 40*ft,
        'secondary_member_span': 45*ft,
        'number_of_joist_spaces': 5,
        'dead_load_uniform': 13.73*psf,
        'dead_load_primary_member': 68*plf,
        'water_density': 62.4*pcf,  
        'snow_density': 15.30*pcf,
        'snow_height': 0.5*0.7*(25.00*psf)/(15.30*pcf),
        'alpha': 1.0,
        'load_factor_dead':    1.0,
        'load_factor_ponding': 1.0,
        'load_factor_snow1':   1.0,
        'load_factor_snow2':   0.0,
        'z_TL':  0*inch,
        'z_TR': 20*inch,
        'z_BL':  0*inch,
        'z_BR': 20*inch,
        'secondary_member_camber': 0.625*inch,
        'primary_member_camber_T': 0.000*inch,
        'primary_member_camber_B': 0.000*inch,
        'edge_condition_L': 'mirrored',
        'edge_condition_R': 'rigid',
        'edge_condition_T': 'mirrored',
        'edge_condition_B': 'mirrored',
        'E': 29000*ksi,
        'As': 100*inch**2,
        'Ap': 100*inch**2,
        'Is': 410.7*inch**4,
        'Ip': 1184*inch**4,
        'analsis_engine': 'FE',
    }

    bay = IdealizedBay(**input)

    zw = 4.0
    results = bay.Run_Analysis(zw)
    
    #from pprint import pprint
    #pprint(vars(bay))
    #pprint(vars(results))
    
    print('\n==== Member Output ====')
    print('  Member         Max Shear      Max Moment')
    print('                  (kips)         (kip-ft) ')
    for i in range(bay.number_of_joist_spaces+1):
        Vmax = max(results.secondary_members_bot_reaction[i],results.secondary_members_top_reaction[i])
        Mmax = max(results.secondary_members_moment[i,:])/12
        print(f'Secondary {i+1:<2d}   {Vmax:>8.2f}         {Mmax:>8.2f}')
    Vmax = max(abs(results.top_primary_member_shear))
    Mmax = max(results.top_primary_member_moment)/12
    print(f'Primary Top    {Vmax:>8.2f}         {Mmax:>8.2f}')
    Vmax = max(abs(results.bot_primary_member_shear))
    Mmax = max(results.bot_primary_member_moment)/12
    print(f'Primary Bot    {Vmax:>8.2f}         {Mmax:>8.2f}')
    
    return


if __name__ == "__main__":
    run_example()

   