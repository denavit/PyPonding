import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, ceil, pi
from PyPonding import FE, PondingLoadCell

#np.set_printoptions(precision=3)
np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

class bay:

    # Geometry
    primary_member_span         = 40*12             # in
    secondary_member_span       = 40*12             # in
    num_spaces                  = 8
    
    # Loads and load factors
    alpha   = 1
    load_factor_dead            = 1.2               # Dead
    load_factor_ponding         = 1.2               # Impounded Water
    load_factor_snow1           = 1.2               # Snow in Ponding Load Cell
    load_factor_snow2           = 0.0               # Snow as Simple Load

    include_ponding_effect      = True
    
    dead_load_uniform           = 10/1000/12**2
    dead_load_primary_member    = 0                 # Self-weight of joist girder

    water_density               = 62.4/1000/12**3   
    
    snow_density                = 20/1000/12**3
    snow_height                 = 12
    
    # Top of roof elevation
    z_TL                        = 0
    z_TR                        = 0
    z_BL                        = 0
    z_BR                        = 0
    
    # Camber
    primary_member_camber_T     = 0.0
    primary_member_camber_B     = 0.0
    secondary_member_camber     = 0.0
    
    # Edge conditions ('mirrored', 'rigid', or '')
    edge_condition_L            = 'mirrored'
    edge_condition_R            = 'mirrored'
    edge_condition_T            = 'mirrored'
    edge_condition_B            = 'mirrored'
    
    # Material and section properties
    E                           = 29000
    Ap                          = 100
    As                          = 100
    Ip                          = 182
    Is                          = 182

    # Analysis options
    num_ele_secondary           = 20
    num_subcell_X               = 4                 # Number of ponding sub-cells along joist direction
    num_subcell_Y               = 4                 # Number of ponding sub-cells along joist girder direction       
    MAX_ITER                    = 50                # Maximum number of ponding analysis iterations
    tol                         = 0.00001           # Ponding analysis tolerance

    def __init__(self):
        pass
        
    def Run_Analysis(self,water_level):
    
        # Build Primary Member Model
        primary_member_model = FE.Model('Primary Member')
        
        # Define Nodes
        for i in range(self.num_spaces+1):
            n = 'n%02i' % i
            x = (i/self.num_spaces)*self.primary_member_span
            primary_member_model.AddNode(n,(x,0.0),('UX','UY','RZ'))

            if i == 0:
                primary_member_model.Nodes[n].dofs['UX'].constrained = True
                primary_member_model.Nodes[n].dofs['UY'].constrained = True
            if i == self.num_spaces:
                primary_member_model.Nodes[n].dofs['UY'].constrained = True
                
            # Only apply self-weight or girder, other loads will come from the reaction of the secondary members
            P_dead = self.dead_load_primary_member*(self.primary_member_span/self.num_spaces)
            if i == 0 or i == self.num_spaces:
                primary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead/2
            else:
                primary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead
                
        # Define Elements
        for i in range(self.num_spaces):
            ni = 'n%02i' % i
            nj = 'n%02i' % (i+1)
            primary_member_model.AddElement('e%02i'%i,'ElasticBeam2d',(ni,nj),self.E,self.Ip,self.Ap)
        
        self.primary_member_model = primary_member_model

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
              
            P_dead = self.dead_load_uniform*(self.primary_member_span/self.num_spaces)*(self.secondary_member_span/self.num_ele_secondary)
            P_snow = self.snow_density*self.snow_height*(self.primary_member_span/self.num_spaces)*(self.secondary_member_span/self.num_ele_secondary)
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
        
        self.secondary_member_model = secondary_member_model

        # Determine profile based on top of roof elevations and camber
        roof_profile = np.zeros((self.num_spaces+1,self.num_ele_secondary+1))
        for i in range(self.num_spaces+1):
            x = i/self.num_spaces
                        
            zB = self.z_BL + x*(self.z_BR - self.z_BL)
            if self.primary_member_camber_B != 0:
                r = self.primary_member_span**2/(8*self.primary_member_camber_B) + self.primary_member_camber_B/2
                zB += sqrt(r**2 - ((x-0.5)*self.primary_member_span)**2) - (r - self.primary_member_camber_B)
                
            zT = self.z_TL + x*(self.z_TR - self.z_TL)
            if self.primary_member_camber_T != 0:
                r = self.primary_member_span**2/(8*self.primary_member_camber_T) + self.primary_member_camber_T/2
                zT += sqrt(r**2 - ((x-0.5)*self.primary_member_span)**2) - (r - self.primary_member_camber_T)
        
            for j in range(self.num_ele_secondary+1):
                x = j/self.num_ele_secondary
                
                z = zT + x*(zB-zT)
                if self.secondary_member_camber != 0: # @todo - don't add camber if joist is rigid
                    r = self.secondary_member_span**2/(8*self.secondary_member_camber) + self.secondary_member_camber/2
                    z += sqrt(r**2 - ((x-0.5)*self.secondary_member_span)**2) - (r - self.secondary_member_camber)    
                
                roof_profile[i,j] = z
        
        
        # print(np.transpose(roof_profile))
        
        
        # Define ponding load cells
        ponding_load_cells = [[0] * self.num_ele_secondary for i in range(self.num_spaces)]
        for i in range(self.num_spaces):
            for j in range(self.num_ele_secondary):
                
                iCell = PondingLoadCell.PondingLoadCell3d()
                iCell.id      = '%i,%i' % (i,j)
                
                iCell.xI      = i/self.num_spaces*self.primary_member_span
                iCell.yI      = -(j+1)/self.num_ele_secondary*self.secondary_member_span
                iCell.zI      = roof_profile[i,j+1]
                iCell.xJ      = (i+1)/self.num_spaces*self.primary_member_span
                iCell.yJ      = -(j+1)/self.num_ele_secondary*self.secondary_member_span
                iCell.zJ      = roof_profile[i+1,j+1]
                iCell.xK      = (i+1)/self.num_spaces*self.primary_member_span
                iCell.yK      = -j/self.num_ele_secondary*self.secondary_member_span
                iCell.zK      = roof_profile[i+1,j]
                iCell.xL      = i/self.num_spaces*self.primary_member_span
                iCell.yL      = -j/self.num_ele_secondary*self.secondary_member_span
                iCell.zL      = roof_profile[i,j]

                iCell.gamma   = self.alpha*self.load_factor_ponding*self.water_density
                iCell.na      = self.num_subcell_X
                iCell.nb      = self.num_subcell_Y

                iCell.gammas  = self.alpha*self.load_factor_snow1*self.snow_density
                iCell.hs      = self.snow_height

                ponding_load_cells[i][j] = iCell
        
        
        
        deflection                      = np.zeros((self.num_spaces+1,self.num_ele_secondary+1))
        ponding_load                    = np.zeros((self.num_spaces+1,self.num_ele_secondary+1))
        ponding_load_last               = np.zeros((self.num_spaces+1,self.num_ele_secondary+1))
        secondary_member_deflection     = np.zeros((self.num_spaces+1,self.num_ele_secondary+1))
        primary_member_deflection_T     = np.zeros((self.num_spaces+1,1))
        primary_member_deflection_B     = np.zeros((self.num_spaces+1,1))
        primary_member_ponding_load_T   = np.zeros((self.num_spaces+1,1))
        primary_member_ponding_load_B   = np.zeros((self.num_spaces+1,1))
        secondary_member_results        = dict()
        
        for iteration in range(self.MAX_ITER):
            
            for i in range(self.num_spaces+1):
                for j in range(self.num_ele_secondary+1):
                    ponding_load_last[i,j] = ponding_load[i,j]
                    ponding_load[i,j] = 0
                    
            # Get water loads
            for i in range(self.num_spaces):
                for j in range(self.num_ele_secondary):
                    ponding_load_cells[i][j].dzI = deflection[i,j+1]
                    ponding_load_cells[i][j].dzJ = deflection[i+1,j+1]
                    ponding_load_cells[i][j].dzK = deflection[i+1,j]
                    ponding_load_cells[i][j].dzL = deflection[i,j]
                    f = ponding_load_cells[i][j].get_load_vector(water_level)
                    ponding_load[i,j+1]   += f[0]
                    ponding_load[i+1,j+1] += f[1]
                    ponding_load[i+1,j]   += f[2]
                    ponding_load[i,j]     += f[3]

            # Perform analysis on secondary members
            for i in range(self.num_spaces+1):
                if i == 0 and self.edge_condition_L == 'rigid':
                    break
                if i == self.num_spaces and self.edge_condition_R == 'rigid':
                    break
                    
                # Apply water loads to the secondary members
                for j in range(self.num_ele_secondary + 1):
                    if (i == 0 and self.edge_condition_L == 'mirrored') or (i == self.num_spaces and self.edge_condition_R == 'mirrored'):
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
                for i in range(self.num_spaces+1):
                    if self.edge_condition_T == 'mirrored':
                        primary_member_model.Nodes['n%02i' % i].dofs['UY'].loads['PONDING'] = 2*primary_member_ponding_load_T[i]
                    else:
                        primary_member_model.Nodes['n%02i' % i].dofs['UY'].loads['PONDING'] = primary_member_ponding_load_T[i]
            
                # Run analyses
                res = FE.LinearAnalysis(primary_member_model)
                res.run({'DEAD':self.alpha*self.load_factor_dead,'PONDING':1.0})
                
                # Get member defelctions
                for i in range(self.num_spaces+1):
                    primary_member_deflection_T[i] = primary_member_model.Nodes['n%02i' % i].dofs['UY'].disp(res)
                    
                # Save results
                top_primary_member_results = res    
                    
            # Perform analysis on bottom primary member
            if self.edge_condition_B != 'rigid':
                # Apply water loads to the primary members
                for i in range(self.num_spaces+1):
                    if self.edge_condition_B == 'mirrored':
                        primary_member_model.Nodes['n%02i' % i].dofs['UY'].loads['PONDING'] = 2*primary_member_ponding_load_B[i]
                    else:
                        primary_member_model.Nodes['n%02i' % i].dofs['UY'].loads['PONDING'] = primary_member_ponding_load_B[i]
            
                # Run analyses
                res = FE.LinearAnalysis(primary_member_model)
                res.run({'DEAD':self.alpha*self.load_factor_dead,'PONDING':1.0})
                
                # Get member defelctions
                for i in range(self.num_spaces+1):
                    primary_member_deflection_B[i] = primary_member_model.Nodes['n%02i' % i].dofs['UY'].disp(res)   
            
                # Save results
                bot_primary_member_results = res  
            
            # Compute roof deflection
            for i in range(self.num_spaces+1):
                for j in range(self.num_ele_secondary+1):            
                    x = j/self.num_ele_secondary
                    primary_member_deflection = primary_member_deflection_T[i] + x*(primary_member_deflection_B[i] - primary_member_deflection_T[i])
                    deflection[i,j] = primary_member_deflection + secondary_member_deflection[i,j]

            # Check for convergence            
            if not self.include_ponding_effect:
                break
            
            sum_of_force = -self.alpha*self.load_factor_snow1*self.primary_member_span*self.secondary_member_span*self.snow_height*self.snow_density
            sum_of_diff  = 0
            for i in range(self.num_spaces+1):
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
        
        # Save results (shears, moments, deflections, reactions, forces)         
        x = np.zeros((2*self.num_ele_secondary,1))
        for i in range(2*self.num_ele_secondary):
            x[i] = self.secondary_member_span*ceil(i/2)/self.num_ele_secondary
        self.results_secondary_x = x
        
        self.results_secondary_M = np.zeros((self.num_spaces+1,2*self.num_ele_secondary))
        self.results_secondary_V = np.zeros((self.num_spaces+1,2*self.num_ele_secondary))
        self.results_secondary_R_top = np.zeros((self.num_spaces+1,1))
        self.results_secondary_R_bot = np.zeros((self.num_spaces+1,1))
        for i in range(self.num_spaces+1):
            for j in range(self.num_ele_secondary):
                ele_force = secondary_member_model.Elements['e%02i'%j].force(secondary_member_results[i])/self.alpha
                self.results_secondary_M[i,2*j+0] = -ele_force.item(2)
                self.results_secondary_M[i,2*j+1] =  ele_force.item(5)
                self.results_secondary_V[i,2*j+0] =  ele_force.item(1)
                self.results_secondary_V[i,2*j+1] = -ele_force.item(4)
            
            self.results_secondary_R_top[i] = secondary_member_model.Nodes['n00'].dofs['UY'].react(secondary_member_results[i])/self.alpha
            self.results_secondary_R_bot[i] = secondary_member_model.Nodes['n%02i' % self.num_ele_secondary].dofs['UY'].react(secondary_member_results[i])/self.alpha          
               
        x = np.zeros((2*self.num_spaces,1))
        for i in range(2*self.num_spaces):
            x[i] = self.primary_member_span*ceil(i/2)/self.num_spaces
        self.results_primary_x = x

        self.results_top_primary_M = np.zeros((2*self.num_spaces,1))
        self.results_top_primary_V = np.zeros((2*self.num_spaces,1))
        self.results_bot_primary_M = np.zeros((2*self.num_spaces,1))
        self.results_bot_primary_V = np.zeros((2*self.num_spaces,1))
        for i in range(self.num_spaces):
            if self.edge_condition_T != 'rigid':
                ele_force = primary_member_model.Elements['e%02i'%i].force(top_primary_member_results)/self.alpha
                self.results_top_primary_M[2*i+0] = -ele_force.item(2)
                self.results_top_primary_M[2*i+1] =  ele_force.item(5)
                self.results_top_primary_V[2*i+0] =  ele_force.item(1)
                self.results_top_primary_V[2*i+1] = -ele_force.item(4)
        
            if self.edge_condition_B != 'rigid':
                ele_force = primary_member_model.Elements['e%02i'%i].force(bot_primary_member_results)/self.alpha
                self.results_bot_primary_M[2*i+0] = -ele_force.item(2)
                self.results_bot_primary_M[2*i+1] =  ele_force.item(5)
                self.results_bot_primary_V[2*i+0] =  ele_force.item(1)
                self.results_bot_primary_V[2*i+1] = -ele_force.item(4)
        
    def Amplification_Factor(self):
        gamma = self.alpha*self.load_factor_ponding*self.water_density
        Cp = gamma*self.secondary_member_span*self.primary_member_span**4/(pi**4*self.E*self.Ip)
        Cs = gamma*(self.primary_member_span/self.num_spaces)*self.secondary_member_span**4/(pi**4*self.E*self.Is)
        AF = 1/(1-1.15*Cp-Cs)
        return AF
          
        