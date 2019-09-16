import numpy as np
from math import sqrt
from PyPonding import FE
from PyPonding.structures import basic_structure

class steel_beam(basic_structure.basic_structure):

    # Geometry
    L   = 40*12
    tw  = 5*12
    zi  = 0
    zj  = 0
    c   = 0
    
    # Material and section properties
    E   = 29000
    A   = 100
    I   = 182

    # Strength properties 
    Z   = float('nan')
    Fy  = float('nan')
    
    # Analysis options
    nele = 20
    dof_types = ('UX','UY','RZ')
    
    def __init__(self):
        pass
        
    def BuildModel(self,stiffness_reduction = 1.0):
        model = FE.Model('2d Steel Beam')
        
        # Define Nodes
        for i in range(self.nele+1):
            xi = (i/self.nele)
            
            # Compute Coordinates
            x = xi*self.L
            y = self.zi + xi*(self.zj-self.zi)
            if self.c != 0:
                r = (self.c**2 + (self.L/2)**2)/(2*self.c)
                y += sqrt(r**2 - (xi*self.L - self.L/2)**2) + (self.c-r)
                
            # Define Nodes
            n = 'n%02i' % i
            model.AddNode(n,(x,y),self.dof_types)

            # Assign Boundary Conditions 
            if i == 0:
                model.Nodes[n].dofs['UX'].constrained = True
                model.Nodes[n].dofs['UY'].constrained = True
            if i == self.nele:
                model.Nodes[n].dofs['UY'].constrained = True
                
            # Assign Dead and Snow Load
            if i == 0 or i == self.nele:
                model.Nodes[n].dofs['UY'].loads['DEAD'] = -self.wd*self.tw*(self.L/self.nele)/2
                model.Nodes[n].dofs['UY'].loads['SNOW'] = -self.gammas*self.hs*self.tw*(self.L/self.nele)/2
            else:
                model.Nodes[n].dofs['UY'].loads['DEAD'] = -self.wd*self.tw*(self.L/self.nele)
                model.Nodes[n].dofs['UY'].loads['SNOW'] = -self.gammas*self.hs*self.tw*(self.L/self.nele)
                
        # Define Elements and Ponding Load Cells
        for i in range(self.nele):
            ni = 'n%02i' % i
            nj = 'n%02i' % (i+1)
            model.AddElement('e%02i'%i,'ElasticBeam2d',(ni,nj),stiffness_reduction*self.E,self.I,self.A)
            model.AddPondingLoadCell('p%02i'%i,'2d',(ni,nj),self.alpha*self.LF_P*self.gamma,self.tw)
            model.PondingLoadCells['p%02i'%i].gammas = self.alpha*self.LF_S1*self.gammas
            model.PondingLoadCells['p%02i'%i].hs = self.hs
        
        self.model = model

    def Reactions(self,results):
        Ri = self.model.Nodes['n00'].dofs['UY'].react(results)/self.alpha
        Rj = self.model.Nodes['n%02i'%self.nele].dofs['UY'].react(results)/self.alpha
        return (Ri,Rj)
        
    def Moment_and_Shear(self,results):
        x = np.empty([2*self.nele,1])
        M = np.empty([2*self.nele,1])
        V = np.empty([2*self.nele,1])
        for i in range(self.nele):
            ele_force = self.model.Elements['e%02i'%i].force(results)/self.alpha
            x[2*i+0] = i*self.L/self.nele
            x[2*i+1] = (i+1)*self.L/self.nele
            M[2*i+0] = -ele_force.item(2)
            M[2*i+1] =  ele_force.item(5)
            V[2*i+0] =  ele_force.item(1)
            V[2*i+1] = -ele_force.item(4)            
        return (x,M,V)
       
    def Maximum_Shear(self,results):
        (Ri,Rj) = self.Reactions(results)
        Vmax = max(Ri,Rj)
        return Vmax   
    
    def Maximum_Moment(self,results):
        (x,M,V) = self.Moment_and_Shear(results)
        Mmax = np.amax(np.absolute(M))
        return Mmax
        
    def Strength_Ratio(self,results):
        (x,M,V) = self.Moment_and_Shear(results)
        
        # Moment 
        Mmax = np.amax(np.absolute(M))
        ind  = np.argmax(np.absolute(M))
        
        SR      = Mmax/(self.Z*self.Fy)
        SR_note = 'Moment at x/L = %0.3f' % (x[ind]/self.L) 
        
        # Shear
        # @todo - check shear also
        
        return (SR,SR_note)
       
    def lowest_point(self):
        return min([self.zi, self.zj])
    
    def determine_stiffness_reduction(self,zw,Mmax=0.0,tol=0.00001):
        if Mmax == 0:
            Mmax = self.Z*self.Fy
        
        # Try without any stiffness reduction
        tau = 1.0
        self.BuildModel(tau)
        PA = FE.PondingAnalysis(self.model,'Constant_Level')
        PA.max_iterations_z = 60
        res = PA.run({'DEAD':self.alpha*self.LF_D,'SNOW':self.alpha*self.LF_S2},zw)
        if res != 0:
            print('Not converged')
            return float('nan')
        MR = self.Maximum_Moment(PA)/Mmax
        print('tau = %7.4f, Moment Ratio = %10.7f' % (tau,MR))
        
        # First Level of Iteration
        if abs(MR-1.0) < tol:
            return 1.0
        elif MR < 1.0:
            while (tau > 0.00):
                if MR < 1.0:
                    above_tau = tau
                    above_MR  = MR
                    if tau > 0.05:
                        tau = tau - 0.05
                    else:
                        tau = tau - 0.01
                else:
                    below_tau = tau
                    below_MR  = MR
                    break
            
                self.BuildModel(tau)
                PA = FE.PondingAnalysis(self.model,'Constant_Level')
                PA.max_iterations_z = 60
                res = PA.run({'DEAD':self.alpha*self.LF_D,'SNOW':self.alpha*self.LF_S2},zw)
                if res != 0:
                    print('Not converged')
                    return float('nan')
                MR = self.Maximum_Moment(PA)/Mmax
                print('tau = %7.4f, Moment Ratio = %10.7f' % (tau,MR))
                
            else:
                print('Minimum stiffness reduction reached')
                return float('nan')
        else:
            print('Not implemented in for initial case of MmaxA > Mmax')
            return float('nan')
        
        
        # Second Level of Iteration
        MR = 0
        while (abs(MR-1) > tol):
            tau = below_tau + (above_tau-below_tau)*(1-below_MR)/(above_MR-below_MR)
            
            self.BuildModel(tau)
            PA = FE.PondingAnalysis(self.model,'Constant_Level')
            PA.max_iterations_z = 60
            res = PA.run({'DEAD':self.alpha*self.LF_D,'SNOW':self.alpha*self.LF_S2},zw)
            
            if res != 0:
                print('Not converged')
                return float('nan')
            MR = self.Maximum_Moment(PA)/Mmax
            print('tau = %7.4f, Moment Ratio = %10.7f' % (tau,MR))
            
            if MR < 1:
                above_tau = tau
                above_MR  = MR
            else:
                below_tau = tau
                below_MR  = MR            

        return tau