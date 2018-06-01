import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from PyPonding import FE
from PyPonding.structures import basic_structure

class joist(basic_structure.basic_structure):

    # Geometry
    L   = 40*12
    tw  = 5*12
    zi  = 0
    zj  = 0
    c   = 5/8
    
    # Material and section properties
    E   = 29000
    A   = 100
    I   = 182

    # Strength properties 
    w   = 277/1000/12
    shear_reversal_strength = 0.125
    
    # Analysis options
    nele = 20
    dof_types = ('UX','UY','RZ')
    
    def __init__(self):
        pass
        
    def BuildModel(self):
        model = FE.Model('Beam Example')
        
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
            model.AddElement('e%02i'%i,'ElasticBeam2d',(ni,nj),self.E,self.I,self.A)
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
        
    def Strength_Ratio(self,results):
        SR = 0
        SR_note = ''
        
        # Reaction
        (Ri,Rj) = self.Reactions(results)
        iR = self.w*self.L/2
        if Ri/iR > SR:
            SR = Ri/iR
            SR_note = 'Reaction at I end'
   
        if Rj/iR > SR:
            SR = Rj/iR
            SR_note = 'Reaction at J end'
        
        for i in range(self.nele):
            ele_force = self.model.Elements['e%02i'%i].force(results)/self.alpha
                
            # Moment at I-end
            M = -ele_force.item(2)
            x = i*self.L/self.nele
            iM = self.w*x*(self.L-x)/2
            if iM > 0:
                if M/iM > SR:
                    SR = M/iM
                    SR_note = 'Moment at x/L = %0.3f' % (x/self.L)            
            
            # Moment at J-end
            M = ele_force.item(5)
            x = (i+1)*self.L/self.nele
            iM = self.w*x*(self.L-x)/2
            if iM > 0:
                if M/iM > SR:
                    SR = M/iM
                    SR_note = 'Moment at x/L = %0.3f' % (x/self.L)            
            
            # Shear at middle of element
            V = ele_force.item(1)        
            x = (i+0.5)*self.L/self.nele
            if x < 0.5*self.L:
                if V >= 0:
                    iV = max(self.w*(self.L/2-x),0.25*iR)
                    if V/iV > SR:
                        SR = V/iV
                        SR_note = 'Shear at x/L = %0.3f' % (x/self.L)
                else:
                    iV = -self.shear_reversal_strength*iR
                    if V/iV > SR:
                        SR = V/iV
                        SR_note = 'Shear Reversal at x/L = %0.3f' % (x/self.L)                    
            else:
                if V >= 0:
                    iV = self.shear_reversal_strength*iR
                    if V/iV > SR:
                        SR = V/iV
                        SR_note = 'Shear Reversal at x/L = %0.3f' % (x/self.L)                                        
                else:
                    iV = min(self.w*(self.L/2-x),-0.25*iR)             
                    if V/iV > SR:
                        SR = V/iV
                        SR_note = 'Shear at x/L = %0.3f' % (x/self.L)
        
        return (SR,SR_note)

    def plot_results(self,PA):
        # Make a plot
        fig, ax = plt.subplots()
        
        # Undeformed Shape
        for iNode in self.model.Nodes:
            ax.scatter(self.model.Nodes[iNode].coords[0], self.model.Nodes[iNode].coords[1], color='0.75')
        for iElement in self.model.Elements:
            x = (self.model.Elements[iElement].nodeI.coords[0],self.model.Elements[iElement].nodeJ.coords[0])
            y = (self.model.Elements[iElement].nodeI.coords[1],self.model.Elements[iElement].nodeJ.coords[1]) 
            ax.plot(x, y, color='0.75')
            
        # Deformed Shape
        for iNode in self.model.Nodes:
            ax.scatter(self.model.Nodes[iNode].coords[0]+PA.d[self.model.Nodes[iNode].dofs['UX'].id], self.model.Nodes[iNode].coords[1]+PA.d[self.model.Nodes[iNode].dofs['UY'].id], color='r')
        for iElement in self.model.Elements:
            x = (self.model.Elements[iElement].nodeI.coords[0]+PA.d[self.model.Elements[iElement].nodeI.dofs['UX'].id],self.model.Elements[iElement].nodeJ.coords[0]+PA.d[self.model.Elements[iElement].nodeJ.dofs['UX'].id])
            y = (self.model.Elements[iElement].nodeI.coords[1]+PA.d[self.model.Elements[iElement].nodeI.dofs['UY'].id],self.model.Elements[iElement].nodeJ.coords[1]+PA.d[self.model.Elements[iElement].nodeJ.dofs['UY'].id]) 
            ax.plot(x, y, color='r')
        for iCell in self.model.PondingLoadCells:
            x =  (self.model.PondingLoadCells[iCell].nodeI.coords[0]+PA.d[self.model.PondingLoadCells[iCell].nodeI.dofs['UX'].id],self.model.PondingLoadCells[iCell].nodeJ.coords[0]+PA.d[self.model.PondingLoadCells[iCell].nodeJ.dofs['UX'].id])
            ys = (self.model.PondingLoadCells[iCell].nodeI.coords[1]+PA.d[self.model.PondingLoadCells[iCell].nodeI.dofs['UY'].id]+self.model.PondingLoadCells[iCell].hs,self.model.PondingLoadCells[iCell].nodeJ.coords[1]+PA.d[self.model.PondingLoadCells[iCell].nodeJ.dofs['UY'].id]+self.model.PondingLoadCells[iCell].hs)
            yw = (PA.z,PA.z) 
            ax.plot(x, ys, color='c')
            ax.plot(x, yw, color='b')
        
        
        plt.show()        
       
    def lowest_point(self):
        return min([self.zi, self.zj])
        