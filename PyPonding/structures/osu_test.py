import numpy as np
import matplotlib.pyplot as plt
import math
from math import pi,pow
from PyPonding import FE
from PyPonding.structures import basic_structure

class osu_test:

    # Plan Dimensions
    L       = 48*12
    Le      = 4
    S       = 67
    Se      = 5

    # Joist and Beam Section Properties
    E       = 29000

    A_J     = 100.
    Iz_Jc   = 309.33 # 24K9
    Iz_Je   = 187.74 # 24KSP
    Iy_J    = 1000.
    GJ_J    = 1000.

    A_W     = 14.4 # W10x49 support beam
    Iz_W    = 272.
    Iy_W    = 93.4
    GJ_W    = 1000.   
    
    # Other Properties
    c       = 0
    gamma   = 62.4/1000/12**3
    
    # Analysis Options
    nele    = 20
    dof_types = ('UX','UY','UZ','RX','RY','RZ')
    na      = 4
    nb      = 4    
    
    def __init__(self,specimen):
        self.specimen = specimen
        
        if specimen == 'flat':
            self.zi  = 0
            self.zj  = 0
        elif specimen == 'pitched':
            self.zi  = 0
            self.zj  = 12
        else:
            raise Exception('Unknown specimen: %s' % specimen)
        
    def lowest_point(self):
        return min([self.zi, self.zj])         
        
    def BuildModel(self):
        model = FE.Model('OSU Test Model, %s' % self.specimen)
        
        # Define Nodes and Constraints
        model.AddNode('A00',(-self.Le,-self.S-self.Se,self.zi),self.dof_types)
        model.AddNode('B00',(-self.Le,        -self.S,self.zi),self.dof_types)
        model.AddNode('C00',(-self.Le,            0.0,self.zi),self.dof_types)
        model.AddNode('D00',(-self.Le,         self.S,self.zi),self.dof_types)
        model.AddNode('E00',(-self.Le, self.S+self.Se,self.zi),self.dof_types)
        
        for i in self.dof_types:
            model.Nodes['A00'].dofs[i].constrained = True
            model.Nodes['B00'].dofs[i].constrained = True
            model.Nodes['C00'].dofs[i].constrained = True
            model.Nodes['D00'].dofs[i].constrained = True
            model.Nodes['E00'].dofs[i].constrained = True        
        
        for i in range(1,self.nele+2):
            xi = (i-1)/self.nele
            x = xi*self.L
            z = self.zi + xi*(self.zj-self.zi)
            if self.c != 0:
                r = (self.c**2 + (self.L/2)**2)/(2*self.c)
                z += math.sqrt(r**2 - (xi*self.L - self.L/2)**2) + (self.c-r)
            model.AddNode('A%02i'%i,(x,-self.S-self.Se,z),self.dof_types)
            model.AddNode('B%02i'%i,(x,        -self.S,z),self.dof_types)
            model.AddNode('C%02i'%i,(x,            0.0,z),self.dof_types)
            model.AddNode('D%02i'%i,(x,         self.S,z),self.dof_types)
            model.AddNode('E%02i'%i,(x, self.S+self.Se,z),self.dof_types)            

            
            model.Nodes['A%02i'%i].dofs['UY'].constrained = True
            model.Nodes['A%02i'%i].dofs['RX'].constrained = True
            model.Nodes['A%02i'%i].dofs['RY'].constrained = True
            model.Nodes['A%02i'%i].dofs['RZ'].constrained = True
                                
            model.Nodes['B%02i'%i].dofs['UY'].constrained = True
            model.Nodes['B%02i'%i].dofs['RX'].constrained = True
            model.Nodes['B%02i'%i].dofs['RZ'].constrained = True
                                
            model.Nodes['C%02i'%i].dofs['UY'].constrained = True
            model.Nodes['C%02i'%i].dofs['RX'].constrained = True
            model.Nodes['C%02i'%i].dofs['RZ'].constrained = True
                                
            model.Nodes['D%02i'%i].dofs['UY'].constrained = True
            model.Nodes['D%02i'%i].dofs['RX'].constrained = True
            model.Nodes['D%02i'%i].dofs['RZ'].constrained = True
                                
            model.Nodes['E%02i'%i].dofs['UY'].constrained = True
            model.Nodes['E%02i'%i].dofs['RX'].constrained = True
            model.Nodes['E%02i'%i].dofs['RY'].constrained = True
            model.Nodes['E%02i'%i].dofs['RZ'].constrained = True            
            
            if i == 1:
                model.Nodes['A%02i'%i].dofs['UX'].constrained = True
                model.Nodes['A%02i'%i].dofs['UZ'].constrained = True
                model.Nodes['B%02i'%i].dofs['UX'].constrained = True
                model.Nodes['B%02i'%i].dofs['UZ'].constrained = True
                model.Nodes['C%02i'%i].dofs['UX'].constrained = True
                model.Nodes['D%02i'%i].dofs['UX'].constrained = True
                model.Nodes['D%02i'%i].dofs['UZ'].constrained = True
                model.Nodes['E%02i'%i].dofs['UX'].constrained = True
                model.Nodes['E%02i'%i].dofs['UZ'].constrained = True
            elif i == self.nele+1:
                model.Nodes['A%02i'%i].dofs['UX'].constrained = model.Nodes['B%02i'%i].dofs['UX'].id
                model.Nodes['A%02i'%i].dofs['UZ'].constrained = True
                model.Nodes['B%02i'%i].dofs['UZ'].constrained = True
                model.Nodes['D%02i'%i].dofs['UZ'].constrained = True
                model.Nodes['E%02i'%i].dofs['UX'].constrained = model.Nodes['D%02i'%i].dofs['UX'].id
                model.Nodes['E%02i'%i].dofs['UZ'].constrained = True
            else:
                model.Nodes['A%02i'%i].dofs['UX'].constrained = model.Nodes['B%02i'%i].dofs['UX'].id
                model.Nodes['A%02i'%i].dofs['UZ'].constrained = model.Nodes['B%02i'%i].dofs['UZ'].id
                model.Nodes['E%02i'%i].dofs['UX'].constrained = model.Nodes['D%02i'%i].dofs['UX'].id
                model.Nodes['E%02i'%i].dofs['UZ'].constrained = model.Nodes['D%02i'%i].dofs['UZ'].id
                
        model.AddNode('A%02i' %(self.nele+2),(self.L+self.Le,-self.S-self.Se,self.zj),self.dof_types)
        model.AddNode('B%02i' %(self.nele+2),(self.L+self.Le,        -self.S,self.zj),self.dof_types)
        model.AddNode('C%02i' %(self.nele+2),(self.L+self.Le,            0.0,self.zj),self.dof_types)
        model.AddNode('D%02i' %(self.nele+2),(self.L+self.Le,         self.S,self.zj),self.dof_types)
        model.AddNode('E%02i' %(self.nele+2),(self.L+self.Le, self.S+self.Se,self.zj),self.dof_types)    
       
        for i in ('UY','UZ','RX','RY','RZ'):
            model.Nodes['A%02i' %(self.nele+2)].dofs[i].constrained = True
            model.Nodes['B%02i' %(self.nele+2)].dofs[i].constrained = True
            model.Nodes['C%02i' %(self.nele+2)].dofs[i].constrained = True
            model.Nodes['D%02i' %(self.nele+2)].dofs[i].constrained = True
            model.Nodes['E%02i' %(self.nele+2)].dofs[i].constrained = True 
        
        model.Nodes['A%02i' %(self.nele+2)].dofs['UX'].constrained = model.Nodes['B%02i'%(self.nele+1)].dofs['UX'].id
        model.Nodes['B%02i' %(self.nele+2)].dofs['UX'].constrained = model.Nodes['B%02i'%(self.nele+1)].dofs['UX'].id
        model.Nodes['C%02i' %(self.nele+2)].dofs['UX'].constrained = model.Nodes['C%02i'%(self.nele+1)].dofs['UX'].id
        model.Nodes['D%02i' %(self.nele+2)].dofs['UX'].constrained = model.Nodes['D%02i'%(self.nele+1)].dofs['UX'].id
        model.Nodes['E%02i' %(self.nele+2)].dofs['UX'].constrained = model.Nodes['D%02i'%(self.nele+1)].dofs['UX'].id 
       
       
        # Define Joists
        vec_xz = (0,1,0)
        for i in range(1,self.nele+1):
            model.AddElement('B%02i'%i,'ElasticBeam3d',('B%02i'%i,'B%02i'%(i+1)),vec_xz,self.E,self.Iz_Je,self.Iy_J,self.A_J,self.GJ_J)
            model.AddElement('C%02i'%i,'ElasticBeam3d',('C%02i'%i,'C%02i'%(i+1)),vec_xz,self.E,self.Iz_Jc,self.Iy_J,self.A_J,self.GJ_J)
            model.AddElement('D%02i'%i,'ElasticBeam3d',('D%02i'%i,'D%02i'%(i+1)),vec_xz,self.E,self.Iz_Je,self.Iy_J,self.A_J,self.GJ_J)

        # Define Beams
        model.AddNode('B01W',(0.0,-self.S,self.zi),self.dof_types)
        model.AddNode('C01W',(0.0,    0.0,self.zi),self.dof_types)
        model.AddNode('D01W',(0.0, self.S,self.zi),self.dof_types)
        
        model.Nodes['B01W'].dofs['UX'].constrained = True
        model.Nodes['B01W'].dofs['UZ'].constrained = True
        model.Nodes['B01W'].dofs['RY'].constrained = True
        model.Nodes['B01W'].dofs['RZ'].constrained = True
               
        model.Nodes['C01W'].dofs['UX'].constrained = True
        model.Nodes['C01W'].dofs['UY'].constrained = True
        model.Nodes['C01W'].dofs['UZ'].constrained = model.Nodes['C01'].dofs['UZ'].id 
        model.Nodes['C01W'].dofs['RY'].constrained = True
        model.Nodes['C01W'].dofs['RZ'].constrained = True
        
        model.Nodes['D01W'].dofs['UX'].constrained = True
        model.Nodes['D01W'].dofs['UZ'].constrained = True
        model.Nodes['D01W'].dofs['RY'].constrained = True
        model.Nodes['D01W'].dofs['RZ'].constrained = True
        
        model.AddNode('B%02iW'%(self.nele+1),(self.L,-self.S,self.zj),self.dof_types)
        model.AddNode('C%02iW'%(self.nele+1),(self.L,    0.0,self.zj),self.dof_types)
        model.AddNode('D%02iW'%(self.nele+1),(self.L, self.S,self.zj),self.dof_types)
        
        model.Nodes['B%02iW'%(self.nele+1)].dofs['UX'].constrained = model.Nodes['C%02i'%(self.nele+1)].dofs['UX'].id 
        model.Nodes['B%02iW'%(self.nele+1)].dofs['UZ'].constrained = True
        model.Nodes['B%02iW'%(self.nele+1)].dofs['RY'].constrained = True
        model.Nodes['B%02iW'%(self.nele+1)].dofs['RZ'].constrained = True
        
        model.Nodes['B%02iW'%(self.nele+1)].dofs['UX'].constrained = model.Nodes['C%02i'%(self.nele+1)].dofs['UX'].id 
        model.Nodes['C%02iW'%(self.nele+1)].dofs['UY'].constrained = True
        model.Nodes['C%02iW'%(self.nele+1)].dofs['UZ'].constrained = model.Nodes['C%02i'%(self.nele+1)].dofs['UZ'].id
        model.Nodes['C%02iW'%(self.nele+1)].dofs['RY'].constrained = True
        model.Nodes['C%02iW'%(self.nele+1)].dofs['RZ'].constrained = True
        
        model.Nodes['D%02iW'%(self.nele+1)].dofs['UX'].constrained = model.Nodes['C%02i'%(self.nele+1)].dofs['UX'].id 
        model.Nodes['D%02iW'%(self.nele+1)].dofs['UZ'].constrained = True
        model.Nodes['D%02iW'%(self.nele+1)].dofs['RY'].constrained = True
        model.Nodes['D%02iW'%(self.nele+1)].dofs['RZ'].constrained = True        
        
        
        model.Nodes['B%02iW'%(self.nele+1)].dofs['UY'].constrained = True
        model.Nodes['B%02iW'%(self.nele+1)].dofs['RX'].constrained = True
        model.Nodes['C%02iW'%(self.nele+1)].dofs['RX'].constrained = True
        model.Nodes['D%02iW'%(self.nele+1)].dofs['UY'].constrained = True
        model.Nodes['D%02iW'%(self.nele+1)].dofs['RX'].constrained = True        
        
        vec_xz = (1,0,0)
        model.AddElement('W1','ElasticBeam3d',('B01W','C01W'),vec_xz,self.E,self.Iz_W ,self.Iy_W,self.A_W,self.GJ_W)
        model.AddElement('W2','ElasticBeam3d',('C01W','D01W'),vec_xz,self.E,self.Iz_W ,self.Iy_W,self.A_W,self.GJ_W)
        model.AddElement('W3','ElasticBeam3d',('B%02iW'%(self.nele+1),'C%02iW'%(self.nele+1)),vec_xz,self.E,self.Iz_W ,self.Iy_W,self.A_W,self.GJ_W)
        model.AddElement('W4','ElasticBeam3d',('C%02iW'%(self.nele+1),'D%02iW'%(self.nele+1)),vec_xz,self.E,self.Iz_W ,self.Iy_W,self.A_W,self.GJ_W)

        
        # Define ponding load cells
        for i in range(self.nele+2):
            ni = 'A%02i' % (  i)
            nj = 'A%02i' % (i+1)
            nk = 'B%02i' % (i+1)
            nl = 'B%02i' % (  i)
            model.AddPondingLoadCell(ni,'3d',(ni,nj,nk,nl),self.gamma,self.na,self.nb)
        
            ni = 'B%02i' % (  i)
            nj = 'B%02i' % (i+1)
            nk = 'C%02i' % (i+1)
            nl = 'C%02i' % (  i)
            model.AddPondingLoadCell(ni,'3d',(ni,nj,nk,nl),self.gamma,self.na,self.nb)
            
            ni = 'C%02i' % (  i)
            nj = 'C%02i' % (i+1)
            nk = 'D%02i' % (i+1)
            nl = 'D%02i' % (  i)
            model.AddPondingLoadCell(ni,'3d',(ni,nj,nk,nl),self.gamma,self.na,self.nb)
            
            ni = 'D%02i' % (  i)
            nj = 'D%02i' % (i+1)
            nk = 'E%02i' % (i+1)
            nl = 'E%02i' % (  i)
            model.AddPondingLoadCell(ni,'3d',(ni,nj,nk,nl),self.gamma,self.na,self.nb)
            
        # Define dead load
        #   no dead load - accounted for in reduced camber
        
        self.model = model
        
        
    def Reaction(self,results):
        R = self.model.Nodes['A00'].dofs['UZ'].react(results) + \
            self.model.Nodes['B00'].dofs['UZ'].react(results) + \
            self.model.Nodes['C00'].dofs['UZ'].react(results) + \
            self.model.Nodes['D00'].dofs['UZ'].react(results) + \
            self.model.Nodes['E00'].dofs['UZ'].react(results) + \
            self.model.Nodes['A01'].dofs['UZ'].react(results) + \
            self.model.Nodes['B01'].dofs['UZ'].react(results) + \
            self.model.Nodes['D01'].dofs['UZ'].react(results) + \
            self.model.Nodes['E01'].dofs['UZ'].react(results) + \
            self.model.Nodes['A%02i' %(self.nele+1)].dofs['UZ'].react(results) + \
            self.model.Nodes['B%02i' %(self.nele+1)].dofs['UZ'].react(results) + \
            self.model.Nodes['D%02i' %(self.nele+1)].dofs['UZ'].react(results) + \
            self.model.Nodes['E%02i' %(self.nele+1)].dofs['UZ'].react(results) + \
            self.model.Nodes['A%02i' %(self.nele+2)].dofs['UZ'].react(results) + \
            self.model.Nodes['B%02i' %(self.nele+2)].dofs['UZ'].react(results) + \
            self.model.Nodes['C%02i' %(self.nele+2)].dofs['UZ'].react(results) + \
            self.model.Nodes['D%02i' %(self.nele+2)].dofs['UZ'].react(results) + \
            self.model.Nodes['E%02i' %(self.nele+2)].dofs['UZ'].react(results) + \
            self.model.Nodes['B01W'].dofs['UZ'].react(results) + \
            self.model.Nodes['D01W'].dofs['UZ'].react(results) + \
            self.model.Nodes['B%02iW'%(self.nele+1)].dofs['UZ'].react(results) + \
            self.model.Nodes['D%02iW'%(self.nele+1)].dofs['UZ'].react(results)
            
        return R

    def CenterJoistMoment_and_Shear(self,results):
        x = np.empty([2*self.nele,1])
        M = np.empty([2*self.nele,1])
        V = np.empty([2*self.nele,1])
        for i in range(1,self.nele+1):
            ele_force = self.model.Elements['C%02i'%i].force(results)
            x[2*(i-1)  ] =  (i-1)*self.L/self.nele
            x[2*(i-1)+1] = i*self.L/self.nele
            M[2*(i-1)  ] =  ele_force.item(4)
            M[2*(i-1)+1] = -ele_force.item(10)
            V[2*(i-1)  ] =  ele_force.item(2)
            V[2*(i-1)+1] = -ele_force.item(8)            
        return (x,M,V)   
    
    def CenterJoistMaxMoment(self,results):
        Mmax = 0;

        for i in range(1,self.nele+1):
            ele_force = self.model.Elements['C%02i'%i].force(results)
            if abs(ele_force.item(4)) > Mmax:
                Mmax = abs(ele_force.item(4))
            if abs(ele_force.item(10)) > Mmax:
                Mmax = abs(ele_force.item(10))    

        return Mmax