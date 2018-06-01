import numpy as np
import matplotlib.pyplot as plt
import math
from math import pi,pow
from PyPonding import FE
from PyPonding.structures import roof2x3

class roof2x3b(roof2x3.roof2x3):

    def __init__(self):
        pass
        
    def BuildModel(self):
        model = FE.Model('Roof Example - Type A')
        
        # Define Nodes at Column Grids
        model.AddNode('A1',(                0.0,                          0.0,self.z_A1),self.dof_types)
        model.AddNode('A2',(          self.L_12,                          0.0,self.z_A2),self.dof_types)
        model.AddNode('A3',(self.L_12+self.L_23,                          0.0,self.z_A3),self.dof_types)
        model.AddNode('B1',(                0.0,                    self.L_AB,self.z_B1),self.dof_types)
        model.AddNode('B2',(          self.L_12,                    self.L_AB,self.z_B2),self.dof_types)
        model.AddNode('B3',(self.L_12+self.L_23,                    self.L_AB,self.z_B3),self.dof_types)
        model.AddNode('C1',(                0.0,          self.L_AB+self.L_BC,self.z_C1),self.dof_types)
        model.AddNode('C2',(          self.L_12,          self.L_AB+self.L_BC,self.z_C2),self.dof_types)
        model.AddNode('C3',(self.L_12+self.L_23,          self.L_AB+self.L_BC,self.z_C3),self.dof_types)
        model.AddNode('D1',(                0.0,self.L_AB+self.L_BC+self.L_CD,self.z_D1),self.dof_types)
        model.AddNode('D2',(          self.L_12,self.L_AB+self.L_BC+self.L_CD,self.z_D2),self.dof_types)
        model.AddNode('D3',(self.L_12+self.L_23,self.L_AB+self.L_BC+self.L_CD,self.z_D3),self.dof_types)
               
        for i in self.dof_types:
            model.Nodes['A1'].dofs[i].constrained = True
            model.Nodes['A2'].dofs[i].constrained = True
            model.Nodes['A3'].dofs[i].constrained = True
            model.Nodes['B1'].dofs[i].constrained = True
            model.Nodes['B2'].dofs[i].constrained = True
            model.Nodes['B3'].dofs[i].constrained = True
            model.Nodes['C1'].dofs[i].constrained = True
            model.Nodes['C2'].dofs[i].constrained = True
            model.Nodes['C3'].dofs[i].constrained = True
            model.Nodes['D1'].dofs[i].constrained = True
            model.Nodes['D2'].dofs[i].constrained = True
            model.Nodes['D3'].dofs[i].constrained = True
     
        # Define Edge Nodes AB1
        for i in range(self.nspaces+1):
            xi = (i/self.nspaces)
            L = self.L_AB
            n = 'AB1_%02i'%i
            x = 0
            y = xi*L
            z = self.z_A1 + xi*(self.z_B1-self.z_A1)
            model.AddNode(n,(x,y,z),self.dof_types)
            for j in self.dof_types:
                model.Nodes[n].dofs[j].constrained = True
        
        # Define Edge Nodes BC1
        for i in range(self.nspaces+1):
            xi = (i/self.nspaces)
            L = self.L_BC
            n = 'BC1_%02i'%i
            x = 0
            y = self.L_AB + xi*L
            z = self.z_B1 + xi*(self.z_C1-self.z_B1)
            model.AddNode(n,(x,y,z),self.dof_types)
            for j in self.dof_types:
                model.Nodes[n].dofs[j].constrained = True
                
        # Define Edge Nodes CD1
        for i in range(self.nspaces+1):
            xi = (i/self.nspaces)
            L = self.L_CD
            n = 'CD1_%02i'%i
            x = 0
            y = self.L_AB + self.L_BC + xi*L
            z = self.z_C1 + xi*(self.z_D1-self.z_C1)
            model.AddNode(n,(x,y,z),self.dof_types)
            for j in self.dof_types:
                model.Nodes[n].dofs[j].constrained = True                
        
        # Define Joist Girder AB2
        for i in range(self.nspaces+1):
            xi = (i/self.nspaces)
            L = self.L_AB
            n = 'AB2_%02i'%i
            c = self.c_JG
            x = self.L_12
            y = xi*L
            z = self.z_A2 + xi*(self.z_B2-self.z_A2)
            if c != 0:
                r = (c**2 + (L/2)**2)/(2*c)
                z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
            model.AddNode(n,(x,y,z),self.dof_types)
            model.Nodes[n].dofs['UX'].constrained = True
            model.Nodes[n].dofs['RY'].constrained = True
            model.Nodes[n].dofs['RZ'].constrained = True
            if i == 0:
                model.Nodes[n].dofs['UY'].constrained = True
                model.Nodes[n].dofs['UZ'].constrained = True
            if i == self.nspaces:
                model.Nodes[n].dofs['UZ'].constrained = True
        
        vec_xz = (1,0,0)
        for i in range(self.nspaces):
            model.AddElement('AB2_%02i'%i,'ElasticBeam3d',('AB2_%02i'%i,'AB2_%02i'%(i+1)),vec_xz,self.E,self.Iz_JG,self.Iy_JG,self.A_JG,self.GJ_JG)

        # Define Joist Girder BC2
        for i in range(self.nspaces+1):
            xi = (i/self.nspaces)
            L = self.L_BC
            n = 'BC2_%02i'%i
            c = self.c_JG
            x = self.L_12
            y = self.L_AB + xi*L
            z = self.z_B2 + xi*(self.z_C2-self.z_B2)
            if c != 0:
                r = (c**2 + (L/2)**2)/(2*c)
                z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
            model.AddNode(n,(x,y,z),self.dof_types)
            model.Nodes[n].dofs['UX'].constrained = True
            model.Nodes[n].dofs['RY'].constrained = True
            model.Nodes[n].dofs['RZ'].constrained = True
            if i == 0:
                model.Nodes[n].dofs['UY'].constrained = True
                model.Nodes[n].dofs['UZ'].constrained = True
            if i == self.nspaces:
                model.Nodes[n].dofs['UZ'].constrained = True
        
        vec_xz = (1,0,0)
        for i in range(self.nspaces):
            model.AddElement('BC2_%02i'%i,'ElasticBeam3d',('BC2_%02i'%i,'BC2_%02i'%(i+1)),vec_xz,self.E,self.Iz_JG,self.Iy_JG,self.A_JG,self.GJ_JG)

        # Define Joist Girder CD2
        for i in range(self.nspaces+1):
            xi = (i/self.nspaces)
            L = self.L_CD
            n = 'CD2_%02i'%i
            c = self.c_JG
            x = self.L_12
            y = self.L_AB + self.L_BC + xi*L
            z = self.z_C2 + xi*(self.z_D2-self.z_C2)
            if c != 0:
                r = (c**2 + (L/2)**2)/(2*c)
                z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
            model.AddNode(n,(x,y,z),self.dof_types)
            model.Nodes[n].dofs['UX'].constrained = True
            model.Nodes[n].dofs['RY'].constrained = True
            model.Nodes[n].dofs['RZ'].constrained = True
            if i == 0:
                model.Nodes[n].dofs['UY'].constrained = True
                model.Nodes[n].dofs['UZ'].constrained = True
            if i == self.nspaces:
                model.Nodes[n].dofs['UZ'].constrained = True
        
        vec_xz = (1,0,0)
        for i in range(self.nspaces):
            model.AddElement('CD2_%02i'%i,'ElasticBeam3d',('CD2_%02i'%i,'CD2_%02i'%(i+1)),vec_xz,self.E,self.Iz_JG,self.Iy_JG,self.A_JG,self.GJ_JG)            
            
        # Define Edge Nodes AB3
        for i in range(self.nspaces+1):
            xi = (i/self.nspaces)
            L = self.L_AB
            n = 'AB3_%02i'%i
            x = self.L_12 + self.L_23
            y = xi*L
            z = self.z_A3 + xi*(self.z_B3-self.z_A3)
            model.AddNode(n,(x,y,z),self.dof_types)
            for j in self.dof_types:
                model.Nodes[n].dofs[j].constrained = True
        
        # Define Edge Nodes BC3
        for i in range(self.nspaces+1):
            xi = (i/self.nspaces)
            L = self.L_BC
            n = 'BC3_%02i'%i
            x = self.L_12 + self.L_23
            y = self.L_AB + xi*L
            z = self.z_B3 + xi*(self.z_C3-self.z_B3)
            model.AddNode(n,(x,y,z),self.dof_types)
            for j in self.dof_types:
                model.Nodes[n].dofs[j].constrained = True
                
        # Define Edge Nodes CD3
        for i in range(self.nspaces+1):
            xi = (i/self.nspaces)
            L = self.L_CD
            n = 'CD3_%02i'%i
            x = self.L_12 + self.L_23
            y = self.L_AB + self.L_BC + xi*L
            z = self.z_C3 + xi*(self.z_D3-self.z_C3)
            model.AddNode(n,(x,y,z),self.dof_types)
            for j in self.dof_types:
                model.Nodes[n].dofs[j].constrained = True               
            
            
            
        # Define joists between 1 and 2
        for i in range(3*self.nspaces+1):
            if i == 0:
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_12
                    n = '12%02i_%02i'%(i,j)
                    x = xi*L
                    y = 0
                    z = self.z_A1 + xi*(self.z_A2-self.z_A1)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    for k in self.dof_types:
                        model.Nodes[n].dofs[k].constrained = True
                        
            elif i < self.nspaces:

                ni = 'AB1_%02i' % i
                nj = 'AB2_%02i' % i
            
                zi = model.Nodes[ni].coords[2]
                zj = model.Nodes[nj].coords[2] 
            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_12
                    n = '12%02i_%02i'%(i,j)
                    c = self.c_J
                    x = xi*L
                    y = (i/self.nspaces)*self.L_AB
                    z = zi + xi*(zj-zi)
                    if c != 0:
                        r = (c**2 + (L/2)**2)/(2*c)
                        z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    model.Nodes[n].dofs['UY'].constrained = True
                    model.Nodes[n].dofs['RX'].constrained = True
                    model.Nodes[n].dofs['RZ'].constrained = True                     
                    if j == 0:
                        model.Nodes[n].dofs['UX'].constrained = True
                        model.Nodes[n].dofs['UZ'].constrained = True
                    if j == self.nele_J:
                        model.Nodes[n].dofs['UZ'].constrained = model.Nodes[nj].dofs['UZ'].id
                
                vec_xz = (0,1,0)
                for j in range(self.nele_J):
                    model.AddElement('12%02i_%02i'%(i,j),'ElasticBeam3d',('12%02i_%02i'%(i,j),'12%02i_%02i'%(i,j+1)),vec_xz,self.E,self.Iz_J,self.Iy_J,self.A_J,self.GJ_J)
                
                
            elif i == self.nspaces:
                ni = 'B1'
                nj = 'B2'
            
                zi = model.Nodes[ni].coords[2]
                zj = model.Nodes[nj].coords[2] 
            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_12
                    n = '12%02i_%02i'%(i,j)
                    c = self.c_J
                    x = xi*L
                    y = self.L_AB
                    z = zi + xi*(zj-zi)
                    if c != 0:
                        r = (c**2 + (L/2)**2)/(2*c)
                        z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    model.Nodes[n].dofs['UY'].constrained = True
                    model.Nodes[n].dofs['RX'].constrained = True
                    model.Nodes[n].dofs['RZ'].constrained = True                     
                    if j == 0:
                        model.Nodes[n].dofs['UX'].constrained = True
                        model.Nodes[n].dofs['UZ'].constrained = True
                    if j == self.nele_J:
                        model.Nodes[n].dofs['UZ'].constrained = True
                
                vec_xz = (0,1,0)
                for j in range(self.nele_J):
                    model.AddElement('12%02i_%02i'%(i,j),'ElasticBeam3d',('12%02i_%02i'%(i,j),'12%02i_%02i'%(i,j+1)),vec_xz,self.E,self.Iz_J,self.Iy_J,self.A_J,self.GJ_J)
                        
            elif i < 2*self.nspaces:
            
                ni = 'BC1_%02i' % (i-self.nspaces)
                nj = 'BC2_%02i' % (i-self.nspaces)
            
                zi = model.Nodes[ni].coords[2]
                zj = model.Nodes[nj].coords[2]
            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_12
                    n = '12%02i_%02i'%(i,j)
                    c = self.c_J
                    x = xi*L
                    y = self.L_AB + ((i-self.nspaces)/self.nspaces)*self.L_BC
                    z = zi + xi*(zj-zi)
                    if c != 0:
                        r = (c**2 + (L/2)**2)/(2*c)
                        z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    model.Nodes[n].dofs['UY'].constrained = True
                    model.Nodes[n].dofs['RX'].constrained = True
                    model.Nodes[n].dofs['RZ'].constrained = True                     
                    if j == 0:
                        model.Nodes[n].dofs['UX'].constrained = True
                        model.Nodes[n].dofs['UZ'].constrained = True
                    if j == self.nele_J:
                        model.Nodes[n].dofs['UZ'].constrained = model.Nodes[nj].dofs['UZ'].id
                
                vec_xz = (0,1,0)
                for j in range(self.nele_J):
                    model.AddElement('12%02i_%02i'%(i,j),'ElasticBeam3d',('12%02i_%02i'%(i,j),'12%02i_%02i'%(i,j+1)),vec_xz,self.E,self.Iz_J,self.Iy_J,self.A_J,self.GJ_J)
            
            elif i == 2*self.nspaces:
                ni = 'C1'
                nj = 'C2'
            
                zi = model.Nodes[ni].coords[2]
                zj = model.Nodes[nj].coords[2] 
            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_12
                    n = '12%02i_%02i'%(i,j)
                    c = self.c_J
                    x = xi*L
                    y = self.L_AB + self.L_BC
                    z = zi + xi*(zj-zi)
                    if c != 0:
                        r = (c**2 + (L/2)**2)/(2*c)
                        z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    model.Nodes[n].dofs['UY'].constrained = True
                    model.Nodes[n].dofs['RX'].constrained = True
                    model.Nodes[n].dofs['RZ'].constrained = True                     
                    if j == 0:
                        model.Nodes[n].dofs['UX'].constrained = True
                        model.Nodes[n].dofs['UZ'].constrained = True
                    if j == self.nele_J:
                        model.Nodes[n].dofs['UZ'].constrained = True
                
                vec_xz = (0,1,0)
                for j in range(self.nele_J):
                    model.AddElement('12%02i_%02i'%(i,j),'ElasticBeam3d',('12%02i_%02i'%(i,j),'12%02i_%02i'%(i,j+1)),vec_xz,self.E,self.Iz_J,self.Iy_J,self.A_J,self.GJ_J)            
            elif i < 3*self.nspaces:
                ni = 'CD1_%02i' % (i-2*self.nspaces)
                nj = 'CD2_%02i' % (i-2*self.nspaces)
            
                zi = model.Nodes[ni].coords[2]
                zj = model.Nodes[nj].coords[2]
            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_12
                    n = '12%02i_%02i'%(i,j)
                    c = self.c_J
                    x = xi*L
                    y = self.L_AB + self.L_BC + ((i-2*self.nspaces)/self.nspaces)*self.L_CD
                    z = zi + xi*(zj-zi)
                    if c != 0:
                        r = (c**2 + (L/2)**2)/(2*c)
                        z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    model.Nodes[n].dofs['UY'].constrained = True
                    model.Nodes[n].dofs['RX'].constrained = True
                    model.Nodes[n].dofs['RZ'].constrained = True                     
                    if j == 0:
                        model.Nodes[n].dofs['UX'].constrained = True
                        model.Nodes[n].dofs['UZ'].constrained = True
                    if j == self.nele_J:
                        model.Nodes[n].dofs['UZ'].constrained = model.Nodes[nj].dofs['UZ'].id
                
                vec_xz = (0,1,0)
                for j in range(self.nele_J):
                    model.AddElement('12%02i_%02i'%(i,j),'ElasticBeam3d',('12%02i_%02i'%(i,j),'12%02i_%02i'%(i,j+1)),vec_xz,self.E,self.Iz_J,self.Iy_J,self.A_J,self.GJ_J)
            
            elif i == 3*self.nspaces:            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_12
                    n = '12%02i_%02i'%(i,j)
                    x = xi*L
                    y = self.L_AB + self.L_BC + self.L_CD
                    z = self.z_D1 + xi*(self.z_D2-self.z_D1)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    for k in self.dof_types:
                        model.Nodes[n].dofs[k].constrained = True
            
            else:
                raise Exception('Should not reach here')
        
        
        # Define joists between 2 and 3
        for i in range(3*self.nspaces+1):
            if i == 0:
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_23
                    n = '23%02i_%02i'%(i,j)
                    x = self.L_12 + xi*L
                    y = 0
                    z = self.z_A2 + xi*(self.z_A3-self.z_A2)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    for k in self.dof_types:
                        model.Nodes[n].dofs[k].constrained = True
                        
            elif i < self.nspaces:

                ni = 'AB2_%02i' % i
                nj = 'AB3_%02i' % i
            
                zi = model.Nodes[ni].coords[2]
                zj = model.Nodes[nj].coords[2] 
            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_23
                    n = '23%02i_%02i'%(i,j)
                    c = self.c_J
                    x = self.L_12 + xi*L
                    y = (i/self.nspaces)*self.L_AB
                    z = zi + xi*(zj-zi)
                    if c != 0:
                        r = (c**2 + (L/2)**2)/(2*c)
                        z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    model.Nodes[n].dofs['UY'].constrained = True
                    model.Nodes[n].dofs['RX'].constrained = True
                    model.Nodes[n].dofs['RZ'].constrained = True                     
                    if j == 0:
                        model.Nodes[n].dofs['UX'].constrained = True
                        model.Nodes[n].dofs['UZ'].constrained = model.Nodes[ni].dofs['UZ'].id
                    if j == self.nele_J:
                        model.Nodes[n].dofs['UZ'].constrained = True
                
                vec_xz = (0,1,0)
                for j in range(self.nele_J):
                    model.AddElement('23%02i_%02i'%(i,j),'ElasticBeam3d',('23%02i_%02i'%(i,j),'23%02i_%02i'%(i,j+1)),vec_xz,self.E,self.Iz_J,self.Iy_J,self.A_J,self.GJ_J)
                
                
            elif i == self.nspaces:
                ni = 'B2'
                nj = 'B3'
            
                zi = model.Nodes[ni].coords[2]
                zj = model.Nodes[nj].coords[2] 
            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_23
                    n = '23%02i_%02i'%(i,j)
                    c = self.c_J
                    x = self.L_12 + xi*L
                    y = self.L_AB
                    z = zi + xi*(zj-zi)
                    if c != 0:
                        r = (c**2 + (L/2)**2)/(2*c)
                        z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    model.Nodes[n].dofs['UY'].constrained = True
                    model.Nodes[n].dofs['RX'].constrained = True
                    model.Nodes[n].dofs['RZ'].constrained = True                     
                    if j == 0:
                        model.Nodes[n].dofs['UX'].constrained = True
                        model.Nodes[n].dofs['UZ'].constrained = True
                    if j == self.nele_J:
                        model.Nodes[n].dofs['UZ'].constrained = True
                
                vec_xz = (0,1,0)
                for j in range(self.nele_J):
                    model.AddElement('23%02i_%02i'%(i,j),'ElasticBeam3d',('23%02i_%02i'%(i,j),'23%02i_%02i'%(i,j+1)),vec_xz,self.E,self.Iz_J,self.Iy_J,self.A_J,self.GJ_J)
                        
            elif i < 2*self.nspaces:
            
                ni = 'BC2_%02i' % (i-self.nspaces)
                nj = 'BC3_%02i' % (i-self.nspaces)
            
                zi = model.Nodes[ni].coords[2]
                zj = model.Nodes[nj].coords[2]
            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_23
                    n = '23%02i_%02i'%(i,j)
                    c = self.c_J
                    x = self.L_12 + xi*L
                    y = self.L_AB + ((i-self.nspaces)/self.nspaces)*self.L_BC
                    z = zi + xi*(zj-zi)
                    if c != 0:
                        r = (c**2 + (L/2)**2)/(2*c)
                        z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    model.Nodes[n].dofs['UY'].constrained = True
                    model.Nodes[n].dofs['RX'].constrained = True
                    model.Nodes[n].dofs['RZ'].constrained = True                     
                    if j == 0:
                        model.Nodes[n].dofs['UX'].constrained = True
                        model.Nodes[n].dofs['UZ'].constrained = model.Nodes[ni].dofs['UZ'].id
                    if j == self.nele_J:
                        model.Nodes[n].dofs['UZ'].constrained = True
                
                vec_xz = (0,1,0)
                for j in range(self.nele_J):
                    model.AddElement('23%02i_%02i'%(i,j),'ElasticBeam3d',('23%02i_%02i'%(i,j),'23%02i_%02i'%(i,j+1)),vec_xz,self.E,self.Iz_J,self.Iy_J,self.A_J,self.GJ_J)
            
            elif i == 2*self.nspaces:
                ni = 'C2'
                nj = 'C3'
            
                zi = model.Nodes[ni].coords[2]
                zj = model.Nodes[nj].coords[2] 
            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_23
                    n = '23%02i_%02i'%(i,j)
                    c = self.c_J
                    x = self.L_12 + xi*L
                    y = self.L_AB + self.L_BC
                    z = zi + xi*(zj-zi)
                    if c != 0:
                        r = (c**2 + (L/2)**2)/(2*c)
                        z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    model.Nodes[n].dofs['UY'].constrained = True
                    model.Nodes[n].dofs['RX'].constrained = True
                    model.Nodes[n].dofs['RZ'].constrained = True                     
                    if j == 0:
                        model.Nodes[n].dofs['UX'].constrained = True
                        model.Nodes[n].dofs['UZ'].constrained = True
                    if j == self.nele_J:
                        model.Nodes[n].dofs['UZ'].constrained = True
                
                vec_xz = (0,1,0)
                for j in range(self.nele_J):
                    model.AddElement('23%02i_%02i'%(i,j),'ElasticBeam3d',('23%02i_%02i'%(i,j),'23%02i_%02i'%(i,j+1)),vec_xz,self.E,self.Iz_J,self.Iy_J,self.A_J,self.GJ_J)            
            elif i < 3*self.nspaces:
                ni = 'CD2_%02i' % (i-2*self.nspaces)
                nj = 'CD3_%02i' % (i-2*self.nspaces)
            
                zi = model.Nodes[ni].coords[2]
                zj = model.Nodes[nj].coords[2]
            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_23
                    n = '23%02i_%02i'%(i,j)
                    c = self.c_J
                    x = self.L_12 + xi*L
                    y = self.L_AB + self.L_BC + ((i-2*self.nspaces)/self.nspaces)*self.L_CD
                    z = zi + xi*(zj-zi)
                    if c != 0:
                        r = (c**2 + (L/2)**2)/(2*c)
                        z += math.sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    model.Nodes[n].dofs['UY'].constrained = True
                    model.Nodes[n].dofs['RX'].constrained = True
                    model.Nodes[n].dofs['RZ'].constrained = True                     
                    if j == 0:
                        model.Nodes[n].dofs['UX'].constrained = True
                        model.Nodes[n].dofs['UZ'].constrained = model.Nodes[ni].dofs['UZ'].id
                    if j == self.nele_J:
                        model.Nodes[n].dofs['UZ'].constrained = True
                
                vec_xz = (0,1,0)
                for j in range(self.nele_J):
                    model.AddElement('23%02i_%02i'%(i,j),'ElasticBeam3d',('23%02i_%02i'%(i,j),'23%02i_%02i'%(i,j+1)),vec_xz,self.E,self.Iz_J,self.Iy_J,self.A_J,self.GJ_J)
            
            elif i == 3*self.nspaces:            
                for j in range(self.nele_J+1):
                    xi = (j/self.nele_J)
                    L = self.L_23
                    n = '23%02i_%02i'%(i,j)
                    x = self.L_12 + xi*L
                    y = self.L_AB + self.L_BC + self.L_CD
                    z = self.z_D2 + xi*(self.z_D3-self.z_D2)
                    model.AddNode(n,(x,y,z),self.dof_types)
                    for k in self.dof_types:
                        model.Nodes[n].dofs[k].constrained = True
            
            else:
                raise Exception('Should not reach here')        

        # Define ponding load cells between 1 and 2
        for i in range(3*self.nspaces):
            for j in range(self.nele_J):
                ni = '12%02i_%02i' % (  i,  j)
                nj = '12%02i_%02i' % (  i,j+1)
                nk = '12%02i_%02i' % (i+1,j+1)
                nl = '12%02i_%02i' % (i+1,  j)
                model.AddPondingLoadCell(ni,'3d',(ni,nj,nk,nl),self.alpha*self.LF_P*self.gamma,self.na,self.nb)
                model.PondingLoadCells[ni].gammas = self.alpha*self.LF_S1*self.gammas
                model.PondingLoadCells[ni].hs = self.hs
        
        # Define ponding load cells between 2 and 3
        for i in range(3*self.nspaces):
            for j in range(self.nele_J):
                ni = '23%02i_%02i' % (  i,  j)
                nj = '23%02i_%02i' % (  i,j+1)
                nk = '23%02i_%02i' % (i+1,j+1)
                nl = '23%02i_%02i' % (i+1,  j)
                model.AddPondingLoadCell(ni,'3d',(ni,nj,nk,nl),self.alpha*self.LF_P*self.gamma,self.na,self.nb)
                model.PondingLoadCells[ni].gammas = self.alpha*self.LF_S1*self.gammas
                model.PondingLoadCells[ni].hs = self.hs        
        
        # Define dead and snow load
        Pi1D = -self.wd*(self.L_AB/self.nspaces)*(self.L_12/self.nele_J)
        Pi2D = -self.wd*(self.L_AB/self.nspaces)*(self.L_23/self.nele_J)
        Pi3D = -self.wd*(self.L_BC/self.nspaces)*(self.L_12/self.nele_J)
        Pi4D = -self.wd*(self.L_BC/self.nspaces)*(self.L_23/self.nele_J)
        Pi5D = -self.wd*(self.L_CD/self.nspaces)*(self.L_12/self.nele_J)
        Pi6D = -self.wd*(self.L_CD/self.nspaces)*(self.L_23/self.nele_J)
        Pi1S = -self.gammas*self.hs*(self.L_AB/self.nspaces)*(self.L_12/self.nele_J)
        Pi2S = -self.gammas*self.hs*(self.L_AB/self.nspaces)*(self.L_23/self.nele_J)
        Pi3S = -self.gammas*self.hs*(self.L_BC/self.nspaces)*(self.L_12/self.nele_J)
        Pi4S = -self.gammas*self.hs*(self.L_BC/self.nspaces)*(self.L_23/self.nele_J)
        Pi5S = -self.gammas*self.hs*(self.L_CD/self.nspaces)*(self.L_12/self.nele_J)
        Pi6S = -self.gammas*self.hs*(self.L_CD/self.nspaces)*(self.L_23/self.nele_J)
        
        for i in range(3*self.nspaces+1):
            for j in range(0,self.nele_J+1):
                if i == 0:
                    if j == 0 or j == self.nele_J:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi1D/4
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi2D/4
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi1S/4
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi2S/4
                    else:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi1D/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi2D/2
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi1S/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi2S/2
                elif i < self.nspaces:               
                    if j == 0 or j == self.nele_J:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi1D/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi2D/2
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi1S/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi2S/2
                    else:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi1D
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi2D
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi1S
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi2S
                elif i == self.nspaces:
                    if j == 0 or j == self.nele_J:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi1D/4 + Pi3D/4
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi2D/4 + Pi4D/4
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi1S/4 + Pi3S/4
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi2S/4 + Pi4S/4
                    else:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi1D/2 + Pi3D/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi2D/2 + Pi4D/2
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi1S/2 + Pi3S/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi2S/2 + Pi4S/2
                elif i < 2*self.nspaces:
                    if j == 0 or j == self.nele_J:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi3D/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi4D/2
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi3S/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi4S/2
                    else:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi3D
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi4D
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi3S
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi4S
                elif i == 2*self.nspaces:
                    if j == 0 or j == self.nele_J:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi3D/4 + Pi5D/4
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi4D/4 + Pi6D/4
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi3S/4 + Pi5S/4
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi4S/4 + Pi6S/4
                    else:                                                                 
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi3D/2 + Pi5D/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi4D/2 + Pi6D/2
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi3S/2 + Pi5S/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi4S/2 + Pi6S/2
                elif i < 3*self.nspaces:
                    if j == 0 or j == self.nele_J:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi5D/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi6D/2
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi5S/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi6S/2
                    else:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi5D
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi6D
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi5S
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi6S
                elif i == 3*self.nspaces:
                    if j == 0 or j == self.nele_J:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi5D/4
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi6D/4
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi5S/4
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi6S/4
                    else:
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi5D/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['DEAD'] = Pi6D/2
                        model.Nodes['12%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi5S/2
                        model.Nodes['23%02i_%02i'%(i,j)].dofs['UZ'].loads['SNOW'] = Pi6S/2
                else:
                    raise Exception('Should not reach here') 
        
        self.model = model
        
        
    def ColumnReaction(self,column,results):
        if column == 'B2':
            R = self.model.Nodes['B2'].dofs['UZ'].react(results) + \
                self.model.Nodes['AB2_%02i'%(self.nele_J)].dofs['UZ'].react(results) + \
                self.model.Nodes['23%02i_00'%(self.nspaces)].dofs['UZ'].react(results) + \
                self.model.Nodes['BC2_00'].dofs['UZ'].react(results) + \
                self.model.Nodes['12%02i_%02i'%(self.nspaces,self.nele_J)].dofs['UZ'].react(results)
        elif column == 'C2':
            R = self.model.Nodes['C2'].dofs['UZ'].react(results) + \
                self.model.Nodes['BC2_%02i'%(self.nele_J)].dofs['UZ'].react(results) + \
                self.model.Nodes['23%02i_00'%(2*self.nspaces)].dofs['UZ'].react(results) + \
                self.model.Nodes['CD2_00'].dofs['UZ'].react(results) + \
                self.model.Nodes['12%02i_%02i'%(2*self.nspaces,self.nele_J)].dofs['UZ'].react(results)
        else:
            raise Exception('Unknown column')
            
        return R/self.alpha
        
    def Strength_Ratio(self,results):
        SR = 0
        SR_note = ''
        
        for i in range(1,3*self.nspaces):
            joist = '12%02i'%(i)
            (iSR,iSR_note) = self.Strength_Ratio_Joist(results,joist)
            if iSR > SR:
                SR = iSR
                SR_note = iSR_note + ' (Joist ' + joist + ')'
    
            joist = '23%02i'%(i)
            (iSR,iSR_note) = self.Strength_Ratio_Joist(results,joist)
            if iSR > SR:
                SR = iSR
                SR_note = iSR_note + ' (Joist ' + joist + ')'
                
        (iSR,iSR_note) = self.Strength_Ratio_Joist_Girder(results,'AB2')
        if iSR > SR:
            SR = iSR
            SR_note = iSR_note + ' (Joist Girder AB2)'
        
        (iSR,iSR_note) = self.Strength_Ratio_Joist_Girder(results,'BC2')
        if iSR > SR:
            SR = iSR
            SR_note = iSR_note + ' (Joist Girder BC2)'
            
        (iSR,iSR_note) = self.Strength_Ratio_Joist_Girder(results,'CD2')
        if iSR > SR:
            SR = iSR
            SR_note = iSR_note + ' (Joist Girder CD2)'                 
        
        return (SR,SR_note) 
    
    def Moment_and_Shear_Joist(self,results,joist):
        if joist[:2]=='12':
            L = self.L_12
        elif joist[:2]=='23':
            L = self.L_23          
        else:
            raise Exception('Unknown joist: %s'%(joist))

        x = np.empty([2*self.nele_J,1])
        M = np.empty([2*self.nele_J,1])
        V = np.empty([2*self.nele_J,1])
        for i in range(self.nele_J):
            ele_force = self.model.Elements['%s_%02i'%(joist,i)].force(results)/self.alpha
            x[2*i+0] = i*L/self.nele_J
            x[2*i+1] = (i+1)*L/self.nele_J
            M[2*i+0] =  ele_force.item(4)
            M[2*i+1] = -ele_force.item(10)
            V[2*i+0] =  ele_force.item(2)
            V[2*i+1] = -ele_force.item(8)            
        return (x,M,V)    
    
    def Strength_Ratio_Joist(self,results,joist):
        SR = 0
        SR_note = ''
               
        if joist[:2]=='12':
            L = self.L_12
        elif joist[:2]=='23':
            L = self.L_12       
        else:
            raise Exception('Unknown joist: %s'%(joist))
        
        # Reaction @todo - can I compare to reactions?
        # (Ri,Rj) = self.Reactions(results)
        iR = self.w_J*L/2
        # if Ri/iR > SR:
        #     SR = Ri/iR
        #     SR_note = 'Reaction at I end'
        # 
        # if Rj/iR > SR:
        #     SR = Rj/iR
        #     SR_note = 'Reaction at J end'
        
        for i in range(self.nele_J):
            ele_force = self.model.Elements['%s_%02i'%(joist,i)].force(results)/self.alpha
            
            # Moment at I-end
            M = ele_force.item(4)
            x = i*L/self.nele_J
            iM = self.w_J*x*(L-x)/2
            if iM > 0:
                if M/iM > SR:
                    SR = M/iM
                    SR_note = 'Moment at x/L = %0.3f' % (x/L)            
            
            # Moment at J-end
            M = -ele_force.item(10)
            x = (i+1)*L/self.nele_J
            iM = self.w_J*x*(L-x)/2
            if iM > 0:
                if M/iM > SR:
                    SR = M/iM
                    SR_note = 'Moment at x/L = %0.3f' % (x/L)            
            
            # Shear at middle of element
            V = ele_force.item(2)        
            x = (i+0.5)*L/self.nele_J
            if x < 0.5*L:
                if V >= 0:
                    iV = max(self.w_J*(L/2-x),0.25*iR)
                    if V/iV > SR:
                        SR = V/iV
                        SR_note = 'Shear at x/L = %0.3f' % (x/L)
                else:
                    iV = -self.shear_reversal_strength*iR
                    if V/iV > SR:
                        SR = V/iV
                        SR_note = 'Shear Reversal at x/L = %0.3f' % (x/L)
            else:
                if V >= 0:
                    iV = self.shear_reversal_strength*iR
                    if V/iV > SR:
                        SR = V/iV
                        SR_note = 'Shear Reversal at x/L = %0.3f' % (x/L)
                else:
                    iV = min(self.w_J*(L/2-x),-0.25*iR)             
                    if V/iV > SR:
                        SR = V/iV
                        SR_note = 'Shear at x/L = %0.3f' % (x/L)
        
        return (SR,SR_note)        
        
    def Moment_and_Shear_Joist_Girder(self,results,Joist_Girder):
        if Joist_Girder == 'AB2':
            L = self.L_AB
        elif Joist_Girder == 'BC2':
            L = self.L_BC
        elif Joist_Girder == 'CD2':
            L = self.L_CD
        else:
            raise Exception('Unknown Joist Girder: %s'%(Joist_Girder))

        x = np.empty([2*self.nspaces,1])
        M = np.empty([2*self.nspaces,1])
        V = np.empty([2*self.nspaces,1])
        for i in range(self.nspaces):
            ele_force = self.model.Elements['%s_%02i'%(Joist_Girder,i)].force(results)/self.alpha
            x[2*i+0] = i*L/self.nspaces
            x[2*i+1] = (i+1)*L/self.nspaces
            M[2*i+0] = -ele_force.item(3)
            M[2*i+1] =  ele_force.item(9)
            V[2*i+0] =  ele_force.item(2)
            V[2*i+1] = -ele_force.item(8)            
        return (x,M,V)         
        
    def Strength_Ratio_Joist_Girder(self,results,Joist_Girder):
        SR = 0
        SR_note = ''
               
        if Joist_Girder == 'AB2':
            L = self.L_AB
        elif Joist_Girder == 'BC2':
            L = self.L_BC
        elif Joist_Girder == 'CD2':
            L = self.L_CD
        else:
            raise Exception('Unknown Joist Girder: %s'%(Joist_Girder))
        
        iR = self.P_JG*(self.nspaces-1)/2
        
        for i in range(self.nspaces):
            ele_force = self.model.Elements['%s_%02i'%(Joist_Girder,i)].force(results)/self.alpha
            
            # Moment at I-end
            M = -ele_force.item(3)
            x = i*L/self.nspaces
            iM = self.P_JG*L*i*(self.nspaces-i)/(2*self.nspaces)
            if iM > 0:
                if M/iM > SR:
                    SR = M/iM
                    SR_note = 'Moment at x/L = %0.3f' % (x/L)            
            
            # Moment at J-end
            M = ele_force.item(9)
            x = (i+1)*L/self.nspaces
            iM = self.P_JG*L*(i+1)*(self.nspaces-(i+1))/(2*self.nspaces)
            if iM > 0:      
                if M/iM > SR:
                    SR = M/iM
                    SR_note = 'Moment at x/L = %0.3f' % (x/L)            
            
            # Shear at I-end (shear capacity and demand are constant so no need to check J-end)
            V = ele_force.item(2)
            x = i*L/self.nspaces
            if x < 0.5*L:
                iV = max(self.P_JG*((self.nspaces-1)/2-i),0.25*iR)
                if V < 0:
                    iV = -0.25*iV
            else:
                iV = min(self.P_JG*((self.nspaces-1)/2-i),-0.25*iR)
                if V > 0:
                    iV = -0.25*iV
            if V/iV > SR:
                SR = V/iV
                SR_note = 'Shear at x/L = %0.3f' % (x/L)

        return (SR,SR_note) 
       