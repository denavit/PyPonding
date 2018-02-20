import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from math import hypot,sqrt

import sys
if(sys.version < '3'):
    raise Exception('This script requires Python 3')

class Model:
    ndof = 0
    use_sparse_matrix_solver = True
    
    def __init__(self,name):
        self.name = name
        self.Nodes = dict()
        self.Elements = dict()
        self.PondingLoadCells = dict()
        
    def AddNode(self,id,coords,dof_types):
        # possible dof_types: 'DX','DY','DZ','RX','RY','RZ'
        if id in self.Nodes:
            raise Exception('Input Error - node "%s" is already defined' % id)
        dofs = dict()
        for i in range(len(dof_types)):
            dofs[dof_types[i]] = dof(self.ndof)
            self.ndof += 1
        self.Nodes[id] = Node(id,coords,dofs)
        
    def AddElement(self,id,type,nodes,*args):
        if id in self.Elements:
            raise Exception('Input Error - element "%s" is already defined' % id)
        if type.lower() == 'elasticbeam2d':
            nodei = self.Nodes[nodes[0]]
            nodej = self.Nodes[nodes[1]]
            self.Elements[id] = ElasticBeam2d(id,nodei,nodej,*args)
        elif type.lower() == 'elasticbeam3d':
            nodei = self.Nodes[nodes[0]]
            nodej = self.Nodes[nodes[1]]
            self.Elements[id] = ElasticBeam3d(id,nodei,nodej,*args)
        else:
            raise Exception('Input Error - unknown element type (%s)' % type)

    def AddPondingLoadCell(self,id,type,nodes,*args):
        if id in self.PondingLoadCells:
            raise Exception('Input Error - ponding load cell "%s" is already defined' % id)
        if type.lower() == '2d':
            nodei = self.Nodes[nodes[0]]
            nodej = self.Nodes[nodes[1]]
            self.PondingLoadCells[id] = PondingLoadCell2d(id,nodei,nodej,*args)
        elif type.lower() == '3d':
            nodei = self.Nodes[nodes[0]]
            nodej = self.Nodes[nodes[1]]
            nodek = self.Nodes[nodes[2]]
            nodel = self.Nodes[nodes[3]]
            self.PondingLoadCells[id] = PondingLoadCell3d(id,nodei,nodej,nodek,nodel,*args)
        else:
            raise Exception('Input Error - unknown ponding load cell type (%s)' % type)
            
    def GetGlobalStiffnessMatrix(self):
        # assemble stiffness matrix
        K = np.zeros((self.ndof,self.ndof))
        for iElement in self.Elements:
            iK    = self.Elements[iElement].StiffnessMatrix()
            idofs = self.Elements[iElement].get_dof_ids()
            for i in range(len(idofs)):
                for j in range(len(idofs)):
                    K[idofs[i],idofs[j]] += iK[i,j]
        return K
        
    def GetNodalForceVector(self,load_factors):
        # assemble force vector
        f = np.zeros(self.ndof)
        for load_pattern in load_factors:
            lf = load_factors[load_pattern]
            for iNode in self.Nodes:
                for idof in self.Nodes[iNode].dofs:
                    if load_pattern in self.Nodes[iNode].dofs[idof].loads:
                        f[self.Nodes[iNode].dofs[idof].id] += lf*self.Nodes[iNode].dofs[idof].loads[load_pattern]
        return f
 
    def GetPondingForceVector(self,d,z):
        # assemble force vector
        f = np.zeros(self.ndof)
        for iCell in self.PondingLoadCells:
            ipf   = self.PondingLoadCells[iCell].get_load_vector(d,z)
            idofs = self.PondingLoadCells[iCell].get_dof_ids()
            for i in range(len(idofs)):
                f[idofs[i]] += ipf[i]
        return f
        
    def GetPondingVolume(self,d,z):
        V = 0 
        dVdz = 0
        for iCell in self.PondingLoadCells:
            ires = self.PondingLoadCells[iCell].get_volume(d,z)
            V    += ires[0]
            dVdz += ires[1]
        return (V,dVdz)
 
    def SolveForDisp(self,K,f):
        # identidy free dofs and equal dof constraints
        free_dofs = list()
        dof_map = dict()
        equal_dof_constraints = dict()
        for iNode in self.Nodes:
            for idof in self.Nodes[iNode].dofs:
                if self.Nodes[iNode].dofs[idof].constrained == True:
                    pass
                elif self.Nodes[iNode].dofs[idof].constrained == False:
                    free_dofs.append(self.Nodes[iNode].dofs[idof].id)
                    dof_map[self.Nodes[iNode].dofs[idof].id] = len(free_dofs) - 1
                else:
                    equal_dof_constraints[self.Nodes[iNode].dofs[idof].id] = self.Nodes[iNode].dofs[idof].constrained
                    
        # Assemble the free dofs
        K_free = np.zeros((len(free_dofs),len(free_dofs)))
        f_free = np.zeros(len(free_dofs))
        for i in range(len(free_dofs)):
            for j in range(len(free_dofs)):
                K_free[i,j] = K[free_dofs[i],free_dofs[j]]
            f_free[i] = f[free_dofs[i]]
                
        # Add in stiffness from equal dof constraints
        for i in equal_dof_constraints:
            for j in free_dofs:
                K_free[dof_map[equal_dof_constraints[i]],dof_map[j]] += K[i,j]       
                K_free[dof_map[j],dof_map[equal_dof_constraints[i]]] += K[j,i] # @todo - Check this. 
            K_free[dof_map[equal_dof_constraints[i]],dof_map[equal_dof_constraints[i]]] += K[i,i]
            f_free[dof_map[equal_dof_constraints[i]]] += f[i]
        
        # Solve the system of equations
        if self.use_sparse_matrix_solver:
            K_free = scipy.sparse.csc_matrix(K_free)
            d_free = scipy.sparse.linalg.spsolve(K_free,f_free)
        else:
            d_free = np.linalg.solve(K_free,f_free)
        
        # assemble the entire deformation vector
        d = np.zeros(self.ndof)
        for i in range(len(free_dofs)):
            d[free_dofs[i]] = d_free[i]
        for i in equal_dof_constraints:
            d[i] = d_free[dof_map[equal_dof_constraints[i]]]
        
        # compute reaction vector
        r = np.dot(K,d) - f
        for i in equal_dof_constraints:
            r[equal_dof_constraints[i]] += r[i]
            r[i] = 0
        
        # Return displacement and reaction vector        
        return (d,r)

    def StoreAnalysis(self):
        # identidy free dofs and equal dof constraints
        free_dofs = list()
        dof_map = dict()
        equal_dof_constraints = dict()
        for iNode in self.Nodes:
            for idof in self.Nodes[iNode].dofs:
                if self.Nodes[iNode].dofs[idof].constrained == True:
                    pass
                elif self.Nodes[iNode].dofs[idof].constrained == False:
                    free_dofs.append(self.Nodes[iNode].dofs[idof].id)
                    dof_map[self.Nodes[iNode].dofs[idof].id] = len(free_dofs) - 1
                else:
                    equal_dof_constraints[self.Nodes[iNode].dofs[idof].id] = self.Nodes[iNode].dofs[idof].constrained

                    
        # Get Global Stiffness Matrix
        K = self.GetGlobalStiffnessMatrix()
                    
        # Assemble the free dofs
        K_free = np.zeros((len(free_dofs),len(free_dofs)))
        for i in range(len(free_dofs)):
            for j in range(len(free_dofs)):
                K_free[i,j] = K[free_dofs[i],free_dofs[j]]
                
        # Add in stiffness from equal dof constraints
        for i in equal_dof_constraints:
            for j in free_dofs:
                K_free[dof_map[equal_dof_constraints[i]],dof_map[j]] += K[i,j]       
                K_free[dof_map[j],dof_map[equal_dof_constraints[i]]] += K[j,i] # @todo - Check this. 
            K_free[dof_map[equal_dof_constraints[i]],dof_map[equal_dof_constraints[i]]] += K[i,i]
        
        # Solve the system of equations
        if self.use_sparse_matrix_solver:
            K_free = scipy.sparse.csc_matrix(K_free)
            inv_K_free = scipy.sparse.linalg.splu(K_free)
        else:
            inv_K_free = np.linalg.tensorinv(K_free)
        
        # Store Analysis Matricies
        self.stored_free_dofs = free_dofs
        self.stored_dof_map = dof_map
        self.stored_equal_dof_constraints = equal_dof_constraints
        self.stored_inv_K_free = inv_K_free
        self.stored_K = K
        
    def SolveForDispWithStored(self,f):
        
        # Assemble the free dofs
        f_free = np.zeros(len(self.stored_free_dofs))
        for i in range(len(self.stored_free_dofs)):
            f_free[i] = f[self.stored_free_dofs[i]]
                
        # Add in stiffness from equal dof constraints
        for i in self.stored_equal_dof_constraints:
            f_free[self.stored_dof_map[self.stored_equal_dof_constraints[i]]] += f[i]
        
        # Solve the system of equations
        d_free = self.stored_inv_K_free.solve(f_free)
        
        # assemble the entire deformation vector
        d = np.zeros(self.ndof)
        for i in range(len(self.stored_free_dofs)):
            d[self.stored_free_dofs[i]] = d_free[i]
        for i in self.stored_equal_dof_constraints:
            d[i] = d_free[self.stored_dof_map[self.stored_equal_dof_constraints[i]]]
        
        # compute reaction vector
        r = self.stored_K.dot(d) - f
        for i in self.stored_equal_dof_constraints:
            r[self.stored_equal_dof_constraints[i]] += r[i]
            r[i] = 0
        
        # Return displacement and reaction vector        
        return (d,r)
        
    def PlotModel3d(self):
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for iNode in self.Nodes:
            ax.scatter(self.Nodes[iNode].coords[0], self.Nodes[iNode].coords[1], self.Nodes[iNode].coords[2])
        for iElement in self.Elements:
            x = (self.Elements[iElement].nodeI.coords[0],self.Elements[iElement].nodeJ.coords[0])
            y = (self.Elements[iElement].nodeI.coords[1],self.Elements[iElement].nodeJ.coords[1]) 
            z = (self.Elements[iElement].nodeI.coords[2],self.Elements[iElement].nodeJ.coords[2])
            ax.plot(x, y, z)
        plt.show()

    def PlotModel2d(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for iNode in self.Nodes:
            ax.scatter(self.Nodes[iNode].coords[0], self.Nodes[iNode].coords[1])
        for iElement in self.Elements:
            x = (self.Elements[iElement].nodeI.coords[0],self.Elements[iElement].nodeJ.coords[0])
            y = (self.Elements[iElement].nodeI.coords[1],self.Elements[iElement].nodeJ.coords[1]) 
            ax.plot(x, y)
        plt.show()
        
    def TotalReaction(self,results):
        react = dict()
        for iNode in self.Nodes:
            for idof in self.Nodes[iNode].dofs:
                if idof in react.keys():
                    react[idof] += self.Nodes[iNode].dofs[idof].react(results)
                else:
                    react[idof] = self.Nodes[iNode].dofs[idof].react(results)
        return react
    
    def print_nodes(self,filename):
        f = open(filename,'w')
        f.write('ID,x,y,z\n')
        for iNode in self.Nodes:
            f.write('%s,%g,%g,%g\n'%(self.Nodes[iNode].id,self.Nodes[iNode].coords[0],self.Nodes[iNode].coords[1],self.Nodes[iNode].coords[2]))
        f.close()

    def print_dofs(self,filename,results):
        f = open(filename,'w')
        f.write('Node ID,dof type,dof id,constrained,dead load,disp,react\n')
        for iNode in self.Nodes:
            for idof in self.Nodes[iNode].dofs:
                if 'DEAD' in self.Nodes[iNode].dofs[idof].loads:
                    p = self.Nodes[iNode].dofs[idof].loads['DEAD']
                else:
                    p = 0
                f.write('%s,%s,%s,%i,%g,%g,%g\n'%(self.Nodes[iNode].id,idof,self.Nodes[iNode].dofs[idof].id,self.Nodes[iNode].dofs[idof].constrained,p, \
                    self.Nodes[iNode].dofs[idof].disp(results),self.Nodes[iNode].dofs[idof].react(results)))
        f.close()
        
class Node:
    def __init__(self,id,coords,dofs):
        self.id = id
        self.coords = coords
        self.dofs = dofs
        
class dof:
    def __init__(self,id):
        self.id = id;
        self.constrained = False;
        self.loads = dict();
    
    def disp(self,results):
        return results.d[self.id]

    def react(self,results):
        return results.r[self.id]
        
class ElasticBeam2d:    
    def __init__(self,id,nodeI,nodeJ,E,I,A):
        self.id = id
        self.nodeI = nodeI
        self.nodeJ = nodeJ
        self.E = E
        self.I = I
        self.A = A
    
    def get_dof_ids(self):
        dofs = [self.nodeI.dofs['UX'].id, 
                self.nodeI.dofs['UY'].id, 
                self.nodeI.dofs['RZ'].id, 
                self.nodeJ.dofs['UX'].id, 
                self.nodeJ.dofs['UY'].id, 
                self.nodeJ.dofs['RZ'].id]
        return dofs
        
    def StiffnessMatrix(self):
        nIx = self.nodeI.coords[0]
        nIy = self.nodeI.coords[1]
        nJx = self.nodeJ.coords[0]
        nJy = self.nodeJ.coords[1]
        
        L  = hypot(nJx-nIx,nJy-nIy);
        lx = (nJx-nIx)/L;
        ly = (nJy-nIy)/L;
        
        k1 =    self.E*self.A/L;
        k2 = 12*self.E*self.I/(L*L*L);
        k3 =  6*self.E*self.I/(L*L);
        k4 =  4*self.E*self.I/L;
        k5 =  2*self.E*self.I/L;
        
        Kp = np.mat([[k1,0,0,-k1,0,0],
                     [0,k2,k3,0,-k2,k3],
                     [0,k3,k4,0,-k3,k5],
                     [-k1,0,0,k1,0,0],
                     [0,-k2,-k3,0,k2,-k3],
                     [0,k3,k5,0,-k3,k4]]);

        T  = np.mat([[ lx,ly,0,  0, 0,0],
                     [-ly,lx,0,  0, 0,0],
                     [  0, 0,1,  0, 0,0],
                     [  0, 0,0, lx,ly,0],
                     [  0, 0,0,-ly,lx,0],
                     [  0, 0,0,  0, 0,1]]);        

        K = np.transpose(T)*Kp*T
        return K
        
    def disp(self,results):
        dofs = self.get_dof_ids()
        d = np.zeros(len(dofs))
        for i in range(len(dofs)):
            d[i] = results.d[dofs[i]]
        return d        
        
    def force(self,results):
        K_ele = self.StiffnessMatrix()
        d_ele = self.disp(results)
        return np.dot(K_ele,d_ele)

class ElasticBeam3d:
    release_Mzi = False
    release_Mzj = False
    release_Myi = False
    release_Myj = False
    release_Ti  = False
    release_Tj  = False
    
    def __init__(self,id,nodeI,nodeJ,vec_xz,E,Iz,Iy,A,GJ):
        self.id = id
        self.nodeI  = nodeI
        self.nodeJ  = nodeJ
        self.vec_xz = vec_xz # A vector in the xz plane of the local coordinate system
        self.E  = E
        self.Iz = Iz
        self.Iy = Iy
        self.A  = A
        self.GJ = GJ
    
    def get_dof_ids(self):
        dofs = [self.nodeI.dofs['UX'].id, 
                self.nodeI.dofs['UY'].id, 
                self.nodeI.dofs['UZ'].id,
                self.nodeI.dofs['RX'].id, 
                self.nodeI.dofs['RY'].id,
                self.nodeI.dofs['RZ'].id,
                self.nodeJ.dofs['UX'].id, 
                self.nodeJ.dofs['UY'].id, 
                self.nodeJ.dofs['UZ'].id,
                self.nodeJ.dofs['RX'].id,
                self.nodeJ.dofs['RY'].id,
                self.nodeJ.dofs['RZ'].id]
        return dofs
        
    def StiffnessMatrix(self,do_releases = True):
        nIx = self.nodeI.coords[0]
        nIy = self.nodeI.coords[1]
        nIz = self.nodeI.coords[2]
        nJx = self.nodeJ.coords[0]
        nJy = self.nodeJ.coords[1]
        nJz = self.nodeJ.coords[2]
        
        L  = sqrt((nJx-nIx)*(nJx-nIx)+(nJy-nIy)*(nJy-nIy)+(nJz-nIz)*(nJz-nIz));
        
        k1  =    self.E*self.A/L
        k2  = 12*self.E*self.Iz/(L*L*L)
        k3  =  6*self.E*self.Iz/(L*L)
        k4  = 12*self.E*self.Iy/(L*L*L)
        k5  =  6*self.E*self.Iy/(L*L)
        k6  =    self.GJ/L
        k7  =  4*self.E*self.Iy/L
        k8  =  2*self.E*self.Iy/L
        k9  =  4*self.E*self.Iz/L
        k10 =  2*self.E*self.Iz/L
        
        Kp = np.mat([[ k1,  0,  0,  0,  0,  0,-k1,  0,  0,  0,  0,  0],
                     [  0, k2,  0,  0,  0, k3,  0,-k2,  0,  0,  0, k3],
                     [  0,  0, k4,  0,-k5,  0,  0,  0,-k4,  0,-k5,  0],
                     [  0,  0,  0, k6,  0,  0,  0,  0,  0,-k6,  0,  0],
                     [  0,  0,-k5,  0, k7,  0,  0,  0, k5,  0, k8,  0],
                     [  0, k3,  0,  0,  0, k9,  0,-k3,  0,  0,  0,k10],
                     [-k1,  0,  0,  0,  0,  0, k1,  0,  0,  0,  0,  0],
                     [  0,-k2,  0,  0,  0,-k3,  0, k2,  0,  0,  0,-k3],
                     [  0,  0,-k4,  0, k5,  0,  0,  0, k4,  0, k5,  0],
                     [  0,  0,  0,-k6,  0,  0,  0,  0,  0, k6,  0,  0],
                     [  0,  0,-k5,  0, k8,  0,  0,  0, k5,  0, k7,  0],
                     [  0, k3,  0,  0,  0,k10,  0,-k3,  0,  0,  0, k9]]);

        if do_releases:
            if self.release_Ti:
                Kp[3,:] = 0;
                Kp[:,3] = 0;
            if self.release_Myi:
                Kp[4,:] = 0;
                Kp[:,4] = 0;
            if self.release_Mzi:
                Kp[5,:] = 0;
                Kp[:,5] = 0;
            if self.release_Tj:
                Kp[9,:] = 0;
                Kp[:,9] = 0; 
            if self.release_Myj:
                Kp[10,:] = 0;
                Kp[:,10] = 0;             
            if self.release_Mzj:
                Kp[11,:] = 0;
                Kp[:,11] = 0;            
 
        local_x  = np.array([(nJx-nIx)/L, (nJy-nIy)/L, (nJz-nIz)/L])
        local_xz = np.array([self.vec_xz[0], self.vec_xz[1], self.vec_xz[2]])
        local_y  = np.cross(local_x,local_xz)
        local_y  = local_y/np.linalg.norm(local_y)
        local_z  = np.cross(local_x,local_y)
        local_z  = local_z/np.linalg.norm(local_z)
        
        lx = local_x[0]
        mx = local_x[1]
        nx = local_x[2]
        
        ly = local_y[0]
        my = local_y[1]
        ny = local_y[2]

        lz = local_z[0]
        mz = local_z[1]
        nz = local_z[2]

        T  = np.mat([[lx,mx,nx, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [ly,my,ny, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [lz,mz,nz, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [ 0, 0, 0,lx,mx,nx, 0, 0, 0, 0, 0, 0],
                     [ 0, 0, 0,ly,my,ny, 0, 0, 0, 0, 0, 0],
                     [ 0, 0, 0,lz,mz,nz, 0, 0, 0, 0, 0, 0],       
                     [ 0, 0, 0, 0, 0, 0,lx,mx,nx, 0, 0, 0],
                     [ 0, 0, 0, 0, 0, 0,ly,my,ny, 0, 0, 0],
                     [ 0, 0, 0, 0, 0, 0,lz,mz,nz, 0, 0, 0],
                     [ 0, 0, 0, 0, 0, 0, 0, 0, 0,lx,mx,nx],
                     [ 0, 0, 0, 0, 0, 0, 0, 0, 0,ly,my,ny],
                     [ 0, 0, 0, 0, 0, 0, 0, 0, 0,lz,mz,nz]]);

        K = np.transpose(T)*Kp*T
        return K
        
    def disp(self,results):
        dofs = self.get_dof_ids()
        d = np.zeros(12)
        for i in range(len(dofs)):
            d[i] = results.d[dofs[i]] 
        return d
        
    def force(self,results):
        K_ele = self.StiffnessMatrix(True)
        d_ele = self.disp(results)
        return np.dot(K_ele,d_ele)
        
class PondingLoadCell2d:        
    def __init__(self,id,nodeI,nodeJ,gamma,tw):
        self.id = id
        self.nodeI = nodeI
        self.nodeJ = nodeJ
        self.gamma = gamma
        self.tw = tw
        self.gammas = 0
        self.hs = 0
    
    def get_dof_ids(self):
        dofs = [self.nodeI.dofs['UY'].id, 
                self.nodeJ.dofs['UY'].id]
        return dofs
        
    def get_load_vector(self,d,z):
        nIx = self.nodeI.coords[0]
        nIy = self.nodeI.coords[1]
        nJx = self.nodeJ.coords[0]
        nJy = self.nodeJ.coords[1]
        L  = nJx-nIx;
        hI = z - (nIy + d[self.nodeI.dofs['UY'].id])
        hJ = z - (nJy + d[self.nodeJ.dofs['UY'].id])
        
        if hI >= 0:
            if hJ >= 0:
                F  = -self.gamma*self.tw*(hI+hJ)*L/2
                x  = L*(2*hJ+hI)/(3*(hI+hJ))
            else:
                Lo = (hI)/(hI-hJ)*L
                F  = -self.gamma*self.tw*hI*Lo/2
                x  = Lo/3
        else:
            if hJ >= 0:
                Lo = (hJ)/(hJ-hI)*L
                F  = -self.gamma*self.tw*hJ*Lo/2
                x  = L-Lo/3
            else:
                F = 0
                x = 0.5*L
        
        if self.gammas > 0 and self.hs > 0:
            # Snow Force      
            Fs = -self.gammas*self.tw*self.hs*L
            xs = L/2
        
            # Overlap Adjustment Force
            gammaoa = min(self.gamma,self.gammas)
            Lxs = (self.hs-hI)*L/(hJ-hI) # length from I-end where the water level crosses the snow level
            Lxb = -hI*L/(hJ-hI)          # length from I-end where the water level crosses the beam 
            
            if hI >= self.hs:
                if hJ >= self.hs:
                    Foa = gammaoa*self.tw*self.hs*L
                    xoa = L/2                
                elif hJ >= 0:
                    F1 = gammaoa*self.tw*self.hs*Lxs
                    x1 = Lxs/2
                    F2 = gammaoa*self.tw*(self.hs+hJ)*(L-Lxs)/2
                    x2 = Lxs + (L-Lxs)*(2*hJ+self.hs)/(3*(hJ+self.hs))
                    Foa = F1 + F2
                    xoa = (x1*F1 + x2*F2)/Foa
                else:
                    F1 = gammaoa*self.tw*self.hs*Lxs
                    x1 = Lxs/2
                    F2 = gammaoa*self.tw*(self.hs)*(Lxb-Lxs)/2
                    x2 = Lxs + (Lxb-Lxs)/3
                    Foa = F1 + F2
                    xoa = (x1*F1 + x2*F2)/Foa
            elif hI >= 0:
                if hJ >= self.hs:
                    F1 = gammaoa*self.tw*(hI+self.hs)*(Lxs)/2
                    x1 = Lxs*(2*self.hs+hI)/(3*(self.hs+hI))
                    F2 = gammaoa*self.tw*self.hs*(L-Lxs)
                    x2 = Lxs + (L-Lxs)/2
                    Foa = F1 + F2
                    xoa = (x1*F1 + x2*F2)/Foa
                elif hJ >= 0:
                    Foa = gammaoa*self.tw*(hI+hJ)*L/2
                    xoa = L*(2*hJ+hI)/(3*(hJ+hI))
                else:
                    Foa = gammaoa*self.tw*hI*Lxb/2
                    xoa = Lxb/3
            else:
                if hJ >= self.hs:
                    F1 = gammaoa*self.tw*self.hs*(Lxs-Lxb)/2
                    x1 = Lxs - (Lxs-Lxb)/3
                    F2 = gammaoa*self.tw*self.hs*(L-Lxs)
                    x2 = L - (L-Lxs)/2
                    Foa = F1 + F2
                    xoa = (x1*F1 + x2*F2)/Foa
                elif hJ >= 0:
                    Foa = gammaoa*self.tw*hJ*(L-Lxb)/2
                    xoa = L - (L-Lxb)/3
                else:
                    Foa = 0
                    xoa = 0

            # Total Force
            x = (x*F + xs*Fs + xoa*Foa)/(F+Fs+Foa)
            F = F + Fs + Foa
            
        f = np.mat([[(1-x/L)*F],
                    [  (x/L)*F]])
        return f
        
    def get_volume(self,d,z):       
        nIx = self.nodeI.coords[0]
        nIy = self.nodeI.coords[1]
        nJx = self.nodeJ.coords[0]
        nJy = self.nodeJ.coords[1]
        L  = nJx-nIx;
        hI = z - (nIy + d[self.nodeI.dofs['UY'].id])
        hJ = z - (nJy + d[self.nodeJ.dofs['UY'].id])
        if hI >= 0:
            if hJ >= 0:
                V    = self.tw*(hI+hJ)*L/2
                dVdz = self.tw*L
            else:
                Lo   = (hI)/(hI-hJ)*L
                V    = self.tw*hI*Lo/2
                dVdz = self.tw*Lo
        else:
            if hJ >= 0:
                Lo   = (hJ)/(hJ-hI)*L
                V    = self.tw*hJ*Lo/2
                dVdz = self.tw*Lo
            else:
                V    = 0
                dVdz = 0
        if self.gammas > 0 and self.hs > 0:
            raise Exception('get_volume not yet implemented for cases with snow')
        return (V,dVdz)
        
class PondingLoadCell3d:        
    def __init__(self,id,nodeI,nodeJ,nodeK,nodeL,gamma,na=1,nb=1):
        self.id = id
        self.nodeI = nodeI
        self.nodeJ = nodeJ
        self.nodeK = nodeK
        self.nodeL = nodeL
        self.gamma = gamma  # Fluid density
        self.na    = na     # Number of sub-cells along IJ
        self.nb    = nb     # Number of sub-cells along JK
        self.gammas = 0
        self.hs = 0        
        
    def get_dof_ids(self):
        dofs = [self.nodeI.dofs['UZ'].id, 
                self.nodeJ.dofs['UZ'].id, 
                self.nodeK.dofs['UZ'].id, 
                self.nodeL.dofs['UZ'].id]
        return dofs
        
    def get_load_vector(self,d,z):
        xI = self.nodeI.coords[0]
        yI = self.nodeI.coords[1]
        zI = self.nodeI.coords[2]
        xJ = self.nodeJ.coords[0]
        yJ = self.nodeJ.coords[1]
        zJ = self.nodeJ.coords[2]
        xK = self.nodeK.coords[0]
        yK = self.nodeK.coords[1]
        zK = self.nodeK.coords[2]
        xL = self.nodeL.coords[0]
        yL = self.nodeL.coords[1]
        zL = self.nodeL.coords[2]

        coords = np.mat([[xI,yI],
                         [xJ,yJ],
                         [xK,yK],
                         [xL,yL]])
        
        hI = z - (zI + d[self.nodeI.dofs['UZ'].id])
        hJ = z - (zJ + d[self.nodeJ.dofs['UZ'].id])
        hK = z - (zK + d[self.nodeK.dofs['UZ'].id])
        hL = z - (zL + d[self.nodeL.dofs['UZ'].id])
        
        # Define numerical integration points and weights
        n_ip   = 4
        xi_ip  = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3),-1/sqrt(3)] 
        eta_ip = [-1/sqrt(3),-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)] 
        w_ip   = [ 1, 1, 1, 1]        

        # Calculate load
        f = np.zeros((4,1))        
        if self.na == 1 and self.nb == 1:
            h = np.array([[max(0,hI)],[max(0,hJ)],[max(0,hK)],[max(0,hL)]])
            
            for iip in range(n_ip):
                j = self.Jacobian(xi_ip[iip],eta_ip[iip],coords)
                N = self.ShapeFunction(xi_ip[iip],eta_ip[iip])                
                f += self.gamma*j*N.dot(np.transpose(N).dot(h))
        else:
            h = np.array([[hI],[hJ],[hK],[hL]])
            
            for ia in range(self.na):
                for ib in range(self.nb):
                    #print('Sub-cell a = %i b = %i \n' % (ia,ib))
                
                    xi_sub  = [-1+2*ia/self.na,-1+2*(ia+1)/self.na,-1+2*(ia+1)/self.na,-1+2*ia/self.na]
                    eta_sub = [-1+2*ib/self.nb,-1+2*ib/self.nb,-1+2*(ib+1)/self.nb,-1+2*(ib+1)/self.nb]
                    
                    #print(xi_sub)
                    #print(eta_sub)
                    
                    # Compute coordinates and height of ponded water in the sub-cell
                    coords_sub = np.zeros((4,2))
                    h_sub = np.zeros((4,1)) 
                    for i in range(4):
                        N = self.ShapeFunction(xi_sub[i],eta_sub[i])
                        coords_sub[i,:] = np.transpose(N).dot(coords)
                        h_sub[i] = max(0,np.transpose(N).dot(h))
                    
                    #print(coords_sub)
                    #print(h_sub)
                    
                    # Compute sub-cell forces
                    f_sub = np.zeros((4,1)) 
                    for iip in range(n_ip):
                        j = self.Jacobian(xi_ip[iip],eta_ip[iip],coords_sub)
                        N = self.ShapeFunction(xi_ip[iip],eta_ip[iip])                
                        f_sub += self.gamma*j*N.dot(np.transpose(N).dot(h_sub))
                    
                    # Convert sub-cell to cell forces
                    for i in range(4):
                        N = self.ShapeFunction(xi_sub[i],eta_sub[i])
                        f = f + N*f_sub[i]
        return f        
        
        
    def get_volume(self,d,z):       
        xI = self.nodeI.coords[0]
        yI = self.nodeI.coords[1]
        zI = self.nodeI.coords[2]
        xJ = self.nodeJ.coords[0]
        yJ = self.nodeJ.coords[1]
        zJ = self.nodeJ.coords[2]
        xK = self.nodeK.coords[0]
        yK = self.nodeK.coords[1]
        zK = self.nodeK.coords[2]
        xL = self.nodeL.coords[0]
        yL = self.nodeL.coords[1]
        zL = self.nodeL.coords[2]
        
        coords = np.mat([[xI,yI],
                         [xJ,yJ],
                         [xK,yK],
                         [xL,yL]])        
        
        hI = z - (zI + d[self.nodeI.dofs['UZ'].id])
        hJ = z - (zJ + d[self.nodeJ.dofs['UZ'].id])
        hK = z - (zK + d[self.nodeK.dofs['UZ'].id])
        hL = z - (zL + d[self.nodeL.dofs['UZ'].id])
        h = np.array([[hI],[hJ],[hK],[hL]])
        
        V = self.get_load_vector(d,z).sum()/self.gamma
        dVdz = 0
        
        for ia in range(self.na):
            for ib in range(self.nb):        
                xi_sub  = [-1+2*ia/self.na,-1+2*(ia+1)/self.na,-1+2*(ia+1)/self.na,-1+2*ia/self.na]
                eta_sub = [-1+2*ib/self.nb,-1+2*ib/self.nb,-1+2*(ib+1)/self.nb,-1+2*(ib+1)/self.nb]
                
                # Compute coordinates and height of ponded water in the sub-cell
                coords_sub = np.zeros((4,2))
                h_sub = np.zeros((4,1)) 
                for i in range(4):
                    N = self.ShapeFunction(xi_sub[i],eta_sub[i])
                    coords_sub[i,:] = np.transpose(N).dot(coords)
                    h_sub[i] = np.transpose(N).dot(h)
                                               
                # Determine pologon where h > 0
                x_coord = np.empty(0)
                y_coord = np.empty(0)
                if h_sub[0] >= 0:
                    x_coord = np.append(x_coord,coords_sub[0,0])
                    y_coord = np.append(y_coord,coords_sub[0,1])
                if (h_sub[0] > 0 and h_sub[1] < 0) or (h_sub[0] < 0 and h_sub[1] > 0):
                    a = h_sub[0]/(h_sub[0]-h_sub[1])
                    x_coord = np.append(x_coord,coords_sub[0,0]+a*(coords_sub[1,0]-coords_sub[0,0]))
                    y_coord = np.append(y_coord,coords_sub[0,1]+a*(coords_sub[1,1]-coords_sub[0,1]))
                if h_sub[1] >= 0:
                    x_coord = np.append(x_coord,coords_sub[1,0])
                    y_coord = np.append(y_coord,coords_sub[1,1])
                if (h_sub[1] > 0 and h_sub[2] < 0) or (h_sub[1] < 0 and h_sub[2] > 0):
                    a = h_sub[1]/(h_sub[1]-h_sub[2])
                    x_coord = np.append(x_coord,coords_sub[1,0]+a*(coords_sub[2,0]-coords_sub[1,0]))
                    y_coord = np.append(y_coord,coords_sub[1,1]+a*(coords_sub[2,1]-coords_sub[1,1]))                 
                if h_sub[2] >= 0:
                    x_coord = np.append(x_coord,coords_sub[2,0])
                    y_coord = np.append(y_coord,coords_sub[2,1])
                if (h_sub[2] > 0 and h_sub[3] < 0) or (h_sub[2] < 0 and h_sub[3] > 0):
                    a = h_sub[2]/(h_sub[2]-h_sub[3])
                    x_coord = np.append(x_coord,coords_sub[2,0]+a*(coords_sub[3,0]-coords_sub[2,0]))
                    y_coord = np.append(y_coord,coords_sub[2,1]+a*(coords_sub[3,1]-coords_sub[2,1]))
                if h_sub[3] >= 0:
                    x_coord = np.append(x_coord,coords_sub[3,0])
                    y_coord = np.append(y_coord,coords_sub[3,1])
                if (h_sub[3] > 0 and h_sub[0] < 0) or (h_sub[3] < 0 and h_sub[0] > 0):
                    a = h_sub[3]/(h_sub[3]-h_sub[0])
                    x_coord = np.append(x_coord,coords_sub[3,0]+a*(coords_sub[0,0]-coords_sub[3,0]))
                    y_coord = np.append(y_coord,coords_sub[3,1]+a*(coords_sub[0,1]-coords_sub[3,1]))               
                
                # Compute area of polygon and add it to dVdz
                if x_coord.size > 0:
                    dVdz += 0.5*np.abs(np.dot(x_coord,np.roll(y_coord,1))-np.dot(y_coord,np.roll(x_coord,1)))
                
        return (V,dVdz)        

    @staticmethod
    def ShapeFunction(xi,eta):
        N = np.array([[(1-xi)*(1-eta)],
                      [(1+xi)*(1-eta)],
                      [(1+xi)*(1+eta)],
                      [(1-xi)*(1+eta)]])/4
        return N
        
    @staticmethod    
    def Jacobian(xi,eta,coords):
        dNd_ = np.array([[-(1-eta), (1-eta), (1+eta),-(1+eta)],
                         [ -(1-xi), -(1+xi),  (1+xi),  (1-xi)]])/4
        jac = np.dot(dNd_,coords)
        j   = np.linalg.det(jac)
        return j
        
class LinearAnalysis:

    def __init__(self,model):
        self.model = model;
        
    def run(self,load_factors):       
        K = self.model.GetGlobalStiffnessMatrix()
        f = self.model.GetNodalForceVector(load_factors)
        (d,r) = self.model.SolveForDisp(K,f)
        
        # store results
        self.load_factors = load_factors;
        self.d = d
        self.r = r


class PondingAnalysis:
    max_iterations_z = 20
    max_iter_const_V = 20
    max_iter_find_z  = 10
    tol_z = 1e-3
    tol_V = 1e-4
    output_level = 0
    use_stored_analysis = False

    def __init__(self,model,type='Constant_Level'):
        self.model = model 
        self.type = type 
        
    def run(self,load_factors,x):
        if self.type == 'Constant_Level':
            z = x
            
            f_nodal = self.model.GetNodalForceVector(load_factors)
            
            # Run Initial Analysis
            if self.use_stored_analysis:
                (d,r) = self.model.SolveForDispWithStored(f_nodal)
            else:
                K = self.model.GetGlobalStiffnessMatrix()
                (d,r) = self.model.SolveForDisp(K,f_nodal)
                
            # Iterate
            for i in range(self.max_iterations_z):
                d_last = d
                f = f_nodal + self.model.GetPondingForceVector(d,z)
                
                if self.use_stored_analysis:
                    (d,r) = self.model.SolveForDispWithStored(f)
                else:
                    (d,r) = self.model.SolveForDisp(K,f)
                
                if self.output_level > 0:
                    print('Min Deflection: %.3f \tNode Disp Incr. %.6f' % (min(d),np.linalg.norm(d-d_last)))
                if np.linalg.norm(d-d_last) < self.tol_z:
                    if self.output_level > 0:
                        print('Converged')
                    # store results
                    self.load_factors = load_factors
                    self.z = z
                    self.d = d
                    self.r = r
                    return 0
                    
        elif self.type == 'Constant_Volume':
            V = x
            
            f_nodal = self.model.GetNodalForceVector(load_factors)
            
            # Run Initial Analysis
            if self.use_stored_analysis:
                (d,r) = self.model.SolveForDispWithStored(f_nodal)
            else:
                K = self.model.GetGlobalStiffnessMatrix()
                (d,r) = self.model.SolveForDisp(K,f_nodal)            
            
            # Iterate
            z = 1
            for i in range(self.max_iter_const_V):
                d_last = d
                # determine z given V
                for j in range(self.max_iter_find_z):
                    z_last = z
                    res = self.model.GetPondingVolume(d,z)
                    V_calc = res[0]
                    dVdz   = res[1]
                    if self.output_level > 0:
                        print('Iteration = %i' % j)
                        print('V         = %g' % V)
                        print('V_calc    = %g' % V_calc)
                        print('z         = %g' % (z + (V-V_calc)/dVdz))
                        print('z_last    = %g' % z_last)
                    if abs(V-V_calc)/V < self.tol_V:
                        if self.output_level > 0:
                            print('found z')
                        break
                    z = z + (V-V_calc)/dVdz
                    if j == self.max_iter_find_z-1:
                        if self.output_level > 0:
                            print('Could not find z')
                        return -1
                
                # solve for d given z
                f = f_nodal + self.model.GetPondingForceVector(d,z)
                
                if self.use_stored_analysis:
                    (d,r) = self.model.SolveForDispWithStored(f)
                else:
                    (d,r) = self.model.SolveForDisp(K,f)  
                
                if self.output_level > 0:
                    print('Min Deflection: %.3f \tNode Disp Incr. %.6f' % (min(d),np.linalg.norm(d-d_last)))
                if np.linalg.norm(d-d_last) < self.tol_z:
                    if self.output_level > 0:
                        print('Converged')
                    # store results
                    self.load_factors = load_factors
                    self.z = z
                    self.d = d
                    self.r = r
                    return 0          
                    
                if i == self.max_iter_const_V-1:
                    if self.output_level > 0:
                        print('Could not find a solution')
                    return -1
        
        elif self.type == 'No_Ponding_Effect':
            z = x
            
            f = self.model.GetNodalForceVector(load_factors) + self.model.GetPondingForceVector(np.zeros(self.model.ndof),z)
            
            if self.use_stored_analysis:
                (d,r) = self.model.SolveForDispWithStored(f)
            else:
                K = self.model.GetGlobalStiffnessMatrix()
                (d,r) = self.model.SolveForDisp(K,f)            

            # store results
            self.load_factors = load_factors
            self.z = z
            self.d = d
            self.r = r
            