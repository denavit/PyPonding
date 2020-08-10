import numpy as np
from math import sqrt

class PondingLoadCell2d:
    id      = ''    # Load cell ID
    
    xI      = 0.0   # X coordinate of Node I
    yI      = 0.0   # Y coordinate of Node I
    xJ      = 0.0   # X coordinate of Node J
    yJ      = 0.0   # Y coordinate of Node J
    
    dyI     = 0.0   # Y deflection of Node I
    dyJ     = 0.0   # Y deflection of Node J
    
    gamma   = 0     # Fluid density
    tw      = 0     # Tributary width
    
    gammas  = 0     # Snow density
    hs      = 0     # Snow height

    def __init__(self):
        pass
           
    def get_load_vector(self,z):
        L  = self.xJ-self.xI;
        
        hI = z - (self.yI + self.dyI)
        hJ = z - (self.yJ + self.dyJ)
        
        if hI >= 0:
            if hJ == 0:
                F = 0
                x = 0.5*L
            elif hJ > 0:
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
        
    def get_volume(self,z):       
        L  = self.xJ-self.xI;
        
        hI = z - (self.yI + self.dyI)
        hJ = z - (self.yJ + self.dyJ)

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
    id      = ''    # Load cell ID
    
    # Define nodes (vertices) in counterclockwise (CCW) direction
    xI      = 0.0   # X coordinate of Node I
    yI      = 0.0   # Y coordinate of Node I
    zI      = 0.0   # Z coordinate of Node I
    xJ      = 0.0   # X coordinate of Node J
    yJ      = 0.0   # Y coordinate of Node J
    zJ      = 0.0   # Z coordinate of Node J
    xK      = 0.0   # X coordinate of Node K
    yK      = 0.0   # Y coordinate of Node K
    zK      = 0.0   # Z coordinate of Node K
    xL      = 0.0   # X coordinate of Node L
    yL      = 0.0   # Y coordinate of Node L
    zL      = 0.0   # Z coordinate of Node L

    dzI     = 0.0   # Z deflection of Node I
    dzJ     = 0.0   # Z deflection of Node J
    dzK     = 0.0   # Z deflection of Node K
    dzL     = 0.0   # Z deflection of Node L
    
    gamma   = 0     # Fluid density
    na      = 1     # Number of sub-cells along IJ
    nb      = 1     # Number of sub-cells along JK
    
    gammas  = 0     # Snow density
    hs      = 0     # Snow height
    
    def __init__(self):
        pass
        
    def get_load_vector(self,z):
        coords = np.mat([[self.xI,self.yI],
                         [self.xJ,self.yJ],
                         [self.xK,self.yK],
                         [self.xL,self.yL]])
        
        hI = z - (self.zI + self.dzI)
        hJ = z - (self.zJ + self.dzJ)
        hK = z - (self.zK + self.dzK)
        hL = z - (self.zL + self.dzL)
        
        # Define numerical integration points and weights
        n_ip   = 4
        xi_ip  = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3),-1/sqrt(3)] 
        eta_ip = [-1/sqrt(3),-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)] 
        w_ip   = [ 1, 1, 1, 1]        

        # Calculate load
        f = np.zeros((4,1))        
        if self.na == 1 and self.nb == 1:
            
            # Compute pressure due to water and snow at each corner of the cell
            if self.gammas > 0 and self.hs > 0:
                if hI <= 0:
                    wpI = self.gammas*self.hs
                elif hI <= self.hs:
                    wpI = max(self.gamma,self.gammas)*hI + self.gammas*(self.hs-hI)
                else:
                    wpI = max(self.gamma,self.gammas)*self.hs + self.gamma*(hI-self.hs)

                if hJ <= 0:
                    wpJ = self.gammas*self.hs
                elif hJ <= self.hs:
                    wpJ = max(self.gamma,self.gammas)*hJ + self.gammas*(self.hs-hJ)
                else:
                    wpJ = max(self.gamma,self.gammas)*self.hs + self.gamma*(hJ-self.hs)                    
                    
                if hK <= 0:
                    wpK = self.gammas*self.hs
                elif hK <= self.hs:
                    wpK = max(self.gamma,self.gammas)*hK + self.gammas*(self.hs-hK)
                else:
                    wpK = max(self.gamma,self.gammas)*self.hs + self.gamma*(hK-self.hs)

                if hL <= 0:
                    wpL = self.gammas*self.hs
                elif hL <= self.hs:
                    wpL = max(self.gamma,self.gammas)*hL + self.gammas*(self.hs-hL)
                else:
                    wpL = max(self.gamma,self.gammas)*self.hs + self.gamma*(hL-self.hs)                      
                    
                wp = np.array([[wpI],[wpJ],[wpK],[wpL]])
            else:
                wp = self.gamma*np.array([[max(0,hI)],[max(0,hJ)],[max(0,hK)],[max(0,hL)]])
            
            # Compute the force vector
            for iip in range(n_ip):
                j = self.Jacobian(xi_ip[iip],eta_ip[iip],coords)
                N = self.ShapeFunction(xi_ip[iip],eta_ip[iip])                
                f += j*N.dot(np.transpose(N).dot(-wp))
                
        else:
            h = np.array([[hI],[hJ],[hK],[hL]])
            
            # Loop over each sub-cell
            for ia in range(self.na):
                for ib in range(self.nb):
                
                    # Define coordinates (in local coordinates) of the corners of the sub-cell
                    xi_sub  = [-1+2*ia/self.na,-1+2*(ia+1)/self.na,-1+2*(ia+1)/self.na,-1+2*ia/self.na]
                    eta_sub = [-1+2*ib/self.nb,-1+2*ib/self.nb,-1+2*(ib+1)/self.nb,-1+2*(ib+1)/self.nb]
                    
                    # Compute for each corner of the sub-cell...
                    coords_sub = np.zeros((4,2))
                    wp_sub = np.zeros((4,1)) 
                    for i in range(4):
                        N = self.ShapeFunction(xi_sub[i],eta_sub[i])
                        
                        # Coordinates (in global coordinates)
                        coords_sub[i,:] = np.transpose(N).dot(coords)
                        
                        # Height of water at corner of sub-cell
                        h_sub = np.transpose(N).dot(h)
                    
                        # Pressure due to water and snow
                        if self.gammas > 0 and self.hs > 0:
                            if h_sub <= 0:
                                wp_sub[i] = self.gammas*self.hs
                            elif h_sub <= self.hs:
                                wp_sub[i] = max(self.gamma,self.gammas)*h_sub + self.gammas*(self.hs-h_sub)
                            else:
                                wp_sub[i] = max(self.gamma,self.gammas)*self.hs + self.gamma*(h_sub-self.hs)
                        else:
                            wp_sub[i] = self.gamma*max(0,h_sub)
                    
                    # Compute sub-cell force vector
                    f_sub = np.zeros((4,1)) 
                    for iip in range(n_ip):
                        j = self.Jacobian(xi_ip[iip],eta_ip[iip],coords_sub)
                        N = self.ShapeFunction(xi_ip[iip],eta_ip[iip])                
                        f_sub += j*N.dot(np.transpose(N).dot(-wp_sub))
                    
                    # Convert sub-cell force vector to cell force vector
                    for i in range(4):
                        N = self.ShapeFunction(xi_sub[i],eta_sub[i])
                        f = f + N*f_sub[i]
        return f        
        
        
    def get_volume(self,z):       
        coords = np.mat([[self.xI,self.yI],
                         [self.xJ,self.yJ],
                         [self.xK,self.yK],
                         [self.xL,self.yL]])
        
        hI = z - (self.zI + self.dzI)
        hJ = z - (self.zJ + self.dzJ)
        hK = z - (self.zK + self.dzK)
        hL = z - (self.zL + self.dzL)
        
        h = np.array([[hI],[hJ],[hK],[hL]])
        
        V = -self.get_load_vector(z).sum()/self.gamma
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
