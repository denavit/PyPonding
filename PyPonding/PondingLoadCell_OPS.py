from math import sin,cos,atan2,sqrt
from PyPonding.PondingLoadCell import PondingLoadCell2d, PondingLoadCell3d

try:
    import opensees as ops
except ImportError:
    try:
        import openseespy.opensees as ops
    except ImportError:
        import warnings
        warnings.warn('OpenSeesPy not found. OpenSees functionality is not available.')


class NodeEnd2d:
    y_offset = 0.

    def __init__(self,node_id):
        self.node_id = node_id
       
    def key(self):
        return self.node_id
        
    def coord(self):
        x = ops.nodeCoord(self.node_id,1)
        y = ops.nodeCoord(self.node_id,2) + self.y_offset
        return (x,y)
    
    def disp(self):
        dx = ops.nodeDisp(self.node_id,1)
        dy = ops.nodeDisp(self.node_id,2)
        return (dx,dy)
    
class ElementEnd2d:
    y_offset = 0.
    
    def __init__(self,element_id,x):
        self.element_id = element_id
        self.nodes = ops.eleNodes(element_id)
        self.x = x
 
    def key(self):
        return (self.element_id,round(self.x,6))
 
    def coord(self):
        xI = ops.nodeCoord(self.nodes[0],1)
        yI = ops.nodeCoord(self.nodes[0],2)
        xJ = ops.nodeCoord(self.nodes[1],1)
        yJ = ops.nodeCoord(self.nodes[1],2)
        x = xI + self.x*(xJ-xI)
        y = yI + self.x*(yJ-yI) + self.y_offset
        return (x,y)
    
    def disp(self):
        dx,dy,dz = ops.cbdiDisplacement(self.element_id,self.x)
        return (dx,dy)

def define_end(arg):
    if isinstance(arg, int):
        return NodeEnd2d(arg)
    elif isinstance(arg, float):
        if arg.is_integer():
            return NodeEnd2d(arg)
        else:
            raise Exception('PondingLoadCell2d_OPS end definition not valid. If numertic, input needs to be an integer.')    
    elif isinstance(arg, tuple):
        if isinstance(arg[0], str):
            if arg[0].lower() == 'node':
                return NodeEnd2d(arg[1])
            elif arg[0].lower() == 'element':
                return ElementEnd2d(arg[1],arg[2])
            else:
                raise Exception('PondingLoadCell2d_OPS end definition not valid. Unknown string: %s' % arg[0]) 
        else:
            raise Exception('PondingLoadCell2d_OPS end definition not valid. If tuple, first item needs to be a string.')  
    else:
        raise Exception('PondingLoadCell2d_OPS end definition not valid. Unknown argument type: %s' % type(arg))  
       
class PondingLoadCell2d_OPS(PondingLoadCell2d):
    def __init__(self,id,endI,endJ,gamma,tw):
        self.id = id
        self.endI = define_end(endI)
        self.endJ = define_end(endJ)
        self.gamma = gamma
        self.tw = tw
        
        # Retreive Coordinates
        self.update_coord()
        
        # Store node ids (if attached to nodes) for backwards compatibility 
        if isinstance(self.endI, NodeEnd2d):
            self.nodeI = self.endI.node_id
        else:
            self.nodeI = None
        if isinstance(self.endJ, NodeEnd2d):
            self.nodeJ = self.endJ.node_id
        else:
            self.nodeJ = None        
        
    def update_coord(self):
        self.xI, self.yI = self.endI.coord()
        self.xJ, self.yJ = self.endJ.coord()
        
    def update(self):
        # Code currently only updates y postion of nodes - @todo maybe update x position also
        dx,dy = self.endI.disp()
        # self.dxI = dx
        self.dyI = dy
        dx,dy = self.endJ.disp()
        # self.dxJ = dx
        self.dyJ = dy

class PondingLoadManager2d:
    def __init__(self):
        self.cells = dict()
    
    def add_cell(self,id,endI,endJ,gamma,tw):
        self.cells[id] = PondingLoadCell2d_OPS(id,endI,endJ,gamma,tw)
        self.build_load_vector()
        
    def update(self):
        for i in self.cells:
            self.cells[i].update()
    
    def get_volume(self,zw):
        V = 0
        dVdz = 0
        for i in self.cells:
            (iV,idVdz) = self.cells[i].get_volume(zw)
            V += iV
            dVdz += idVdz
        return (V,dVdz)
        
    def build_load_vector(self):
        nodal_loads = dict()
        element_loads = dict()
        for i in self.cells:
            if isinstance(self.cells[i].endI , NodeEnd2d):
                if not self.cells[i].endI.key() in nodal_loads:
                    nodal_loads[self.cells[i].endI.key()] = 0.0    
            elif isinstance(self.cells[i].endI , ElementEnd2d):
                if not self.cells[i].endI.key() in element_loads:
                    element_loads[self.cells[i].endI.key()] = 0.0
            else:
                raise Exception('Unknown endI: %s' % type(self.cells[i].endI))
            
            if isinstance(self.cells[i].endJ , NodeEnd2d):
                if not self.cells[i].endJ.key() in nodal_loads:
                    nodal_loads[self.cells[i].endJ.key()] = 0.0    
            elif isinstance(self.cells[i].endJ , ElementEnd2d):
                if not self.cells[i].endJ.key() in element_loads:
                    element_loads[self.cells[i].endJ.key()] = 0.0
            else:
                raise Exception('Unknown endJ: %s' % type(self.cells[i].endJ))            
                    
        self.empty_nodal_loads = nodal_loads;
        self.empty_element_loads = element_loads;
        self.committed_nodal_loads = nodal_loads.copy();
        self.committed_element_loads = element_loads.copy();
        
    def compute_current_load_vector(self,zw):
        nodal_loads = self.empty_nodal_loads.copy()
        element_loads = self.empty_element_loads.copy()
        
        for i in self.cells:  
            f = self.cells[i].get_load_vector(zw)
            
            if isinstance(self.cells[i].endI , NodeEnd2d):
                nodal_loads[self.cells[i].endI.key()] += f.item(0)
            elif isinstance(self.cells[i].endI , ElementEnd2d):
                element_loads[self.cells[i].endI.key()] += f.item(0)
            else:
                raise Exception('Unknown endI: %s' % type(self.cells[i].endI))
            
            if isinstance(self.cells[i].endJ , NodeEnd2d):
                nodal_loads[self.cells[i].endJ.key()] += f.item(1)
            elif isinstance(self.cells[i].endJ , ElementEnd2d):
                element_loads[self.cells[i].endJ.key()] += f.item(1)
            else:
                raise Exception('Unknown endJ: %s' % type(self.cells[i].endJ))
                
        self.current_nodal_loads = nodal_loads
        self.current_element_loads = element_loads
        
    def commit_current_load_vector(self):
        self.committed_nodal_loads = self.current_nodal_loads.copy();
        self.committed_element_loads = self.current_element_loads.copy();

    def apply_load_increment(self):
        for i in self.committed_nodal_loads:        
            fy = self.current_nodal_loads[i] - self.committed_nodal_loads[i]
            ops.load(i, 0.0, fy, 0.0)
        for i in self.committed_element_loads:
            fy = self.current_element_loads[i] - self.committed_element_loads[i]
            iNodes = ops.eleNodes(i[0])
            xI,yI = ops.nodeCoord(iNodes[0])
            xJ,yJ = ops.nodeCoord(iNodes[1])
            ele_angle = atan2(yJ-yI, xJ-xI)
            ops.eleLoad('-ele', i[0], '-type', '-beamPoint', fy*cos(ele_angle), i[1], fy*sin(ele_angle))


class NodeVertex3d:
    z_offset = 0.
    
    def __init__(self,node_id):
        self.node_id = node_id
       
    def key(self):
        return self.node_id
        
    def coord(self):
        x = ops.nodeCoord(self.node_id,1)
        y = ops.nodeCoord(self.node_id,2)
        z = ops.nodeCoord(self.node_id,3) + self.z_offset
        return (x,y,z)
    
    def disp(self):
        dx = ops.nodeDisp(self.node_id,1)
        dy = ops.nodeDisp(self.node_id,2)
        dz = ops.nodeDisp(self.node_id,3)
        return (dx,dy,dz)
    
class ElementVertex3d:
    z_offset = 0.

    def __init__(self,element_id,x):
        self.element_id = element_id
        self.nodes = ops.eleNodes(element_id)
        self.x = x
 
    def key(self):
        return (self.element_id,round(self.x,6))
 
    def coord(self):
        xI = ops.nodeCoord(self.nodes[0],1)
        yI = ops.nodeCoord(self.nodes[0],2)
        zI = ops.nodeCoord(self.nodes[0],3)
        xJ = ops.nodeCoord(self.nodes[1],1)
        yJ = ops.nodeCoord(self.nodes[1],2)
        zJ = ops.nodeCoord(self.nodes[1],3)
        x = xI + self.x*(xJ-xI)
        y = yI + self.x*(yJ-yI)
        z = zI + self.x*(zJ-zI) + self.z_offset
        return (x,y,z)
    
    def disp(self):
        dx,dy,dz = ops.cbdiDisplacement(self.element_id,self.x)
        return (dx,dy,dz)

class FixedVertex3d:
    z_offset = 0.
    
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
       
    def key(self):
        return (round(self.x,6),round(self.y,6),round(self.z,6))
        
    def coord(self):
        return (self.x, self.y, self.z + self.z_offset)
    
    def disp(self):
        return (0.,0.,0.)

def define_vertex(arg):
    if isinstance(arg, int):
        return NodeVertex3d(arg)
    elif isinstance(arg, float):
        if arg.is_integer():
            return NodeVertex3d(arg)
        else:
            raise Exception('PondingLoadCell3d_OPS end definition not valid. If numertic, input needs to be an integer.')    
    elif isinstance(arg, tuple):
        if isinstance(arg[0], str):
            if arg[0].lower() == 'node':
                return NodeVertex3d(arg[1])
            elif arg[0].lower() == 'element':
                return ElementVertex3d(arg[1],arg[2])
            elif arg[0].lower() == 'fixed':
                return FixedVertex3d(arg[1],arg[2],arg[3])
            else:
                raise Exception('PondingLoadCell3d_OPS end definition not valid. Unknown string: %s' % arg[0]) 
        else:
            raise Exception('PondingLoadCell3d_OPS end definition not valid. If tuple, first item needs to be a string.')  
    else:
        raise Exception('PondingLoadCell3d_OPS end definition not valid. Unknown argument type: %s' % type(arg))  
            

class PondingLoadCell3d_OPS(PondingLoadCell3d):
    # Define vertices in counterclockwise (CCW) direction

    def __init__(self, id, vertexI, vertexJ, vertexK, vertexL, gamma, na=1, nb=1):
        self.id = id
        self.vertexI = define_vertex(vertexI)
        self.vertexJ = define_vertex(vertexJ)
        self.vertexK = define_vertex(vertexK)
        self.vertexL = define_vertex(vertexL)
        self.gamma = gamma
        self.na = na
        self.nb = nb

        # Retreive Coordinates
        self.update_coord()
                
        # Store node ids (if attached to nodes) for backwards compatibility 
        if isinstance(self.vertexI, NodeVertex3d):
            self.nodeI = self.vertexI.node_id
        else:
            self.nodeI = None
        if isinstance(self.vertexJ, NodeVertex3d):
            self.nodeJ = self.vertexJ.node_id
        else:
            self.nodeJ = None
        if isinstance(self.vertexK, NodeVertex3d):
            self.nodeK = self.vertexK.node_id
        else:
            self.nodeK = None
        if isinstance(self.vertexL, NodeVertex3d):
            self.nodeL = self.vertexL.node_id
        else:
            self.nodeL = None            
        
    def update_coord(self):
        self.xI, self.yI, self.zI = self.vertexI.coord()
        self.xJ, self.yJ, self.zJ = self.vertexJ.coord()
        self.xK, self.yK, self.zK = self.vertexK.coord()
        self.xL, self.yL, self.zL = self.vertexL.coord()    
    
    def update(self):
        # Code currently only updates z postion of nodes - @todo maybe update x and y positions also
        self.dzI = self.vertexI.disp()[2]
        self.dzJ = self.vertexJ.disp()[2]
        self.dzK = self.vertexK.disp()[2]
        self.dzL = self.vertexL.disp()[2]
    
class PondingLoadManager3d:
    def __init__(self):
        self.cells = dict()
       
    def add_cell(self,id, vertexI, vertexJ, vertexK, vertexL, gamma, na=1, nb=1):
        self.cells[id] = PondingLoadCell3d_OPS(id,vertexI,vertexJ,vertexK,vertexL,gamma,na,nb)
        self.build_load_vector()
       
    def update(self):
        for i in self.cells:
            self.cells[i].update()
    
    def get_volume(self,zw):
        V = 0
        dVdz = 0
        for i in self.cells:
            (iV,idVdz) = self.cells[i].get_volume(zw)
            V += iV
            dVdz += idVdz
        return (V,dVdz)
        
    def build_load_vector(self):
        nodal_loads = dict()
        element_loads = dict()
        for i in self.cells:
            if isinstance(self.cells[i].vertexI , NodeVertex3d):
                if not self.cells[i].vertexI.key() in nodal_loads:
                    nodal_loads[self.cells[i].vertexI.key()] = 0.0    
            elif isinstance(self.cells[i].vertexI , ElementVertex3d):
                if not self.cells[i].vertexI.key() in element_loads:
                    element_loads[self.cells[i].vertexI.key()] = 0.0
            elif isinstance(self.cells[i].vertexI , FixedVertex3d):
                pass
            else:
                raise Exception('Unknown type for vertexI: %s' % type(self.cells[i].vertexI))

            if isinstance(self.cells[i].vertexJ , NodeVertex3d):
                if not self.cells[i].vertexJ.key() in nodal_loads:
                    nodal_loads[self.cells[i].vertexJ.key()] = 0.0    
            elif isinstance(self.cells[i].vertexJ , ElementVertex3d):
                if not self.cells[i].vertexJ.key() in element_loads:
                    element_loads[self.cells[i].vertexJ.key()] = 0.0
            elif isinstance(self.cells[i].vertexJ , FixedVertex3d):
                pass
            else:
                raise Exception('Unknown type for vertexJ: %s' % type(self.cells[i].vertexJ))

            if isinstance(self.cells[i].vertexK , NodeVertex3d):
                if not self.cells[i].vertexK.key() in nodal_loads:
                    nodal_loads[self.cells[i].vertexK.key()] = 0.0    
            elif isinstance(self.cells[i].vertexK , ElementVertex3d):
                if not self.cells[i].vertexK.key() in element_loads:
                    element_loads[self.cells[i].vertexK.key()] = 0.0
            elif isinstance(self.cells[i].vertexK , FixedVertex3d):
                pass
            else:
                raise Exception('Unknown type for vertexK: %s' % type(self.cells[i].vertexK))

            if isinstance(self.cells[i].vertexL , NodeVertex3d):
                if not self.cells[i].vertexL.key() in nodal_loads:
                    nodal_loads[self.cells[i].vertexL.key()] = 0.0    
            elif isinstance(self.cells[i].vertexL , ElementVertex3d):
                if not self.cells[i].vertexL.key() in element_loads:
                    element_loads[self.cells[i].vertexL.key()] = 0.0
            elif isinstance(self.cells[i].vertexL , FixedVertex3d):
                pass
            else:
                raise Exception('Unknown type for vertexL: %s' % type(self.cells[i].vertexL))

        self.empty_nodal_loads = nodal_loads;
        self.empty_element_loads = element_loads;
        self.committed_nodal_loads = nodal_loads.copy();
        self.committed_element_loads = element_loads.copy();
        
    def compute_current_load_vector(self,zw):
        nodal_loads = self.empty_nodal_loads.copy()
        element_loads = self.empty_element_loads.copy()
        
        for i in self.cells:  
            f = self.cells[i].get_load_vector(zw)
            
            if isinstance(self.cells[i].vertexI , NodeVertex3d):
                nodal_loads[self.cells[i].vertexI.key()] += f.item(0)
            elif isinstance(self.cells[i].vertexI , ElementVertex3d):
                element_loads[self.cells[i].vertexI.key()] += f.item(0)
            elif isinstance(self.cells[i].vertexI , FixedVertex3d):
                pass
            else:
                raise Exception('Unknown type for vertexI: %s' % type(self.cells[i].vertexI))
            
            if isinstance(self.cells[i].vertexJ , NodeVertex3d):
                nodal_loads[self.cells[i].vertexJ.key()] += f.item(1)
            elif isinstance(self.cells[i].vertexJ , ElementVertex3d):
                element_loads[self.cells[i].vertexJ.key()] += f.item(1)
            elif isinstance(self.cells[i].vertexJ , FixedVertex3d):
                pass
            else:
                raise Exception('Unknown type for vertexJ: %s' % type(self.cells[i].vertexJ))
            
            if isinstance(self.cells[i].vertexK , NodeVertex3d):
                nodal_loads[self.cells[i].vertexK.key()] += f.item(2)
            elif isinstance(self.cells[i].vertexK , ElementVertex3d):
                element_loads[self.cells[i].vertexK.key()] += f.item(2)
            elif isinstance(self.cells[i].vertexK , FixedVertex3d):
                pass
            else:
                raise Exception('Unknown type for vertexK: %s' % type(self.cells[i].vertexK))
            
            if isinstance(self.cells[i].vertexL , NodeVertex3d):
                nodal_loads[self.cells[i].vertexL.key()] += f.item(3)
            elif isinstance(self.cells[i].vertexL , ElementVertex3d):
                element_loads[self.cells[i].vertexL.key()] += f.item(3)
            elif isinstance(self.cells[i].vertexL , FixedVertex3d):
                pass
            else:
                raise Exception('Unknown type for vertexL: %s' % type(self.cells[i].vertexL))
                
        self.current_nodal_loads = nodal_loads
        self.current_element_loads = element_loads
        
    def commit_current_load_vector(self):
        self.committed_nodal_loads = self.current_nodal_loads.copy();
        self.committed_element_loads = self.current_element_loads.copy();

    def apply_load_increment(self):
        for i in self.committed_nodal_loads:        
            fy = self.current_nodal_loads[i] - self.committed_nodal_loads[i]
            ops.load(i, 0.0, 0.0, fy, 0.0, 0.0, 0.0)
        for i in self.committed_element_loads:
            # @todo - this code assumed local y axis is upward, update with axis call
            fy = self.current_element_loads[i] - self.committed_element_loads[i]
            iNodes = ops.eleNodes(i[0])
            xI,yI,zI = ops.nodeCoord(iNodes[0])
            xJ,yJ,zJ = ops.nodeCoord(iNodes[1])
            Lz = zJ-zI
            Lxy = sqrt((xJ-xI)**2 + (yJ-yI)**2)
            ele_angle = atan2(Lz, Lxy)
            ops.eleLoad('-ele', i[0], '-type', '-beamPoint', fy*cos(ele_angle), 0.0, i[1], fy*sin(ele_angle))
    