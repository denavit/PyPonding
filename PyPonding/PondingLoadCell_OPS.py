from math import sin,cos,atan2
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
    def __init__(self,node_id):
        self.node_id = node_id
       
    def key(self):
        return self.node_id
        
    def coord(self):
        x = ops.nodeCoord(self.node_id,1)
        y = ops.nodeCoord(self.node_id,2)
        return (x,y)
    
    def disp(self):
        dx = ops.nodeDisp(self.node_id,1)
        dy = ops.nodeDisp(self.node_id,2)
        return (dx,dy)
    
class ElementEnd2d:
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
        y = yI + self.x*(yJ-yI)
        return (x,y)
    
    def disp(self):
        import warnings
        warnings.warn('PondingLoadCell not fully implemented for ends within and element. Displacement will be returned as zero.')
        dx = 0.
        dy = 0.
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
        x,y = self.endI.coord()
        self.xI = x
        self.yI = y
        x,y = self.endJ.coord()
        self.xJ = x
        self.yJ = y
        
        # Store node ids (if attached to nodes) for backwards compatibility 
        if isinstance(self.endI, NodeEnd2d):
            self.nodeI = self.endI.node_id
        else:
            self.nodeI = None
        if isinstance(self.endJ, NodeEnd2d):
            self.nodeJ = self.endJ.node_id
        else:
            self.nodeJ = None        
        
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
            iNodes = ops.eleNodes(0)
            xI,yI = ops.nodeCoord(iNodes[0])
            xJ,yJ = ops.nodeCoord(iNodes[1])
            ele_angle = atan2(yJ-yI, xJ-xI)
            ops.eleLoad('-ele', i[0], '-type', '-beamPoint', fy*cos(ele_angle), i[1], fy*sin(ele_angle))
            
class PondingLoadCell3d_OPS(PondingLoadCell3d):
    def __init__(self, id, vertexI, nodeJ, nodeK, nodeL, gamma, na=1, nb=1):
        self.id = id
        self.nodeI = nodeI
        self.nodeJ = nodeJ
        self.nodeK = nodeK
        self.nodeL = nodeL
        self.gamma = gamma
        self.na = na
        self.nb = nb

        self.xI, self.yI, self.zI = ops.nodeCoord(self.nodeI)
        self.xJ, self.yJ, self.zJ = ops.nodeCoord(self.nodeJ)
        self.xK, self.yK, self.zK = ops.nodeCoord(self.nodeK)
        self.xL, self.yL, self.zL = ops.nodeCoord(self.nodeL)

    def update(self):
        self.dzI = ops.nodeDisp(self.nodeI, 3)
        self.dzJ = ops.nodeDisp(self.nodeJ, 3)
        self.dzK = ops.nodeDisp(self.nodeK, 3)
        self.dzL = ops.nodeDisp(self.nodeL, 3)
