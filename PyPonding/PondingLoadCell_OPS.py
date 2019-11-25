from PyPonding.PondingLoadCell import PondingLoadCell2d, PondingLoadCell3d
import openseespy.opensees as ops

class PondingLoadCell2d_OPS(PondingLoadCell2d):
    def __init__(self,id,nodeI,nodeJ,gamma,tw):
        self.id = id
        self.nodeI = nodeI
        self.nodeJ = nodeJ
        self.gamma = gamma
        self.tw = tw
        
        # Retreive Node Coordinates
        self.xI = ops.nodeCoord(self.nodeI,1)
        self.yI = ops.nodeCoord(self.nodeI,2)
        self.xJ = ops.nodeCoord(self.nodeJ,1)
        self.yJ = ops.nodeCoord(self.nodeJ,2)
        
    def update(self):
        # Code currently only updates y postion of nodes - @todo maybe update x position also
        # self.dxI = ops.nodeDisp(self.nodeI,1)
        self.dyI = ops.nodeDisp(self.nodeI,2)
        # self.dxJ = ops.nodeDisp(self.nodeJ,1)
        self.dyJ = ops.nodeDisp(self.nodeJ,2)


class PondingLoadCell3d_OPS(PondingLoadCell3d):
    def __init__(self, id, nodeI, nodeJ, nodeK, nodeL, gamma, na=1, nb=1):
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
