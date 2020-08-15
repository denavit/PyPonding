from PyPonding import PondingLoadCell3d
  
# Define units
inch = 1.0
kip = 1.0
lb = kip/1000.0
ft = 12.0*inch
pcf = lb/ft**3

# Define Load Cell
cell = PondingLoadCell3d()

# Define vertices in counterclockwise direction
cell.xI =  0.0*inch
cell.yI =  0.0*inch
cell.zI =  2.0*inch
cell.xJ = 10.0*inch
cell.yJ =  0.0*inch
cell.zJ =  0.0*inch
cell.xK = 10.0*inch
cell.yK = 20.0*inch
cell.zK =  0.0*inch
cell.xL =  0.0*inch
cell.yL = 20.0*inch
cell.zL =  2.0*inch

cell.dzI = 0.0*inch
cell.dzJ = 0.0*inch
cell.dzK = 0.0*inch
cell.dzL = 0.0*inch

cell.gamma = 62.4*pcf

cell.na = 4
cell.nb = 4

print(cell.get_load_vector(1.0))

print(cell.get_volume(1.0))