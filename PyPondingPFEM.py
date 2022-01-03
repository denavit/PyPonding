import opensees as ops
from math import ceil
import os
import matplotlib.pyplot as plt

# Define units
m = 1.0
kg = 1.0
sec = 1.0
Newton = kg*m/sec**2
Pa = Newton / m**2
inch = 0.0254 * m
ft = 12.0*inch
in_per_ft = inch/ft
ksi = 6.895e+6 * Pa
psf = ksi / 1000.0
pcf = psf/ft

# structure parameters
Fy = 50.0*ksi
E = 29000.0*ksi
Hk = 1.0e-4*E
TW = 10*ft
shape_name = 'W14X22'
L = 40*ft
H = 20*ft
slope = 0.0*in_per_ft
qD = 20.0*psf
gamma = 62.4*pcf
wD = qD*TW
shape = {'W': 22, 'A': 6.49, 'd': 13.7, 'ddet': 13.75, 'bf': 5, 'bfdet': 5, 'tw': 0.23, 'twdet': 0.25, 'twdet/2': 0.125, 'tf': 0.335, 'tfdet': 0.3125, 'kdes': 0.735, 'kdet': 1.0625, 'k1': 0.75, 'bf/2tf': 7.46, 'h/tw': 53.3, 'Ix': 199, 'Zx': 33.2, 'Sx': 29,
         'rx': 5.54, 'Iy': 7, 'Zy': 4.39, 'Sy': 2.8, 'ry': 1.04, 'J': 0.208, 'Cw': 314, 'Wno': 16.7, 'Sw1': 7, 'Qf': 5.34, 'Qw': 16.1, 'rts': 1.27, 'ho': 13.4, 'PA': 41.3, 'PB': 46.3, 'PC': 32.4, 'PD': 37.4, 'T': 11.625, 'WGi': 2.75, 'metric_label': 'W360X32.9'}

num_elements = 20
num_fiber = 20
eleLen = L / num_elements

# A = shape["A"] * inch**2
d = shape["d"] * inch
tw = shape["tw"] * inch
bf = shape["bf"] * inch
tf = shape["tf"] * inch

# fluid parameters
numx = 3.0
numy = 3.0

rho = 1000.0 * kg/m**3
mu = 0.0001 * Newton * sec/m**2
b1 = 0.0
b2 = -9.81 * m / sec**2
thk = TW
kappa = -1.0

bmass = wD * eleLen / abs(b2)

# analysis parameters
dtmax = 1e-3
dtmin = 1e-3
totaltime = 1.0
filename = "ponding"

# model
ops.wipe()

ops.model('basic', '-ndm', 2, '-ndf', 3)

# recorder
ops.recorder('PVD', filename, 'disp', 'vel', 'pressure', '-dT', 1e-3)
if not os.path.exists(filename):
  os.makedirs(filename)

# fluid mesh
ndf = 3

# nx = round(L / eleLen * numx)
# ny = round(H / h * numy)

# wall mesh
ops.node(1, 0.0, H)
ops.node(2, 0.0, 0.0)
ops.node(3, L, 0.0)
ops.node(4, L, H)

sid = 1
walltag1 = 1
ops.mesh('line', walltag1, 2, 1, 2, sid, ndf, eleLen)

walltag2 = 2
ops.mesh('line', walltag2, 2, 3, 4, sid, ndf, eleLen)

wallNodes1 = ops.getNodeTags('-mesh', walltag1)
wallNodes2 = ops.getNodeTags('-mesh', walltag2)

for nd in wallNodes1:
  ops.fix(nd, 1, 1, 1)
for nd in wallNodes2:
  ops.fix(nd, 1, 1, 1)

# geom transf
transfTag = 1
ops.geomTransf('Corotational', transfTag)

# material
matTag = 1
ops.uniaxialMaterial('Elastic', matTag, E)

# section
dw = d - 2*tf
Nfw = ceil(dw*(num_fiber/d))
Nff = ceil(tf*(num_fiber/d))

secTag = 1
ops.section('WFSection2d', secTag, matTag,
            d, tw, bf, tf, Nfw, Nff)

# beam integration
inteTag = 1
Npts = 3
ops.beamIntegration('Lobatto', inteTag, secTag, Npts)

# beam mesh
beamtag = 3
eleArgs = ["dispBeamColumn", transfTag, inteTag]
ops.mesh('line', beamtag, 2, 2, 3, sid, ndf, eleLen, *eleArgs)

sNodes = ops.getNodeTags('-mesh', beamtag)

sEles = ops.getEleTags('-mesh', beamtag)

tsTag = 1
ops.timeSeries('Constant', tsTag)

patternTag = 1
ops.pattern('Plain', patternTag, tsTag)
for ele in sEles:
  ops.eleLoad('-ele', ele, '-type', 'beamUniform', -wD)

for nd in sNodes:
  ops.mass(nd, bmass, bmass, 0.0)

# create constraint object
ops.constraints('Plain')

# create numberer object
ops.numberer('Plain')

# create algorithm object
ops.algorithm('Newton')

# create integrator object
ops.integrator("LoadControl", 0.0)

# create SOE object
ops.system("Mumps")
# system('PFEM')
# ops.system('PFEM', '-mumps')

# create analysis object
ops.analysis('Static')
# ops.analysis('PFEM', dtmax, dtmin, b2)

# run dead load analysis
ops.analyze(1)

for nd in sNodes:
  print(f'node {nd} disp = {ops.nodeDisp(nd)}')
