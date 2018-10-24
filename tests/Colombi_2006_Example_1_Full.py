from PyPonding import FE

# Colombi, P. (2006). “The Ponding Problem on Flat Steel Roof Grids.” 
# Journal of Constructional Steel Research, 62(7), 647–655.

## Example 1 - Full ponding ##
E  = 2.1e8    # modulus of elasticity, kN/m^2
A  = 1000     # cross-sectional area, m^2
I  = 8.356e-5 # moment of inertia, m^2 (IPE 300)
L1 = 10       # member length, m
L2 = 5        # tributary width, m
h0 = 0.4      # water level, m
gamma = 9.8   # rainwater specific weight, kN/m^2

nele  = 40
nnode = nele+1

mdl = FE.Model('Example 1')

for i in range(nnode):
    mdl.AddNode('n%i' % (i+1),(L1*i/nele,0.0),('UX','UY','RZ'))

mdl.Nodes['n1'].dofs['UX'].constrained = True
mdl.Nodes['n1'].dofs['UY'].constrained = True
mdl.Nodes['n%i' % nnode].dofs['UY'].constrained = True

for i in range(nele):
    mdl.AddElement('e%i'%(i+1),'ElasticBeam2d',('n%i'%(i+1),'n%i'%(i+2)),E,I,A)
    mdl.AddPondingLoadCell('p%i'%(i+1),'2d',('n%i'%(i+1),'n%i'%(i+2)),gamma,L2)


res = FE.PondingAnalysis(mdl,'Constant_Level')
res.run({},h0)
    
res0 = FE.PondingAnalysis(mdl,'No_Ponding_Effect')
res0.run({},h0)

midnode = 'n%i' % (nele/2+1)

f  = mdl.Nodes[midnode].dofs['UY'].disp(res)
f0 = mdl.Nodes[midnode].dofs['UY'].disp(res0)

print('\nExample 1 (Full Ponding) from Colombi 2006')
print('   f = %.3f' % f)
print('  f0 = %.3f' % f0)
print('f/f0 = %.3f' % (f/f0))
print('f/f0 = 1.41 (analytical from paper)')
print('f/f0 = 1.40 (numerical from paper)\n')

