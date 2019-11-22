import numpy as np
from scipy.stats import norm

df = np.loadtxt('trials.csv',skiprows=1,delimiter=',')

[Ntrials,c] = df.shape

for method in range(3):
    x = df[:,method] + df[:,3] > df[:,4]

    Nfail = np.count_nonzero(x == True)
    pf = Nfail/Ntrials
    beta = -norm.ppf(pf)
    
    print('Method %d, pf=%g, beta=%g' % (method,pf,beta))
