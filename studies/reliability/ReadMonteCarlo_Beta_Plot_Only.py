import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
#from matplotlib import rc

import sys
sys.path.append('/home/mhscott/PyPonding/')
sys.path.append('/home/mhscott/PyPonding/PyPonding/')
sys.path.append('/home/mhscott/PyPonding/PyPonding/structures/')
from PyPonding.structures import wide_flange
from wide_flange import wf,wf_shapes

import openseespy.opensees as ops

from math import pi,log

y_value_for_infty = 4.4
ylim_max = y_value_for_infty + 0.2

# Define units
inch = 1.0
kip = 1.0
minute = 1.0

lb = kip/1000.0
ft = 12.0*inch
hr = 60*minute
sec = minute/60.0

g = 386.4*inch/sec**2

ksi = kip/inch**2
psf = lb/ft**2
pcf = psf/ft

# Parametric Set
list_spans_and_sections = [
    (20*ft, 'W8X10'),(20*ft,'W10X12'),(20*ft,'W12X16'),(20*ft,'W14X22'),
    (30*ft,'W12X16'),(30*ft,'W14X22'),(30*ft,'W16X26'),(30*ft,'W18X35'),
    (40*ft,'W16X26'),(40*ft,'W18X35'),(40*ft,'W21X44'),(40*ft,'W24X55'),
    (50*ft,'W21X44'),(50*ft,'W24X55'),(50*ft,'W27X84'),(50*ft,'W30X90')]
list_slope = [0.0*inch/ft,0.25*inch/ft,0.5*inch/ft,1.0*inch/ft]
list_qD = [10.0*psf,20.0*psf,30.0*psf]

methods = ['AISC Appendix 2','DAMP','Proposed for ASCE 7','Neglect Ponding']
Nmethods = len(methods)

colorMethod = {}
colorMethod['AISC Appendix 2'] = 'r1'
colorMethod['DAMP'] = 'b2'
colorMethod['Proposed for ASCE 7'] = 'g_'
colorMethod['Neglect Ponding'] = 'k|'

locations = ['Denver','New York','New Orleans']
Ncities = len(locations)

#
# A fix to correct incorrect input to OpenSees Type I largest value PDF without
# having to re-run all 192 Monte Carlo simulations
#
rate = {}
perc95 = {}
rate['Denver'] = 1.26*inch/(0.25*hr)
perc95['Denver'] = 1.72*inch/(0.25*hr)
rate['New York'] = 1.67*inch/(0.25*hr)
perc95['New York'] = 2.33*inch/(0.25*hr)
rate['New Orleans'] = 2.33*inch/(0.25*hr)
perc95['New Orleans'] = 3.11*inch/(0.25*hr)

scale = {}
loc = {}
icity = 0; Ncities = len(locations)
for city in locations:
        icity += 1
        eulergamma = 0.57721566490153286061
        scale[city] = (rate[city]-perc95[city])/(log(-log(0.95))+eulergamma)
        loc[city] = rate[city] - scale[city]*eulergamma
        # Old random variable
        ops.randomVariable(icity,'type1LargestValue','-mean',rate[city],'-stdv',(scale[city]*6**0.5)/pi)
        # New random variable
        ops.randomVariable(icity+Ncities,'type1LargestValue','-parameters',loc[city],1/scale[city])

        
betaMethod = {}
pfMethod = {}
for city in locations:
        for method in methods:
                betaMethod[method,city] = []
                pfMethod[method,city] = []
spanToDepth = []

icase = 0
# Run Analysis
for span_and_section in list_spans_and_sections:
        L = span_and_section[0]
        shape_name = span_and_section[1]
        wfs = wf_shapes[shape_name]
        d = wfs['d']*inch

        for slope in list_slope:
                for qD in list_qD:

                        spanToDepth = np.append(spanToDepth,L/d)

                        icase += 1
                        print(icase)
       
                        df = np.loadtxt(f'../../tests/trials{icase}.csv',skiprows=2,delimiter=',')
                        [Ntrials,c] = df.shape

                        #
                        # Modify the data as a workaround for the OpenSees random variable input issue
                        icity = 0
                        for city in locations:
                                for j in range(Ntrials):
                                        # Read old dh, then calculate intensity
                                        dhold = df[j,Nmethods+Ncities+icity]
                                        q = (0.6*12*inch*(2*g)**0.5*dhold**1.5)/1.5
                                        As = (40*ft)*L
                                        intensity = q/As

                                        # Convert to new intensity based on CDFs
                                        p = ops.getCDF(icity+1, intensity)
                                        intensity = ops.getInverseCDF(icity+1+Ncities, p)

                                        # Calculate new flow and dh
                                        q = intensity*As
                                        dhnew = (1.5*q/(0.6*12*inch*(2*g)**0.5))**(2.0/3)
                                        df[j,Nmethods+Ncities+icity] = dhnew
                                icity += 1
                        # End modification
                        #
                        
                        icity = 0
                        for city in locations:
                                print(city)
                                imethod = 0
                                for method in methods:
                                        if method == 'AISC Appendix 2' and slope != 0.0:
                                                pfMethod[method,city] = np.append(pfMethod[method,city],1.0)
                                                betaMethod[method,city] = np.append(betaMethod[method,city],-5.0)
                                                imethod += 1
                                                continue
                                        ds = df[:,imethod] - df[:,Nmethods+icity]
                                        x = ds + df[:,Nmethods+Ncities+icity] > df[:,Nmethods+Ncities+Ncities]

                                        Nfail = np.count_nonzero(x == True)
                                        pf = Nfail/Ntrials
                                        pfMethod[method,city] = np.append(pfMethod[method,city],pf)
                                        if pf > 0.0:
                                                beta = -norm.ppf(pf)
                                        else:
                                                beta = y_value_for_infty
                                        betaMethod[method,city] = np.append(betaMethod[method,city],beta)
                                        print(f'  Method {imethod+1}, pf={pf}, beta={beta}')
                                        imethod += 1
                                icity += 1

# Plot Results
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=8)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)

fig = plt.figure(figsize=(7.00,5.00))


method = methods[0]
ax = fig.add_axes([0.075,0.58,0.40,0.37])
for city in locations:
        plt.plot(spanToDepth,betaMethod[method,city],'k2',label=method)
plt.plot([0,32],[2.6,2.6],'k--')
plt.xlim([16,32])
plt.ylim([0,ylim_max])
ax.set_yticks([0,1,2,3,3.7,y_value_for_infty])
ax.set_yticklabels(['0','1','2','3','3.7','$>$3.7'])
ax.set_xlabel('Span-to-Depth Ratio ($L/d$)')
ax.set_ylabel('Reliability Index ($\\beta$)')
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.set_title('(a) AISC Appendix 2')

method = methods[1]
ax = fig.add_axes([0.575,0.58,0.40,0.37])
for city in locations:
        plt.plot(spanToDepth,betaMethod[method,city],'k2',label=method)
plt.plot([0,32],[2.6,2.6],'k--')
plt.xlim([16,32])
plt.ylim([0,ylim_max])
ax.set_yticks([0,1,2,3,3.7,y_value_for_infty])
ax.set_yticklabels(['0','1','2','3','3.7','$>$3.7'])
ax.set_xlabel('Span-to-Depth Ratio ($L/d$)')
ax.set_ylabel('Reliability Index ($\\beta$)')
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.set_title('(b) DAMP')

method = methods[2]
ax = fig.add_axes([0.075,0.08,0.40,0.37])
for city in locations:
        plt.plot(spanToDepth,betaMethod[method,city],'k2',label=method)
plt.plot([0,32],[2.6,2.6],'k--')
plt.xlim([16,32])
plt.ylim([0,ylim_max])
ax.set_yticks([0,1,2,3,3.7,y_value_for_infty])
ax.set_yticklabels(['0','1','2','3','3.7','$>$3.7'])
ax.set_xlabel('Span-to-Depth Ratio ($L/d$)')
ax.set_ylabel('Reliability Index ($\\beta$)')
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.set_title('(c) Modified Rain Load')

method = methods[3]
ax = fig.add_axes([0.575,0.08,0.40,0.37])
for city in locations:
        plt.plot(spanToDepth,betaMethod[method,city],'k2',label=method)
plt.plot([0,32],[2.6,2.6],'k--')
plt.xlim([16,32])
plt.ylim([0,ylim_max])
ax.set_yticks([0,1,2,3,3.7,y_value_for_infty])
ax.set_yticklabels(['0','1','2','3','3.7','$>$3.7'])
ax.set_xlabel('Span-to-Depth Ratio ($L/d$)')
ax.set_ylabel('Reliability Index ($\\beta$)')
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.set_title('(d) Neglect Ponding')

plt.savefig(f'betaSpanToDepth.png',dpi=300)
plt.savefig(f'betaSpanToDepth.pdf')

plt.show()
