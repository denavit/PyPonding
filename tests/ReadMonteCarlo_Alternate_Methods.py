import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
#from matplotlib import rc
from wide_flange import wf,wf_shapes

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

ksi = kip/inch**2
psf = lb/ft**2
pcf = psf/ft

# Constants
Fy = 50.0*ksi # yield stress
E  = 29000.0*ksi # elastic modulus
Hk = 29.0*ksi # kinematic hardening modulus
tw = 10.0*ft # Tributary width
gamma = 62.4*pcf

# Parametric Set
list_spans_and_sections = [
    (20*ft, 'W8X10'),(20*ft,'W10X12'),(20*ft,'W12X16'),(20*ft,'W14X22'),
    (30*ft,'W12X16'),(30*ft,'W14X22'),(30*ft,'W16X26'),(30*ft,'W18X35'),
    (40*ft,'W16X26'),(40*ft,'W18X35'),(40*ft,'W21X44'),(40*ft,'W24X55'),
    (50*ft,'W21X44'),(50*ft,'W24X55'),(50*ft,'W27X84'),(50*ft,'W30X90')]
list_slope = [0.0*inch/ft,0.25*inch/ft,0.5*inch/ft,1.0*inch/ft]
list_qD = [10.0*psf,20.0*psf,30.0*psf]
#methods = ['AISC Appendix 2','DAMP','Modified Rain Load','Neglect Ponding']
#methods = ['DAMP without 0.8','Modified Rain Load with 0.8','DAMP 1.4','DAMP 1.6']
methods = ['DAMP without 0.8','Modified Rain Load with 0.8']

locations = ['Denver','New York','New Orleans']


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
        
                        df = np.loadtxt(f'trials{icase}.csv',skiprows=2,delimiter=',')
                        [Ntrials,c] = df.shape

                        icity = 0
                        Ncities = len(locations)
                        for city in locations:
                                print(city)
                                imethod = 0
                                Nmethods = len(methods)
                                for method in methods:
                                        if method == 'AISC Appendix 2' and slope != 0.0:
                                                pfMethod[method,city] = np.append(pfMethod[method,city],1.0)
                                                betaMethod[method,city] = np.append(betaMethod[method,city],-5.0)
                                                imethod += 1
                                                continue
                                        
                                        # Read from file
                                        # zwlim = df[:,imethod];
                                        
                                        # Compute new 
                                        wf_section = wf(wfs['d'],wfs['tw'],wfs['bf'],wfs['tf'],Fy,E,Hk)
                                        wf_section.gamma = gamma
                                        wf_section.TW    = tw
                                        wf_section.wD    = qD*tw
                                        wf_section.L     = L
                                        wf_section.zi    = 0.0
                                        wf_section.zj    = slope*L                                   
                                        zwlim = wf_section.maximum_permitted_zw(method)
                                        # print(zwlim)                                        
                                        
                                        ds = zwlim - df[:,4+icity]
                                        x = ds + df[:,4+Ncities+icity] > df[:,4+Ncities+Ncities]

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


fig = plt.figure(figsize=(7.00,2.50))

method = methods[0]
ax = fig.add_axes([0.075,0.18,0.40,0.73])
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
ax.set_title('(a) DAMP (without stiffness reduction)')

method = methods[1]
ax = fig.add_axes([0.575,0.18,0.40,0.73])
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
ax.set_title('(b) Modified Rain Load (with stiffness reduction)')

plt.savefig('betaSpanToDepth_alt_stiff.png',dpi=300)
plt.savefig('betaSpanToDepth_alt_stiff.pdf')

plt.show()
