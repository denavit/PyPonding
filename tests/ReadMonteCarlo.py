import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text',usetex=True)
rc('font',family='serif')

from wide_flange import wf,wf_shapes

from pylatex import *

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

# Parametric Set
list_spans_and_sections = [
    (20*ft, 'W8X10'),(20*ft,'W10X12'),(20*ft,'W12X16'),(20*ft,'W14X22'),
    (30*ft,'W12X16'),(30*ft,'W14X22'),(30*ft,'W16X26'),(30*ft,'W18X35'),
    (40*ft,'W16X26'),(40*ft,'W18X35'),(40*ft,'W21X44'),(40*ft,'W24X55'),
    (50*ft,'W21X44'),(50*ft,'W24X55'),(50*ft,'W27X84'),(50*ft,'W30X90')]
list_slope = [0.0*inch/ft,0.25*inch/ft,0.5*inch/ft,1.0*inch/ft]
list_qD = [10.0*psf,20.0*psf,30.0*psf]
methods = ['AISC Appendix 2','DAMP','Proposed for ASCE 7','Neglect Ponding']
colorMethod = {}
colorMethod['AISC Appendix 2'] = 'r1'
colorMethod['DAMP'] = 'b2'
colorMethod['Proposed for ASCE 7'] = 'g_'
colorMethod['Neglect Ponding'] = 'k|'

doc = Document('MonteCarloRoofPonding',documentclass='report')

betaMethod = {}
pfMethod = {}
for method in methods:
        betaMethod[method] = []
        pfMethod[method] = []
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

                        doc.append(NoEscape(r'\clearpage'))
        
                        plt.figure()
                        xy = np.loadtxt(f'meanHV{icase}.csv',delimiter=',')
                        plt.plot(xy[:,0],xy[:,1],'-k')
                        plt.xlabel('Volume (in$^3$)')
                        plt.ylabel('Height (in)')
                        plt.title(f'{shape_name}, {L/ft} ft, {slope/(inch/ft)} in/ft, {qD/psf} psf')
                        plt.grid()
                        plt.ylim(bottom=-5)
                        plt.savefig(f'meanHV{icase}.pdf')
                        with doc.create(Figure(position="htpb")) as plot:
                                plot.add_plot(width=NoEscape(r'0.8\textwidth'))
                        plt.close()
        
                        df = np.loadtxt(f'trials{icase}.csv',skiprows=2,delimiter=',')
                        [Ntrials,c] = df.shape

                        imethod = 0
                        for method in methods:
                                x = df[:,imethod] + df[:,4] > df[:,5]

                                Nfail = np.count_nonzero(x == True)
                                pf = Nfail/Ntrials
                                pfMethod[method] = np.append(pfMethod[method],pf)
                                beta = -norm.ppf(pf)
                                betaMethod[method] = np.append(betaMethod[method],beta)
                                print(f'Method {imethod+1}, pf={pf}, beta={beta}')
                                doc.append(NoEscape(f'{method}, pf={pf}, beta={beta}\n\n'))
                                imethod += 1


plt.figure()
plt.subplot(2,1,1)
for method in methods:
        plt.plot(spanToDepth,betaMethod[method],f'{colorMethod[method]}',label=method)
#plt.xlabel('Span to Depth (L/d)')
plt.ylabel('Reliability Index')
plt.xlim(left=16)
plt.ylim([0,3.5])
plt.legend()
#plt.grid()
plt.subplot(2,1,2)
for method in methods:
        plt.plot(spanToDepth,pfMethod[method],f'{colorMethod[method]}',label=method)
plt.xlabel('Span to Depth ($L/d$)')
plt.ylabel('Probability of Failure')
plt.xlim(left=16)
plt.ylim(top=0.4)
plt.legend()
#plt.grid()

plt.savefig(f'allSpanToDepth.pdf')
with doc.create(Figure(position="htpb")) as plot:
        plot.add_plot(width=NoEscape(r'\textwidth'))
plt.close()



doc.append(NoEscape(r'\clearpage'))



plt.figure()
im = 1
for method in methods:
        plt.subplot(2,2,im)
        plt.plot(spanToDepth,betaMethod[method],f'{colorMethod[method]}',label=method)
        plt.xlim(left=16)
        plt.ylim([0,3.5])
        if im == 3 or im == 4:
                plt.xlabel('Span to Depth ($L/d$)')
        if im == 1 or im == 3:
                plt.ylabel('Reliability Index')
        plt.title(method)
        im += 1

plt.savefig(f'betaSpanToDepth.pdf')
with doc.create(Figure(position="htpb")) as plot:
        plot.add_plot(width=NoEscape(r'\textwidth'))
plt.close()


plt.figure()
im = 1
for method in methods:
        plt.subplot(2,2,im)
        plt.plot(spanToDepth,pfMethod[method],f'{colorMethod[method]}',label=method)
        plt.xlim(left=16)
        plt.ylim(top=0.01)
        if im == 4:
                plt.ylim(top=0.4)
        if im == 3 or im == 4:
                plt.xlabel('Span to Depth ($L/d$)')
        if im == 1 or im == 3:
                plt.ylabel('Probability of Failure')
        plt.title(method)
        im += 1

plt.savefig(f'pfSpanToDepth.pdf')
with doc.create(Figure(position="htpb")) as plot:
        plot.add_plot(width=NoEscape(r'\textwidth'))
plt.close()

doc.generate_pdf(clean_tex=False)
