import openseespy.opensees as ops
from math import pi,cos,cosh,ceil
import numpy as np
import matplotlib.pyplot as plt
from wide_flange import wf,wf_shapes

  
# Define units
inch = 1.0
kip = 1.0

lb = kip/1000.0
ft = 12.0*inch

in_per_ft = inch/ft

ksi = kip/inch**2
psf = lb/ft**2
pcf = psf/ft


# Input parameters
Fy      = 50.0*ksi
E       = 29000.0*ksi
Hk      = E/10000.0

TW      = 10.0*ft

list_spans_and_sections = [
    (20*ft, 'W8X10'),(20*ft,'W10X12'),(20*ft,'W12X16'),(20*ft,'W12X22'),(20*ft,'W14X22'),
    (30*ft,'W12X16'),(30*ft,'W12X22'),(30*ft,'W14X22'),(30*ft,'W16X26'),(30*ft,'W18X35'),
    (40*ft,'W14X22'),(40*ft,'W16X26'),(40*ft,'W18X35'),(40*ft,'W21X44'),(40*ft,'W24X55'),
    (50*ft,'W18X35'),(50*ft,'W21X44'),(50*ft,'W24X55'),(50*ft,'W27X84'),(50*ft,'W30X90')]
list_slope = [0.0*in_per_ft,0.25*in_per_ft,0.5*in_per_ft,1.0*in_per_ft]
list_qD = [10.0*psf,20.0*psf,30.0*psf]


# Run Analysis
f = open('StrengthEvaluationOutput.csv','w')
f.write('shape_name,L (ft),slope (in/ft),dead load (psf),zmax_inelastic,zmax_elastic,tau,zmax_elastic_reduced\n')
for span_and_section in list_spans_and_sections:
    L = span_and_section[0]
    shape_name = span_and_section[1]
    for slope in list_slope:
        for qD in list_qD:
    
            # Lookup shape data
            shape_data = wf_shapes[shape_name]
            d  = shape_data['d']*inch
            bf = shape_data['bf']*inch
            tf = shape_data['tf']*inch
            tw = shape_data['tw']*inch

            # Create wide-flange object
            wf_section = wf(d,tw,bf,tf,Fy,E,Hk)

            # Set additional properties
            wf_section.L        = L 
            wf_section.gamma    = 62.4*pcf
            wf_section.TW       = TW
            wf_section.zi       = 0.0*inch
            wf_section.zj       = L*slope
            wf_section.wD       = qD*TW # Dead load per length

            # Set OpenSees analysis options
            wf_section.material_type = 'Hardening'
            wf_section.max_volume = (500*inch)*wf_section.L*wf_section.TW
            wf_section.num_steps = 10000
            wf_section.vol_tol = wf_section.max_volume/wf_section.num_steps/10000.

            # Run OpenSees analysis 
            (data_volume,data_height) = wf_section.perform_OpenSees_analysis();
            zmax_inelastic = np.max(data_height)


            # Elastic Analyses
            elastic_beam = wf_section.steel_beam_object()
            
            elastic_beam.alpha = 1.0
            elastic_beam.LF_D  = 1.0
            elastic_beam.LF_S2 = 0.0
            zmax_elastic = elastic_beam.Run_To_Strength_Limit(start_level = 1.0,max_level=500.0)
            tau = elastic_beam.determine_stiffness_reduction(zmax_inelastic)
                        
            elastic_beam.E = 0.8*elastic_beam.E
            zmax_elastic_reduced = elastic_beam.Run_To_Strength_Limit(start_level = 1.0,max_level=500.0)
            
            # Print Data
            f.write('%s,%i,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n' % (shape_name,L/ft,slope/in_per_ft,qD/psf,zmax_inelastic,zmax_elastic,tau,zmax_elastic_reduced))
f.close()
