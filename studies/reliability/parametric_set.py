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


# Parametric Set
list_spans_and_sections = [
    (20*ft, 'W8X10'),(20*ft,'W10X12'),(20*ft,'W12X16'),(20*ft,'W14X22'),
    (30*ft,'W12X16'),(30*ft,'W14X22'),(30*ft,'W16X26'),(30*ft,'W18X35'),
    (40*ft,'W16X26'),(40*ft,'W18X35'),(40*ft,'W21X44'),(40*ft,'W24X55'),
    (50*ft,'W21X44'),(50*ft,'W24X55'),(50*ft,'W27X84'),(50*ft,'W30X90')]
list_slope = [0.0*in_per_ft,0.25*in_per_ft,0.5*in_per_ft,1.0*in_per_ft]
list_qD = [10.0*psf,20.0*psf,30.0*psf]

parametric_set = []

# Run Analysis
for span_and_section in list_spans_and_sections:
    L = span_and_section[0]
    shape_name = span_and_section[1]
    for slope in list_slope:
        for qD in list_qD:
            
            # Add case to list
            parametric_set.append({
                "shape": shape_name,
                "L": L,
                "slope": slope,
                "qD": qD
            })


# Print
f = open('Parametric_Set.csv','w')
f.write('shape_name,L (ft),slope (in/ft),dead load (psf)\n')
for i in parametric_set: 
    f.write('%s,%i,%.4f,%.4f,%.4f,%.4f\n' % (i['shape'],i['L']/ft,i['slope']/in_per_ft,i['qD']/psf))
f.close
    