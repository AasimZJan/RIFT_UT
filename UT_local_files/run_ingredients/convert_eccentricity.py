import sys
import RIFT.lalsimutils as lsu
import numpy as np

grid = sys.argv[1]
eccentricity = float(sys.argv[2])
Plist = lsu.xml_to_ChooseWaveformParams_array(grid)

eccentricity_previous = []
for i in np.arange(len(Plist)):
    P = Plist[i]
    eccentricity_previous.append(P.eccentricity)
eccentricity_previous=np.array(eccentricity_previous)
med_ecc, std = np.median(eccentricity_previous), np.std(eccentricity_previous)
print('Median eccentricity = {med_ecc}, std = {std}')
n = 0
for i in np.arange(len(Plist)):
    P = Plist[i]
    ecc_new = P.eccentricity - med_ecc + eccentricity
    if ecc_new < 0.0:
        scaling_factor = 6 # just to be safe
        ecc_new = np.abs(np.random.normal(0, std*scaling_factor, 1))[0]
    print(f'Converting {P.eccentricity} to {ecc_new}')
    P.eccentricity = ecc_new
    n+=1

lsu.ChooseWaveformParams_array_to_xml(Plist,"overlap-grid.xml.gz")
