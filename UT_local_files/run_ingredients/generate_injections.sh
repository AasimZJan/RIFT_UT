#!/bin/bash
# write mdc.xml.gz
util_WriteInjectionFile.py  --parameter m1 --parameter-value 8.741082482004446 --parameter m2 --parameter-value 2.162761279562819  --parameter s1x --parameter-value 0.0 --parameter s1y --parameter-value 0.0 --parameter s1z --parameter-value -0.07443819135153053 --parameter s2x --parameter-value 0.0 --parameter s2y --parameter-value 0.0 --parameter s2z --parameter-value -0.03167082572099747  --parameter eccentricity --parameter-value 0.1 --approx SEOBNRv5EHM --parameter dist --parameter-value 502.4800559609168 --parameter fmin --parameter-value 19 --parameter incl --parameter-value 0.0  --parameter tref --parameter-value 1000000000.0 --parameter phiref --parameter-value 3.1771242031241766  --parameter theta --parameter-value -0.12906498250155196 --parameter phi --parameter-value 1.7021502910723427  --parameter psi --parameter-value 1.1864755949271302

# generate injections, will create plots with verbose flag
util_GWSignalWriteFrame.py --inj mdc.xml.gz --event 0 --instrument H1 --start 999999850 --stop 1000000150 --approx "SEOBNRv5EHM" --l-max 4 --verbose --srate 8192 --seglen 32 --fref 19 --verbose 
util_GWSignalWriteFrame.py --inj mdc.xml.gz --event 0 --instrument L1 --start 999999850 --stop 1000000150 --approx "SEOBNRv5EHM" --l-max 4 --verbose --srate 8192 --seglen 32 --fref 19 --verbose 
util_GWSignalWriteFrame.py --inj mdc.xml.gz --event 0 --instrument V1 --start 999999850 --stop 1000000150 --approx "SEOBNRv5EHM" --l-max 4 --verbose --srate 8192 --seglen 32 --fref 19 --verbose 

# calculate SNR, expects the psds to be outside one directory
ls *.gwf |lal_path2cache >local.cache; util_FrameZeroNoiseSNR.py --cache local.cache  --psd-file H1=../H1-psd.xml.gz --fmin-snr 20 --psd-file L1=../L1-psd.xml.gz --psd-file V1=../V1-psd.xml.gz --fmax-snr 4095 

# generate coinc.xml
util_SimInspiralToCoinc.py --sim-xml mdc.xml.gz --event 0





