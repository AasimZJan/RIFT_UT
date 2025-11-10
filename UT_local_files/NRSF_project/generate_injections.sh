#!/bin/bash
# write mdc.xml.gz
util_WriteInjectionFile.py --parameter m1 --parameter-value 90.9090909090909 --parameter m2 --parameter-value 9.090909090909092 --parameter s1x --parameter-value 0.0 --parameter s1y --parameter-value 0.0 --parameter s1z --parameter-value 0.0 --parameter s2x --parameter-value 0.0 --parameter s2y --parameter-value 0.0 --parameter s2z --parameter-value 0.0 --parameter eccentricity --parameter-value 0 --approx SEOBNRv4 --parameter dist --parameter-value 530.6990386393425 --parameter fmin --parameter-value 15 --parameter incl --parameter-value 0.7 --parameter tref --parameter-value 1000000000 --parameter phiref --parameter-value 0.0 --parameter theta --parameter-value 0.1 --parameter phi --parameter-value 0.57 --parameter psi --parameter-value 0.0

# generate injections, will create plots with verbose flag
util_ROMWriteFrame.py --inj mdc.xml.gz --event 0 --instrument H1 --start 999999850 --stop 1000000150 --lmax 4 --verbose --srate 8192 --seglen 4 --verbose --group 'surrogate_downloads/' --param 'NRHybSur2dq15.h5'
util_ROMWriteFrame.py --inj mdc.xml.gz --event 0 --instrument L1 --start 999999850 --stop 1000000150 --lmax 4 --verbose --srate 8192 --seglen 4 --verbose --group 'surrogate_downloads/' --param 'NRHybSur2dq15.h5' 
util_ROMWriteFrame.py --inj mdc.xml.gz --event 0 --instrument V1 --start 999999850 --stop 1000000150 --lmax 4 --verbose --srate 8192 --seglen 4 --verbose --group 'surrogate_downloads/' --param 'NRHybSur2dq15.h5'

# calculate SNR, expects the psds to be outside one directory
ls *.gwf |lal_path2cache >local.cache; util_FrameZeroNoiseSNR.py --cache local.cache  --psd-file H1=../H1-psd.xml.gz --fmin-snr 20 --psd-file L1=../L1-psd.xml.gz --psd-file V1=../V1-psd.xml.gz --fmax-snr 4095 

# generate coinc.xml
util_SimInspiralToCoinc.py --sim-xml mdc.xml.gz --event 0





