#! /usr/bin/env python
#  
#
# Tool: dumps a candidate waveform to a frame file.
# Default GPS time used is boring
# Tries to be as close to C as possible -- no interface via pylal/glue
#
# EXAMPLE
#    python util_LALWriteFrame.py; FrDump -i dummy.gwf
#
#  WARNING: My version does NOT interpolate the signal to a synchronized set of sample times.
#                       This may cause problems for some applications, particularly at low sample rates.


import argparse
import numpy as np
import RIFT.lalsimutils as lalsimutils
import lalsimulation as lalsim
import lalframe
import lal
import RIFT.physics.GWSignal as gwsignal
#import pyseobnr

parser = argparse.ArgumentParser()
parser.add_argument("--fname", default=None, help = "Base name for output frame file. Otherwise auto-generated ")
parser.add_argument("--background-cache",default=None,help="If nonzero, loads a frame cache, and overlays with result. Used to add noise. Must specify background-channel ")
parser.add_argument("--background-channel",default=None,help="If nonzero, loads a frame cache, and overlays with result. Used to add noise. Must specify background-cache ")
parser.add_argument("--background-add-from-gaussian-psd",action="store_true",help="Adds recolored gaussian noise, based on a specific PSD file")
parser.add_argument("--psd-file",action='append')
parser.add_argument("--instrument", default="H1",help="Use H1, L1,V1")
parser.add_argument("--inj", dest='inj', default=None,help="inspiral XML file containing injection information.")
parser.add_argument("--event",type=int, dest="event_id", default=None,help="event ID of injection XML to use.")
parser.add_argument("--single-ifo",default=False,action='store_true')
parser.add_argument("--approx",type=str,default=None)
parser.add_argument("--srate",type=int,default=16384,help="Sampling rate")
parser.add_argument("--seglen", type=float,default=16., help="Default window size for processing.")
parser.add_argument("--start", type=int,default=None)
parser.add_argument("--stop", type=int,default=None)
parser.add_argument("--fref", dest='fref', type=float, default=0.0, help="Waveform reference frequency [template]. Required, default is 0 (coalescence).")
parser.add_argument("--l-max",default=4,type=int,help='L-max modes in injection')
parser.add_argument("--incl",default=None,help="Set the inclination of L (at fref). Particularly helpful for aligned spin tests")
parser.add_argument("--mass1",default=10,type=float,help='Mass 1 (solar masses)')
parser.add_argument("--mass2",default=1.4,type=float,help='Mass 2 (solar masses)')
parser.add_argument("--verbose", action="store_true",default=False)
parser.add_argument("--use-hlms-as-injections", action="store_true",default=False)
opts=  parser.parse_args()


# Generate signal
P = lalsimutils.ChooseWaveformParams()
P.deltaT = 1./opts.srate
P.radec = True  # use a real source with a real instrument
if not opts.inj:
    P.randomize(aligned_spin_Q=True,default_inclination=opts.incl)
    P.m1 = opts.mass1*lalsimutils.lsu_MSUN
    P.m2 = opts.mass2*lalsimutils.lsu_MSUN
    P.taper = lalsimutils.lsu_TAPER_START
    P.tref =1000000000  # default
    if opts.approx:
        P.approx = opts.approx #lalsim.GetApproximantFromString(str(opts.approx))
    else:
        P.approx = lalsim.GetApproximantFromString("SpinTaylorT2")
else:
    from igwn_ligolw import lsctables, table, utils # check all are needed

    filename = opts.inj
    event = opts.event_id
    xmldoc = utils.load_filename(filename, verbose = True, contenthandler =lalsimutils.cthdler)
    sim_inspiral_table = lsctables.SimInspiralTable.get_table(xmldoc)
    P.copy_sim_inspiral(sim_inspiral_table[int(event)])
    P.taper = lalsimutils.lsu_TAPER_START
    if opts.approx:
        P.approx = opts.approx # lalsim.GetApproximantFromString(str(opts.approx))
#        P.approx = pyseobnr.models.SEOBNRv5HM
#        P.approx = lalsim.GetApproximantFromString(str(SEOBNRv5PHM))
print("Approx =",P.approx)
P.taper = lalsimutils.lsu_TAPER_START  # force taper
P.detector = opts.instrument
if opts.approx == "EccentricTD":
    P.phaseO = 3
P.fref = opts.fref
P.print_params()


T_est = lalsimutils.estimateWaveformDuration(P)
T_est = P.deltaT*lalsimutils.nextPow2(T_est/P.deltaT)
if T_est < opts.seglen:
    T_est =opts.seglen
P.deltaF = 1./T_est
print(" Duration ", T_est)
if T_est < opts.seglen:
    print(" Buffer length too short, automating retuning forced ")


# Generate signal
#hoft = lalsimutils.hoft(P)   # include translation of source, but NOT interpolation onto regular time grid
#Option to use hlms as injection instead of hoft directly (motivation: saw slight differences between hoft and hlmoft for TEOBResumSDALI)
if not (opts.use_hlms_as_injections):
    print("Injecting with hoft")
    hoft = gwsignal.hoft(P,approx_string=opts.approx)
else:
    print("Injecting with hlms")
    hlm = gwsignal.hlmoft(P,Lmax=opts.l_max, approx_string = str(opts.approx))
    hoft = lalsimutils.hoft_from_hlm(hlm,P)
# zero pad to be opts.seglen long, if necessary
if opts.seglen/hoft.deltaT > hoft.data.length:
    TDlenGoal = int(opts.seglen/hoft.deltaT)
    hoft = lal.ResizeREAL8TimeSeries(hoft, 0, TDlenGoal)

# zero pad some more on either side, to make sure the segment covers start to stop
if opts.start and hoft.epoch > opts.start:
    nToAddBefore = int((float(hoft.epoch)-opts.start)/hoft.deltaT)
    # hoft.epoch - nToAddBefore*hoft.deltaT  # this is close to the epoch, but not quite ... we are adjusting it to be within 1 time sample
    print(nToAddBefore, hoft.data.length)
    ht = lal.CreateREAL8TimeSeries("Template h(t)", 
            opts.start , 0, hoft.deltaT, lalsimutils.lsu_DimensionlessUnit, 
            hoft.data.length+nToAddBefore)
    ht.data.data = np.zeros(ht.data.length)  # clear
    ht.data.data[nToAddBefore:nToAddBefore+hoft.data.length] = hoft.data.data
    hoft = ht

if opts.stop and hoft.epoch+hoft.data.length*hoft.deltaT < opts.stop:
    nToAddAtEnd = int( (-(hoft.epoch+hoft.data.length*hoft.deltaT)+opts.stop)/hoft.deltaT)
    print("Padding end ", nToAddAtEnd, hoft.data.length)
    hoft = lal.ResizeREAL8TimeSeries(hoft,0, int(hoft.data.length+nToAddAtEnd))

# Import background data, if needed, and add it
if not ( opts.background_cache is None  or opts.background_channel is None and opts.start is None and opts.stop is None):
    # Don't use specific times
    start = opts.start
    stop = opts.stop
    hoft_bg = lalsimutils.frame_data_to_hoft(opts.background_cache,opts.background_channel,start,stop)
    hoft.data.data += hoft_bg.data.data  # add noise
# elif opts.background_add_from_gaussian_psd:
#     # Retrieve PSDs
#     psd_here=None
#     for inst, psdf in map(lambda c: c.split("="), opts.psd_file):
#         if inst == opts.instrument:
#             print "Reading PSD for instrument %s from %s" % (inst, psdf)
#             psd_here = lalsimutils.get_psd_series_from_xmldoc(psdf, inst)

#     # Generate white noise
#     ht_noise = lal.CreateREAL8TimeSeries("Template h(t)", 
#             hoft.epoch, 0, hoft.deltaT, lalsimutils.lsu_DimensionlessUnit, 
#             hoft.data.length)
#     ht_noise.data.data = np.random.normal(0,1,size=ht_noise.data.length)

#     # Fourier transform
#     hF_noise = lalsimutils.DataFourierREAL8(ht_noise)


channel = opts.instrument+":FAKE-STRAIN"

tstart = int(hoft.epoch)
duration = int(round(hoft.data.length*hoft.deltaT))
if not opts.fname:
    fname = opts.instrument.replace("1","")+"-fake_strain-"+str(tstart)+"-"+str(duration)+".gwf"

print("Writing signal with ", hoft.data.length*hoft.deltaT, " to file ", fname)
lalsimutils.hoft_to_frame_data(fname,channel,hoft)

# TEST: Confirm it works by reading the frame
if opts.verbose:
    print(" -----  Plotting data ------ ")
    import os
    from matplotlib import pyplot as plt
    # First must create corresponding cache file
    os.system("echo "+ fname+ " | lal_path2cache   > test.cache")
    # Now I can read it
    # Beware that the results are OFFSET FROM ONE ANOTHER due to PADDING,
    #    but that the time associations are correct
    hoft2 = lalsimutils.frame_data_to_hoft("test.cache", channel)
    tvals2 = (float(hoft2.epoch) - float(P.tref)) +  np.arange(hoft2.data.length)*hoft2.deltaT
    tvals = (float(hoft.epoch) - float(P.tref)) +  np.arange(hoft.data.length)*hoft.deltaT
    
    fig, ax = plt.subplots(2, 1)
    index_non_zero = np.argwhere(hoft.data.data >0).flatten()[0]
    ax[0].plot(tvals[index_non_zero:index_non_zero + 10000], hoft.data.data[index_non_zero:index_non_zero + 10000])

    ax[1].plot(tvals2,hoft2.data.data,label='Fr')
    ax[1].plot(tvals,hoft.data.data,label='orig')
    ax[1].legend(); #plt.show()
    fig.savefig("injected-data_"+opts.instrument +".png")

    hf = lalsimutils.DataFourierREAL8(hoft2)
    fvals = np.arange(hf.data.length) * hf.deltaF
    fig2, ax2 =  plt.subplots()
    ax2.set_xlim([P.fmin - 2, fvals[-1]])
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.plot(fvals, 2*fvals*np.abs(hf.data.data))
    fig2.savefig("injected-data_"+opts.instrument +"_FD.png")

