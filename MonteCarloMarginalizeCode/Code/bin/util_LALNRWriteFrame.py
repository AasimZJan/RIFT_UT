#! /usr/bin/env python
#  
# To generate NR injections using LALSimulations's function. This bypasses RIFT's NR infrastructure, the base code is utiL_LALWriteFrame.py with added function to take in NRhdf5 file and generate polarizations.


import argparse
import numpy as np
import RIFT.lalsimutils as lalsimutils
import lalsimulation as lalsim
import lalframe
import lal

parser = argparse.ArgumentParser()
parser.add_argument("--fname", default=None, help = "Base name for output frame file. Otherwise auto-generated ")
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
parser.add_argument("--incl",default=None,help="Set the inclination of L (at fref). Particularly helpful for aligned spin tests")
parser.add_argument("--mass1",default=10,type=float,help='Mass 1 (solar masses)')
parser.add_argument("--mass2",default=1.4,type=float,help='Mass 2 (solar masses)')
parser.add_argument("--l-max",default=None,type=float,help='Inclusion of modes in injection')
parser.add_argument("--path-to-NRhdf5", help='Path to NRhdf5 file. This needs to be in the LVK format')
parser.add_argument("--modes-list", default=None, help="List of specific modes you want to use. Set l-max to None if you want to use this option.")
parser.add_argument("--verbose", action="store_true",default=False)
opts=  parser.parse_args()

def generate_polarizations_from_NRhdf5(P, opts.path_to_NRhdf5):

    print(f"Reading waveform from {opts.path_to_NRhdf5}")

    # get mtotal based on user input. This is in kgs.
    mtotal= (P.m1 + P.m2) 

    # load in hdf5 file to get masses, mass ratio and fmin
    data_1 = h5py.File(path_to_hdf5,"r")
    m1 = data_1.attrs["mass1"] * mtotal
    m2 = data_1.attrs["mass2"] * mtotal
    fmin = data_1.attrs["f_lower_at_1MSUN"] * lal.MSUN_SI/mtotal
    fref = 0.0 # set to zero to avoid errors
    print(f"Smallest possible fmin for this waveform {fmin} Hz. fmin at 1 solar mass is {data_1.attrs['f_lower_at_1MSUN']}. fref is set to 0.0")

    # if provided fmin is lower than the waveform can actually have
    if P.fmin < fmin and P.fmin != 0.0:
        fmin = fmin + 0.5*10**(-2)*fmin
        print(f"WARNING: Can't have fmin less than that of the NR waveform. Provided fmin is {P.fmin} Hz, defaulting to fmin={fmin} Hz.")
    else:
        fmin = P.fmin
        print(f"Generating waveform with fmin is {P.fmin} Hz.")

    # get spins, useful for precessing case
    s1x, s1y, s1z, s2x, s2y, s2z = lalsim.SimInspiralNRWaveformGetSpinsFromHDF5File(fref, mtotal/lal.MSUN_SI, opts.path_to_hdf5)

    # extract modes, either using l-max or modes-list option
    params = lal.CreateDict()
    modes = []
    lmax = opts.l_max
    if opts.modes_list == None and lmax is not None:
        for l in range(2,lmax+1):
            for m in range(-l,0):
                modes.append((l,m))
            for m in range(1,l+1):
                modes.append((l,m))
    elif opts.modes_list is not None and lmax is None:
        for j in only_mode:
            modes.append(j)
    elif only_mode is not None and lmax is not None:
        print("Inconsistent input, use either l-max or modes-list.")
        sys.exit()
    print(f"Modes used = {modes}")
    ma = lalsim.SimInspiralCreateModeArray()
    for l,m in modes:
        lalsim.SimInspiralModeArrayActivateMode(ma, l, m)
    
    # pass modes and hdf5 path to lalDict object 
    lalsim.SimInspiralWaveformParamsInsertModeArray(params, ma)
    lalsim.SimInspiralWaveformParamsInsertNumRelData(params, path_to_hdf5)

    # sanity print statement
    print(f"Generating waveform with m1 = {m1/lal.MSUN_SI:0.4f} MSUN, m2 = {m2/lal.MSUN_SI:0.4f} MSUN \n s1 = {s1x, s1y, s1z}, s2 = {s2x, s2y, s2z}\n fmin = {fmin} Hz, fref= {fref}")

    # generate polarizations
    h_p, h_c = lalsim.SimInspiralChooseTDWaveform(m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, P.dist, P.incl, \
                P.phiref, P.psi, P.eccentricity, P.meanPerAno, P.deltaT, fmin, fref, params, lalsim.NR_hdf5 )

    # find epoch: based on the apporach in GWSignal.py
    amplitude = np.sqrt(h_p.data.data**2 + h_c.data.data**2)
    max_amp_index = np.argmax(amplitude)
    print(f"Initial epoch = {h_p.epoch}, updated epoch = {-max_amp_index * h_p.deltaT}")
    h_p.epoch, h_c.epoch = -max_amp_index * h_p.deltaT, -max_amp_index * h_c.deltaT

    return h_p, h_c


def hoft(P, opts.path_to_NRhdf5):

    P_copy = P.manual_copy()
    hp, hc = generate_polarizations_from_NRhdf5(P_copy, opts.path_to_NRhdf5) 

    # Apply detector response
    if Fp!=None and Fc!=None:
        hp.data.data *= Fp
        hc.data.data *= Fc
        hp = lal.AddREAL8TimeSeries(hp, hc)
        ht = hp
    elif P.radec==False:
        fp = Fplus(P.theta, P.phi, P.psi)
        fc = Fcross(P.theta, P.phi, P.psi)
        hp.data.data *= fp
        hc.data.data *= fc
        hp = lal.AddREAL8TimeSeries(hp, hc)
        ht = hp
    else:
        # If astropy Time function, overwrite with GPS time, otherwise use normal addition
        if isinstance(hp.epoch, Time):
            dT = hp.epoch.to_value('gps','long')  # pull out the time
            hp.epoch = P.tref + dT
            hc.epoch = P.tref +dT
        else:
            hp.epoch = hp.epoch + P.tref
            hc.epoch = hc.epoch + P.tref
        ht = lalsim.SimDetectorStrainREAL8TimeSeries(hp, hc,
                P.phi, P.theta, P.psi,
                lalsim.DetectorPrefixToLALDetector(str(P.detector)))

    # Resize such that TDlen = 1/deltaF
    if P.deltaF is not None:
        TDlen = int(1./P.deltaF * 1./P.deltaT)
        assert TDlen >= ht.data.length, f"TDlen = {TDlen}, data_length = {ht.data.length}, 1/deltaT = {1/P.deltaT}, 1/deltaF = {1/P.deltaF}"
        if TDlen < ht.data.length:
            print(f'Data removed from {TDlen}:{ht.data.length}.') 
        ht = lal.ResizeREAL8TimeSeries(ht, 0, TDlen)

    # Match lalsimutils tapering
    try:
        taper=True 
        if taper :
            ntaper = int(0.01*TDlen)
            if P.fmin > 0: # avoid failure if waveform start frequency 0 is nominally specified
                ntaper = np.max([ntaper, int(1./(P.fmin*P.deltaT))])
            vectaper= 0.5 - 0.5*np.cos(np.pi*np.arange(ntaper)/(1.*ntaper))
            # Taper at the start of the segment
            ht.data.data[:ntaper]*=vectaper
    except Exception as e:
        print("Couldn't apply tapering", e)
    return ht

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
        P.approx = lalsim.GetApproximantFromString(str(opts.approx))
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
        P.approx = lalsim.GetApproximantFromString(str(opts.approx))

P.taper = lalsimutils.lsu_TAPER_START  # force taper
P.detector = opts.instrument
if opts.approx == "EccentricTD":
    P.phaseO = 3
P.print_params()

T_est = lalsimutils.estimateWaveformDuration(P)
T_est = P.deltaT*lalsimutils.nextPow2(T_est/P.deltaT)
if T_est > opts.seglen:
    print(" WARNING: THE SIGNAL WILL LIKELY BE TRUNCATED when writing the frame, which is VERY BAD ")
T_est =opts.seglen
P.deltaF = 1./T_est
print(" Duration ", T_est)
if T_est < opts.seglen:
    print(" Buffer length too short, automating retuning forced ")

# Generate signal
hoft = hoft(P, opts.path_to_hdf5)   # include translation of source, but NOT interpolation onto regular time grid
epoch_orig = hoft.epoch
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
    plt.plot(tvals2,hoft2.data.data,label='Fr')
    plt.plot(tvals,hoft.data.data,label='orig')
    plt.xlim(float(epoch_orig)- float(P.tref), 0.2)
    plt.xlabel('t - tref')
    plt.legend(); #plt.show()
    plt.savefig("injected-data_"+opts.instrument +".png")
