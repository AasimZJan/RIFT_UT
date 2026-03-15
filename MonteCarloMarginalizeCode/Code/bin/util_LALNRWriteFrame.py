#! /usr/bin/env python
#  
# To generate NR injections using LALSimulations's function. This bypasses RIFT's NR infrastructure, the base code is utiL_LALWriteFrame.py with added function to take in NRhdf5 file and generate polarizations.

import argparse
import numpy as np
import RIFT.lalsimutils as lalsimutils
import lalsimulation as lalsim
import lal
import h5py
from astropy.time import Time
try:
    import sxs
except:
    print('SXS package not installed.')

parser = argparse.ArgumentParser()
parser.add_argument("--fname", default=None, help = "Base name for output frame file. Otherwise auto-generated ")
parser.add_argument("--instrument", default="H1",help="Use H1, L1,V1")
parser.add_argument("--inj", dest='inj', default=None,help="inspiral XML file containing injection information.")
parser.add_argument("--event",type=int, dest="event_id", default=None,help="event ID of injection XML to use.")
parser.add_argument("--approx",type=str,default=None)
parser.add_argument("--srate",type=int,default=16384,help="Sampling rate")
parser.add_argument("--seglen", type=float,default=16., help="Default window size for processing.")
parser.add_argument("--start", type=int,default=None)
parser.add_argument("--stop", type=int,default=None)
parser.add_argument("--incl",default=None,help="Set the inclination of L (at fref). Particularly helpful for aligned spin tests")
parser.add_argument("--mass1",default=10,type=float,help='Mass 1 (solar masses)')
parser.add_argument("--mass2",default=1.4,type=float,help='Mass 2 (solar masses)')
parser.add_argument("--l-max",default=None,type=int,help='Inclusion of modes in injection')
parser.add_argument("--path-to-hdf5", help='Path to NRhdf5 file. This needs to be in the LVK format')
parser.add_argument("--modes-list", type=str, default=None, help="List of specific modes you want to use. Set l-max to None if you want to use this option.")
parser.add_argument("--sxs-simulation-name", type=str, default=None, help="SXS simulation name you want to use to generate the injection")
parser.add_argument("--verbose", action="store_true",default=False)
opts=  parser.parse_args()

def generate_polarizations_from_NRhdf5(P, path_to_hdf5):
    print(f"Reading waveform from {path_to_hdf5}")
    # get mtotal based on user input. This is in kgs.
    mtotal= P.m1 + P.m2 

    # load in hdf5 file to get masses and fmin
    data_1 = h5py.File(path_to_hdf5,"r")
    m1, m2 = data_1.attrs["mass1"] * mtotal, data_1.attrs["mass2"] * mtotal
    fmin = data_1.attrs["f_lower_at_1MSUN"] * lal.MSUN_SI/mtotal
    fref = 0.0 # set to zero to avoid errors
    print(f"Smallest possible frequency for this waveform {fmin} Hz. Frequency at 1 solar mass is {data_1.attrs['f_lower_at_1MSUN']}.\nReference frequency is set to {fref} Hz since the lalsimulation function does not take non-zero reference frequency.")

    # if provided fmin is lower than the waveform can actually have
    if P.fmin < fmin and P.fmin != 0.0:
        fmin = fmin + 0.5*10**(-2)*fmin
        print(f"WARNING: Can't have fmin less than that of the NR waveform. Provided fmin is {P.fmin} Hz, defaulting to fmin={fmin} Hz.")
    else:
        fmin = P.fmin
        print(f"Generating waveform with fmin is {P.fmin} Hz.")

    # get spins, useful for precessing case
    s1x, s1y, s1z, s2x, s2y, s2z = lalsim.SimInspiralNRWaveformGetSpinsFromHDF5File(fref, mtotal/lal.MSUN_SI, path_to_hdf5)

    # extract modes, either using l-max or modes-list option
    params = lal.CreateDict()
    modes = []
    lmax = opts.l_max
    if opts.modes_list == None and lmax is not None:
        modes = [(l, m) for l in range(2, lmax + 1) for m in range(-l, l + 1) if m != 0]
    elif opts.modes_list is not None and lmax is None:
        modes = list(eval(opts.modes_list))
    elif opts.modes_list is not None and lmax is not None:
        raise ValueError("Use either l_max or modes_list, not both.")
    else:
        raise ValueError("One of l_max or modes_list must be provided.")
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
                P.phiref, P.psi, P.eccentricity, P.meanPerAno, P.deltaT, fmin, fref, params, lalsim.NR_hdf5)

    # find epoch: based on the apporach in GWSignal.py
    amplitude = np.sqrt(h_p.data.data**2 + h_c.data.data**2)
    max_amp_index = np.argmax(amplitude)
    print(f"Initial epoch = {h_p.epoch}, updated epoch = {-max_amp_index * h_p.deltaT}")
    h_p.epoch, h_c.epoch = -max_amp_index * h_p.deltaT, -max_amp_index * h_c.deltaT

    return h_p, h_c

def get_lvk_modes_from_SXS(P, simulation_name, modes_list=opts.modes_list, l_max=opts.l_max, set_fref_equal_to_flow=True):
    print(f"Loading SXS waveform {simulation_name}")
    # Load metadata
    sim = sxs.load(simulation_name)

    # Converstion factors using lal
    kgs_to_sec = lal.G_SI/lal.C_SI**3
    code_units_to_sec = lal.MTSUN_SI 
    meters_to_sec = 1/lal.C_SI

    # Conversion of mass and distance to seconds
    mtotal = P.m1 + P.m2 # in kgs
    mtot_in_sec = mtotal * kgs_to_sec
    dist_in_sec = P.dist * meters_to_sec # distance in m

    # Get mass in real units based on simulation parameters
    m1, m2 = sim.metadata['reference_mass1'] * mtotal, sim.metadata['reference_mass2'] * mtotal 

    # Get minimum frequency. Note: initial_orbital_frequency gives frequency that is almost of f_low
    # This did not match with NRHybSur3dq8 comparisons
    # _, _, dyn_test = sim.to_lvk(t_ref=0.0, ell_max=8)
    # flow_at_1MSUN = dyn_test['f_low'] / code_units_to_sec 
    flow_at_1MSUN = sim.metadata['initial_orbital_frequency']/ (np.pi*code_units_to_sec)
    flow = flow_at_1MSUN * lal.MSUN_SI/mtotal 
    
    # Get fref. f_22 = |omega|/pi
    Omega_vec = sim.metadata["reference_orbital_frequency"]
    Omega = np.linalg.norm(Omega_vec)
    fref_at_1MSUN = np.array(Omega) / (np.pi*code_units_to_sec)
    fref = fref_at_1MSUN * lal.MSUN_SI/mtotal
    print(f"Smallest possible frequency for this waveform {flow} Hz. Frequency at 1 solar mass is {flow_at_1MSUN} Hz.\nReference frequency for this waveform is {fref} Hz. Reference frequency at 1 solar mass in {fref_at_1MSUN} Hz.")

    # Sanity checks for fmin. Don't go below what the waveform can provide.
    if P.fmin <= flow and P.fmin != 0.0:
        fmin = flow #+ 0.5*10**(-2)*flow
        print(f"WARNING: Can't have fmin less than that of the NR waveform. Provided fmin is {P.fmin} Hz, defaulting to fmin={fmin} Hz.")
    else:
        fmin = P.fmin
        print(f"Generating waveform with fmin is {P.fmin} Hz.")
    print(f"Generating waveform with m1 = {m1/lal.MSUN_SI:0.4f} MSUN, m2 = {m2/lal.MSUN_SI:0.4f} MSUN \n s1 = {sim.metadata['reference_dimensionless_spin1']}, s2 = {sim.metadata['reference_dimensionless_spin2']}, eccentricity = {sim.metadata['reference_eccentricity']}, meanAnomaly = {sim.metadata['reference_mean_anomaly']}\n fmin = {fmin} Hz")

    # Collect modes based on input
    modes = []
    if modes_list == None and l_max is not None:
        modes = [(l, m) for l in range(2, l_max + 1) for m in range(-l, l + 1) if m != 0]
    elif modes_list is not None and l_max is None:
        modes = list(eval(modes_list))
    elif modes_list is not None and l_max is not None:
        raise ValueError("Use either l_max or modes_list, not both.")
    else:
        raise ValueError("One of l_max or modes_list must be provided.")
    print(f"Modes used = {modes}")
    
    # Now we use dt and fmin to generate modes
    dt_in_code_units = P.deltaT / code_units_to_sec * lal.MSUN_SI/mtotal
    print(f"WARNING: The provided fmin has no effect; the waveform starts at the lowest available frequency of the NR simulation, which is {flow} Hz.")
    fmin_in_code_units = fmin * (mtotal/lal.MSUN_SI) * code_units_to_sec
    fref_in_code_units = fref * (mtotal/lal.MSUN_SI) * code_units_to_sec
    # This function gives weird results for eccentric waveform if f_ref is different than f_low
    if set_fref_equal_to_flow:
        fref_in_code_units = fmin_in_code_units
    times, hlms_dict, dyn = sim.to_lvk(f_ref=fref_in_code_units, ell_max=l_max, dt=dt_in_code_units, f_low=fmin_in_code_units)

    # collect modes
    hlm = {}
    # remove junk radiation, which I have found to be usually within 150M of the start.
    tvals = np.arange(len(times))*dt_in_code_units
    index = np.argwhere(tvals >= 150).flatten()[0]
    print(f"Removing content from 0 M upto {tvals[index]} M")
    for i, mode in enumerate(modes):
        # extract modes from the SXS object
        mode_content_here = hlms_dict[mode][index:]

        # Scale based on distance and mtotal
        mode_content_at_distance = mtot_in_sec/dist_in_sec * mode_content_here

        # Save as a lal object
        wf = lal.CreateCOMPLEX16TimeSeries("hlm", 0, 0, P.deltaT, lal.DimensionlessUnit, len(mode_content_at_distance))
        wf.data.data = mode_content_at_distance

        # resize
        if P.deltaF:
            TDlen = int(1./P.deltaF * 1./P.deltaT)
            if TDlen < wf.data.length:   # Truncate the series to the desired length, removing data at the *start* (left)
                wf = lal.ResizeCOMPLEX16TimeSeries(wf, wf.data.length-TDlen, TDlen)
            elif TDlen > wf.data.length:   # Zero pad, extend at end
                wf = lal.ResizeCOMPLEX16TimeSeries(wf, 0, TDlen)
        
        # tapering
        taper = True
        if taper:
            ntaper = int(0.01*TDlen)
            if P.fmin > 0: # avoid failure if waveform start frequency 0 is nominally specified
                ntaper = np.max([ntaper, int(1./(P.fmin*P.deltaT))])
            vectaper= 0.5 - 0.5*np.cos(np.pi*np.arange(ntaper)/(1.*ntaper))
            # Taper at the start of the segment
            wf.data.data[:ntaper]*=vectaper
        hlm[modes[i][0],modes[i][1]] = wf
    
    # set epoch based on GWsignal approach
    rhosq = np.zeros(TDlen)
    for mode in hlm:
        rhosq += np.abs(hlm[mode].data.data)**2
    indx_max = np.argmax(rhosq)
    new_epoch = - indx_max*P.deltaT
    for mode in hlm:
        hlm[mode].epoch = new_epoch
    return hlm
    

def generate_hoft(P, path_to_hdf5, Fp=None, Fc=None):

    P_copy = P.manual_copy()
    hp, hc = generate_polarizations_from_NRhdf5(P_copy, path_to_hdf5) 

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
    #if opts.approx:
    #    P.approx = lalsim.GetApproximantFromString(str(opts.approx))

P.taper = lalsimutils.lsu_TAPER_START  # force taper
P.detector = opts.instrument
if opts.approx == "EccentricTD":
    P.phaseO = 3
P.print_params()

T_est = lalsimutils.estimateWaveformDuration(P)
T_est = P.deltaT*lalsimutils.nextPow2(T_est/P.deltaT)
if T_est > opts.seglen:
    print(" WARNING: THE SIGNAL WILL LIKELY BE TRUNCATED when writing the frame, which is VERY BAD ")
T_est = opts.seglen
P.deltaF = 1./T_est
print(" Duration ", T_est)
if T_est < opts.seglen:
    print(" Buffer length too short, automating retuning forced ")

# Generate signal
if opts.sxs_simulation_name is not None:
    hlms = get_lvk_modes_from_SXS(P, opts.sxs_simulation_name)
    hoft = lalsimutils.hoft_from_hlm(hlms, P)
else:
    hoft = generate_hoft(P, opts.path_to_hdf5)   # include translation of source, but NOT interpolation onto regular time grid
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
