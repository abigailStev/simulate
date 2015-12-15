#!/usr/bin/env python
"""
Makes fake QPO energy-dependent light curves and does spectral-timing (power
spectrum, ccf, and lag-energy spectra) from a set of SED parameters (either tied
value of parameters of best-fit function to the parameter variation with QPO-
phase), from energy_spectra/multifit_plots.py.

Notes: HEASOFT 6.11.* and Python 2.7.* (with supporting libraries) must be
installed in order to run this script.

WARNING: Location of 'simpler' XSPEC model is hardwired in.

Example call:
    python fake_qpo_spectra.py ./cygx1_BBPL ./cygx1_BBPL_sines.txt \
            --prefix cygx1 --rsp ./cygx1_PCU2.rsp

"""
import numpy as np
import subprocess
import os
import warnings
import argparse
from astropy.io import fits

## These are things I've written.
## Their directories are in my PYTHONPATH bash environment variable.
import simulate_lightcurves as simlc  ## in simulate
import powerspec as psd  ## in power_spectra
import rebin_powerspec as rb_psd  ## in power_spectra
import ccf as xcor  ## in cross_correlation
import plot_ccf as p_xcor  ## in cross_correlation
import plot_2d as p2D_xcor  ## in cross_correlation
import get_lags as lags  ## in lag_spectra
import tools  ## in whizzy_scripts
import ccf_lightcurves as ccf_lc  ## in cross_correlation

__author__ = 'Abigail Stevens <A.L.Stevens at uva.nl>'
__year__ = "2015"

class Parameter(object):
    """
    This Parameter class is slightly different from the one in multfit_plots.py.
    This one doesn't have any strings for model name or plot label, and requires
    n_spectra.
    """
    def __init__(self, n_spectra):
        self.value = np.zeros(n_spectra)
        self.error = 0
        self.lo_v = 0
        self.hi_v = 0
        self.pos_err = 0
        self.neg_err = 0
        self.varying = False
        self.funcfit = np.zeros(n_spectra)
        self.best_fit = np.zeros(5)
        self.phase = None
        self.phase_err = None


################################################################################
def fit_function(t, p):
    """
    Computing a function to fit to the SED parameter variations. This is exactly
    the same as what's in energy_spectra/multifit_plots.py.

    Parameters
    ----------
    t : np.array of floats
        1-D array of time steps for the fit function.

    p : np.array of floats
        1-D array of the function parameters.

    Returns
    -------
    np.array of floats
        1-D array of the function fit to the data at steps t with parameters p.
    """
    return p[0] * np.sin(2. * np.pi * t + p[1]) + \
           p[2] * np.sin(2. * 2. * np.pi * t + p[3]) + p[4]


################################################################################
def average_of_parameter(temp_stack):
    """
    Takes the average of a 2-D array of values along the 0 axis. Excludes the
    0th row from this average because reasons. (because the arrays are defined
    with zeros going down the zeroth row, and those get removed later once stuff
    is appended to that)

    Parameters
    ----------
    temp_stack : np.array of floats
        2-D array of parameter values.

    Returns
    -------
    np.array of floats
        1-D array of averaged parameter values.
    """
    # print "Before:", np.shape(temp_stack)
    temp_stack = np.average(temp_stack[1:,:], axis=0)
    # print "After:", np.shape(temp_stack)
    # print temp_stack

    return temp_stack


################################################################################
def read_parameters(funcfit_file, n_spectra, n_params):
    """
    Reads in the fit function parameters that were put out by multifit_plots.py.

    Parameters
    ----------
    funcfit_file : str
        The full path of the file with the best-fitting function parameters.

    n_spectra : int
        The number of QPO phase-resolved SEDs that were simultaneously fit,
        i.e. the number of x-axis points for the fit function.

    n_params : int
        The number of SED parameters for all models for one data set.

    Returns
    -------
    np.array of Parameter objects
        1-D array (size=n_params) of all the SED parameters for one data set.

    str
        The SED model as fed into XSPEC, with no spaces.
    """

    tinybins = np.arange(0, 1.0, 0.01)
    model = []
    parameters = np.array([Parameter(n_spectra) for i in range(n_params)])
    # print np.shape(parameters)

    with open(funcfit_file, 'r') as f:

        for line in f:
            line_element = line.strip().split("    ")
            model.append(line_element[0])
            assert line_element[0] == model[0]

            for (element, i) in zip(line_element[1:], range(n_params)):
                if "[" in element or "]" in element or "," in element:
                    parameters[i].varying = True
                    # print "Best fit par before:", parameters[i].best_fit
                    temp = element.replace('[','').replace(']','').split(',')
                    # print "Temp:", temp
                    parameters[i].best_fit = np.vstack((parameters[i].best_fit,
                            np.array(temp, dtype=np.float)))
                    # print "Best fit par after:", parameters[i].best_fit
                else:
                    temp = float(element)
                    parameters[i].value = np.vstack((parameters[i].value,
                            np.repeat(temp, n_spectra) ))

                    # print "Value shape:", np.shape(parameters[k].value)

    # print np.shape(parameters)
    # for par in parameters:
    #     print par.value
    #     print par.best_fit
    #
    # print np.shape(parameters)

    for component in parameters:
        if component.varying:
            component.best_fit = average_of_parameter(component.best_fit)

            smooth_func = fit_function(tinybins, component.best_fit)
            func_split = np.array_split(smooth_func, n_spectra)
            component.value = np.array([np.mean(element) for element in \
                    func_split])
            # print component.best_fit
        else:
            component.value = average_of_parameter(component.value)
            # print component.value

    print("SED model: %s" % model[0])
    return parameters, model[0]


################################################################################
def run_fakeit(parameters, model, exposure, n_spectra, out_root, rsp_matrix):
    """
    Makes fake SED for each time bin of the QPO phase with the HEAsoft FTOOL
    "fakeit".

    Parameters
    ----------
    parameters : list of Parameter objects
        List of SED parameters.

    model : str
        The full SED model to fake.

    exposure : float
        Total exposure time of the observation, in seconds.

    n_spectra : int
        Total number of SEDs per QPO phase to make/use.

    out_root : str
        Root of name of the output file (i.e., without extension).

    Returns
    -------
    fake_spec_list : list of strings
        The file names of the fake SEDs created.
    """

    fake_spec_list = []
    current_dir = os.getcwd()
    # print "Current dir:", current_dir
    out_dir = os.path.dirname(out_root)
    # print out_dir
    os.chdir(out_dir)
    out_file = os.path.basename(out_root)

    for spec_num in np.arange(n_spectra):
        # subprocess.call(["cd", out_dir])

        # out_file = "%s_%d" % (out_root, spec_num)
        fakeit_script = "%s_%d_fakeit.xcm" % (out_root, spec_num)
        fake_spec = "%s_%d.fak" % (out_file, spec_num)  ## Need local name,
                ## otherwise get "***XSPEC Error: Invalid chars in file name."
        mod_par = ""
        for component in parameters:
            mod_par += "& %.6g " % (component.value[spec_num])

        with open(fakeit_script, mode='w') as out:
            out.write("lmod simpler /Users/abigailstevens/Dropbox/Research/"
                    "xspecmods\n")
            out.write("query yes\n")
            out.write("abund wilm\n")
            out.write("xsect vern\n")
            out.write("mod "+model+"  "+mod_par+"\n")
            out.write("fakeit 1 & %s & & N & & %s & %.3f, 1.0, 1.0 \n" \
                      % (rsp_matrix, fake_spec, exposure))

        # subprocess.call(["open", fakeit_script])
        # subprocess.call(["open", "dump.txt"])
        # print os.getcwd()
        # print "Script:", os.path.exists(fakeit_script)

        # p = subprocess.Popen("xspec %s" % (fakeit_script), shell=True)
            ## Setting shell=True allows you to run non- standard shell
            ## commands, and Popen lets us redirect the output
        p = subprocess.Popen("xspec %s > %s" % (fakeit_script, "dump.txt"), \
                shell=True)  ## Setting shell=True allows you to run non-
                ## standard shell commands, and Popen lets us redirect the
                ## output
        p.communicate()  ## Waits for the command to finish running.


        # print os.getcwd()
        # subprocess.call(["open", "dump.txt"])
        # print "Spectrum:", os.path.exists("%s/%s" % (out_dir, fake_spec))

        if os.path.exists("%s/%s" % (out_dir, fake_spec)):
            tools.replace_key_val(fake_spec, 1, "POISSERR", False)
            fake_spec_list.append("%s/%s" % (out_dir, fake_spec))
            # print tools.get_key_val(fake_spec, 1, "POISSERR")
        else:
            warnings.warn("Fake SED #%d doesn't exist. Exiting." % \
                    spec_num, RuntimeWarning)
            exit()

    os.chdir(current_dir)

    return fake_spec_list


################################################################################
def make_powerspectrum(power_array, mean_rate_array, meta_dict, prefix,
        out_root):
    """
    Makes a power spectrum of the simulated data and plots it.

    Parameters
    ----------
    power_array : np.array of floats
        2-D array (size = n_bins x n_seg) of the un-normalized power spectrum.

    mean_rate_array : np.array of floats
        1-D array (size = n_seg) of the mean count rate.

    meta_dict : dict
        Dictionary of important parameters.

    prefix : str
        Identifying prefix of the simulated data (should start with 'FAKE').

    out_root : str
        Root of the output file name, including full directory path.

    Returns
    -------
    nothing

    Files created
    -------------
    *_psd.fits :
        Data file of the segment-averaged power spectrum (PSD) of all the data.

    *_psd_rb.eps :
        Plot of the geometrically re-binned averaged power spectrum.

    """
    # print np.shape(power_array)
    mean_rate = np.mean(mean_rate_array)
    power = np.mean(power_array, axis=-1)

    out_file = out_root + "_psd.fits"
    plot_file = out_root + "_psd_rb.eps"

    freq, power, leahy_power, fracrms_power, fracrms_err = psd.normalize(power,
            meta_dict, mean_rate, False)

    psd.fits_out(out_file, prefix, meta_dict, mean_rate, freq, fracrms_power,
            fracrms_err, leahy_power, "Simulated PSD")

    rebin_const = 1.01

    rb_freq, rb_power, rb_err, freq_min, freq_max = rb_psd.geometric_rebinning(\
            freq, fracrms_power, fracrms_err, rebin_const)

    rb_psd.plot_rb(plot_file, rebin_const, prefix, rb_freq, rb_freq*rb_power,
            rb_freq*rb_err)

    subprocess.call(["open", plot_file])


################################################################################
def make_crosscorrelation(cross_spec_array, ci, ref, meta_dict, prefix,
        out_root):
    """
    Makes the ccf of the simulated data and plots it.

    Parameters
    ---------
    cross_spec_array : np.array of complex numbers
        3-D array (size = n_bins x detchans x n_seg) of the cross spectrum.

    ci : ccf_lc.Lightcurve object
        Channel of interest.

    ref : ccf_lc.Lightcurve object
        Reference band.

    meta_dict : dict
        Dictionary of important parameters.

    prefix : string
        Identifying prefix of the simulated data (should start with 'FAKE').

    out_root : string
        Root of the output file name, including full directory path.

    Returns
    -------
    nothing

    Files created
    -------------
    *_ccf.fits :
        Data file of the segment-averaged cross-correlation function.

    *_ccf.eps :
        Plot of the CCF in energy channel 15.

    *_ccf_2D.eps :
        Plot of the 2-dimensional CCF.

    """

    ci.mean_rate /= np.float(meta_dict['n_seg'])
    ci.power /= np.float(meta_dict['n_seg'])
    ref.power /= np.float(meta_dict['n_seg'])
    ref.mean_rate /= np.float(meta_dict['n_seg'])
    avg_cross_spec = np.mean(cross_spec_array, axis=-1)
    ci.pos_power = ci.power[0:meta_dict['n_bins']/2+1, :]
    ref.pos_power = ref.power[0:meta_dict['n_bins']/2+1]

    ## Compute the variance and rms of the absolute-rms-normalized reference
    ## band power spectrum. Remember that this has no Poisson noise!
    absrms_ref_pow = xcor.raw_to_absrms(ref.pos_power, ref.mean_rate,
            meta_dict['n_bins'], meta_dict['dt'], noisy=False)

    ref.var, ref.rms = xcor.var_and_rms(absrms_ref_pow, meta_dict['df'])

    ## Save the cross spectrum and power spectra to compute lags
    out_file = out_root + ".fits"
    xcor.save_for_lags(out_file, prefix, meta_dict, avg_cross_spec, ci, ref)

    ## Make ccf output file names
    out_file = out_root + "_ccf.fits"
    plot_file = out_root + "_ccf.eps"

    meta_dict['dt'] = np.repeat(meta_dict['dt'], meta_dict['n_seg'])

    ## Compute average ccf and error on average ccf
    ccf_avg, ccf_error = xcor.unfilt_cs_to_ccf_w_err(cross_spec_array,
            meta_dict, ref)

    ## Save ccf to a fits file
    print "Rms float:", float(ref.rms)
    print "Rms:", ref.rms
    print type(ref.rms)
    xcor.fits_out(out_file, prefix, " ", meta_dict, ci.mean_rate, ref.mean_rate,
            float(ref.rms), ccf_avg, ccf_error, -1, -1, "Simulated CCF")

    ## Make a 1-D ccf plot for energy channel 15
    time_bins = np.arange(meta_dict['n_bins'])
    half = meta_dict['n_bins']/2
    pos_time_bins = time_bins[0:half]
    neg_time_bins = time_bins[half:] - meta_dict['n_bins']
    time_bins = np.append(neg_time_bins, pos_time_bins)
    ccf_15 = ccf_avg[:,15]
    pos_time_ccf = ccf_15[0:half]
    neg_time_ccf = ccf_15[half:]
    ccf_15 = np.append(neg_time_ccf, pos_time_ccf)
    ccf_error = ccf_error[:,15]
    pos_time_ccf_err = ccf_error[0:half]
    neg_time_ccf_err = ccf_error[half:]
    ccf_err = np.append(neg_time_ccf_err, pos_time_ccf_err)

    p_xcor.make_plot(time_bins, ccf_15, ccf_err, meta_dict['n_bins'], prefix,
            plot_file, 15, int(1.0 / np.mean(meta_dict['dt'])))

    ## Make a 2-D ccf plot
    t_length = 30
    ccf = ccf_avg[half - t_length:half + t_length + 1, :].T
    t_bins = time_bins[half - t_length:half + t_length + 1]
    # energies = np.loadtxt("/Users/abigailstevens/Dropbox/Academic/Conferences_and_Talks/DC_talks/NICER-energies.dat")
    energies = np.loadtxt("/Users/abigailstevens/Reduced_data/" + prefix +
            "/energies.txt")
    plot_file = out_root + "ccf_2D.eps"

    print np.shape(ccf)
    print np.shape(t_bins)
    p2D_xcor.make_plot(ccf, t_bins, ci.mean_rate, t_length,
            int(1.0 / np.mean(meta_dict['dt'])), plot_file, energies)

    subprocess.call(["open", plot_file])


################################################################################
def chisquared_lagspectrum(data_file, sim_file):
    """
    Computes the chisquared statistic of the simulated lag-energy spectrum with
    the data's lag-energy spectrum.

    S. Vaughan, "Scientific Inference: Learning from data" 1st ed, p134 eq 6.12

    Parameters
    ----------
    data_file : string
        Name of the FITS file with the lags of the data.

    sim_file : string
        Name of the FITS file with the lags of the simulation.

    Returns
    -------
    chisquared : float
        The chisquared fit statistic of the simulation's lag-energy spectrum
        with the data's lag-energy spectrum.

    dof : float
        The relevant degrees of freedom (number of data points + number of
        simulation data points).

    """

    try:
        data_hdu = fits.open(data_file)
    except IOError:
        print "\tERROR: File does not exist: %s" % data_file
        exit()

    try:
        sim_hdu = fits.open(sim_file)
    except IOError:
        print "\tERROR: File does not exist: %s" % sim_file
        exit()

    data_lag = data_hdu[2].data.field('TIME_LAG')
    data_err = data_hdu[2].data.field('TIME_LAG_ERR')
    sim_lag = sim_hdu[2].data.field('TIME_LAG')

    data_lag = np.delete(data_lag, 10)
    data_err = np.delete(data_err, 10)
    sim_lag = np.delete(sim_lag, 10)

    ## Only fitting between 3-20 keV (incl)
    data_lag = data_lag[2:25]
    data_err = data_err[2:25]
    sim_lag = sim_lag[2:25]

    resids = np.square(data_lag - sim_lag) / np.square(data_err)
    chisquared = np.sum(resids)

    return chisquared, len(data_lag) + len(sim_lag)


################################################################################
def make_lagspectrum(out_root, prefix):
    """
    Makes lag-energy and lag-frequency spectra of the simulated data, and
    computes the chisquared fit statistic of the simulation with the data.

    WARNING: data file name is hard-coded!

    Parameters
    ---------

    out_root : str
        The root of the out_root, including full directory path.

    prefix : str
        Identifying prefix of the simulated data (should start with 'FAKE').

    Returns
    -------
    nothing

    """
    in_file = out_root + "_cs.fits"
    out_file = out_root + "_lag.fits"

    energies_file = "/Users/abigailstevens/Reduced_data/GX339-BQPO/energies.txt"

    lags.main(in_file, out_file, energies_file, out_root, prefix, "eps", 4.0, \
            7.0, 3, 20)
    subprocess.call(["open", out_root + "_lag-energy.eps"])

    data_file = "/Users/abigailstevens/Dropbox/Research/lag_spectra/out_lags/"\
            "GX339-BQPO/GX339-BQPO_151204_t64_64sec_adj_lag.fits"

    chisquared, dof = chisquared_lagspectrum(data_file, out_file)
    print chisquared, dof

    ensp_model_file = "/Users/abigailstevens/Dropbox/Research/CCF_paper1/ensp_"\
            "models.txt"

    with open(ensp_model_file, 'a') as out:
        out.write("$%.3f / %d$  & \\ ref{ } \\\\ \n" % (chisquared, dof))


################################################################################
def main(out_root, funcfit_file, prefix="GX339-BQPO", n_bins=8192, dt=0.0078125,
        n_seg=198, n_chans=64, exposure=13224.3984375, n_spectra=24, n_params=9,
        epoch=5, rsp_matrix="PCU2.rsp", test=False, psd_flag=False):
    """
    Main of fake_qpo_spectra.py

    Parameters
    ----------
    out_root : str
        The full path of the file to write the FITS output to.

    funcfit_file : str
        The full path of the text file containing the fit function parameters
        of the SED parameter variations.

    prefix : str
        The identifying prefix of the data. Should have 'FAKE' in it for
        simulated data. [GX339-BQPO]

    n_bins : int
        Number of time bins in a segment of light curve. [8192]

    dt : float
        Timestep between time bins. [0.0078125]

    n_seg : int
        Number of segments to compute (i.e., length of light curve). [198]

    n_chans : int
        Number of effective energy channels for the detector's data mode. [64]

    exposure : float
        The exposure time of the fake observation. [13224.3984375]

    n_spectra : int
        Number of phase-resolved SEDs per QPO phase. [24]

    n_params : int
        Total number of parameters in one SED. [9]

    epoch : int
        RXTE observation epoch (affects the keV energy of the detector channels,
        which affects which ones we bin up for the reference band, and the
        conversion from channel to keV for plotting). [5]

    rsp_matrix : str
        The file name (local path) of the XSPEC response matrix. ["PCU2.rsp"]

    test : bool
        If true, runs one segment for testing. [False]

    psd_flag : bool
        If true, computes and prints the power spectral density of the fake QPO.
        [False]

    Returns
    -------
    nothing
    """

    if test:
        n_seg = 10

    adjust_seg = 0  ## No need to adjust perfectly simulated data.
    n_seconds = n_bins * dt
    nyquist_freq = 1.0 / (2.0 * dt)
    df = 1.0 / float(n_seconds)

    meta_dict = {'dt': dt,
                 'n_seconds': n_seconds,
                 'df': df,
                 'nyquist': nyquist_freq,
                 'n_bins': n_bins,
                 'detchans': n_chans,
                 'adjust_seg': adjust_seg,
                 'n_seg': n_seg,
                 'obs_epoch': epoch,
                 'exposure': exposure,
                 'filter': False}
    # print(meta_dict)
    print_seg = 20

    ############################################################
    ## Reading in SED parameters from the text file
    ############################################################
    parameters, model = read_parameters(funcfit_file, n_spectra, n_params)

    #############################################
    ## Making fake SED for one period
    #############################################
    fake_spec_list = run_fakeit(parameters, model, meta_dict['exposure'], \
            n_spectra, out_root, rsp_matrix)

    spectra = np.zeros((n_spectra, meta_dict['detchans']))

    for (spec_file, n_spec) in zip(fake_spec_list, np.arange(n_spectra)):

        single_spectrum = simlc.read_fakeit_spectra(spec_file)
        assert len(single_spectrum) == n_chans, "ERROR: Channels of fake "\
                "SED do not match the expected detector energy"\
                " channels."
        spectra[n_spec, :] = single_spectrum

    n_repeat = np.ceil(float(meta_dict['n_bins'] + meta_dict['n_seg']) / \
            float(n_spectra))
    spectra = np.tile(spectra, (n_repeat, 1))

    #############################################
    ## Looping through the segments to do timing
    #############################################
    ref = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ref')
    ci = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ci')
    mean_rate_array_1D = 0
    power_array = np.zeros((meta_dict['n_bins'], 1))
    cross_spec_array = np.zeros((meta_dict['n_bins'], meta_dict['detchans'], \
            1), dtype=np.complex128)

    # mean = 5.42089871724  ## Hz
    # std_dev = 0.352722903769
    # frequencies = np.random.normal(loc=mean, scale=std_dev, size=n_seg)
    # n_spectra = 1 / frequencies / dt
    # n_spectra = n_spectra.astype(int)
    # print n_spectra

    print("Segments computed:")
    # for (segment, n_spectra) in zip(np.arange(n_seg), n_spectra):
    for segment in np.arange(n_seg):

        if segment % print_seg == 0:
            print("\t%d" % segment)

        ## Change lightcurve into count RATE units in here by dividing by the
        ## exposure time.
        ## This part I can move up per segment to change it slightly.
        lightcurve_2D = spectra[segment:meta_dict['n_bins']+segment,:] / \
                meta_dict['exposure']
        lightcurve_1D = np.sum(lightcurve_2D, axis=1)
        lightcurve_ref = xcor.stack_reference_band(lightcurve_2D,
                instrument='PCA', obs_epoch=meta_dict['obs_epoch'])

        power_segment, mean_rate_segment = psd.make_ps(lightcurve_1D)

        cs_seg, ci_seg, ref_seg = xcor.make_cs(lightcurve_2D, lightcurve_ref, \
                meta_dict)

        ## Adding this segment to the total arrays
        power_array = np.hstack((power_array, np.reshape(power_segment,
                (meta_dict['n_bins'], 1)) ))
        mean_rate_array_1D = np.vstack((mean_rate_array_1D, mean_rate_segment))
        # # print "Mean count rate per segment:", mean_rate_segment
        # cross_spec_array = np.dstack((cross_spec_array, cs_seg))
        # ci.power_array = np.dstack((ci.power_array, np.reshape(ci_seg.power, \
        #         (meta_dict['n_bins'], meta_dict['detchans'], 1)) ))
        # ref.power_array = np.hstack((ref.power_array, np.reshape(ref_seg.power,
        #         (meta_dict['n_bins'], 1)) ))
        # ci.mean_rate_array = np.hstack((ci.mean_rate_array,
        #         np.reshape(ci_seg.mean_rate, (meta_dict['detchans'], 1)) ))
        # ref.mean_rate_array = np.vstack((ref.mean_rate_array, \
        #         ref_seg.mean_rate))

        ## Compute the variance and rms of each segment's reference band power
        ## spectrum. Remember that there's no Poisson noise, so noisy=False!
        absrms_pow = xcor.raw_to_absrms(ref_seg.power[0:meta_dict['n_bins'] \
                / 2 + 1], ref_seg.mean_rate, meta_dict['n_bins'],
                meta_dict['dt'], noisy=False)

        var, rms = xcor.var_and_rms(absrms_pow, meta_dict['df'])

        ## Append segment to arrays
        cross_spec_array = np.dstack((cross_spec_array, cs_seg))
        ci.mean_rate_array = np.hstack((ci.mean_rate_array,
                np.reshape(ci_seg.mean_rate, (meta_dict['detchans'],
                1))))
        ref.power_array = np.hstack((ref.power_array,
                np.reshape(ref_seg.power, (meta_dict['n_bins'], 1))))
        ref.mean_rate_array = np.append(ref.mean_rate_array,
                ref_seg.mean_rate)
        ref.var_array = np.append(ref.var_array, var)

        ## Sum across segments -- arrays, so it adds by index
        ci.mean_rate += ci_seg.mean_rate
        ref.mean_rate += ref_seg.mean_rate
        ci.power += ci_seg.power
        ref.power += ref_seg.power

    print("Number of segments computed: %d" % n_seg)

    ## Cutting off the initializing zeros
    mean_rate_array_1D = mean_rate_array_1D[1:]
    power_array = power_array[:,1:]
    cross_spec_array = cross_spec_array[:,:,1:]
    ci.mean_rate_array = ci.mean_rate_array[:,1:]
    ref.power_array = ref.power_array[:,1:]
    ref.mean_rate_array = ref.mean_rate_array[1:]
    ref.var_array = ref.var_array[1:]

    print ref.var_array

    # print mean_rate_array_1D[1:4]
    # print ref.mean_rate_array[1:4]

    ###########################
    ## Making a power spectrum
    ###########################
    if psd_flag:
        make_powerspectrum(power_array, mean_rate_array_1D, meta_dict, prefix,
                out_root)

    #######################################
    ## Making a cross correlation function
    #######################################

    make_crosscorrelation(cross_spec_array, ci, ref, meta_dict, prefix,
            out_root)

    ##################################################
    ## Making a lag-energy and lag_frequency spectrum
    ##################################################

    make_lagspectrum(out_root, prefix)


################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python fake_qpo_spectra.py outroot"\
            " funcfitfile [OTHER ARGUMENTS]", description=__doc__, \
            epilog="For optional arguments, default values are given in "\
            "brackets at end of description.")

    parser.add_argument('outroot', help="The root of the file names to write "\
            "the different outputs to.")

    parser.add_argument('funcfitfile', help="Full path of the text file with "\
            "the best-fitting function paramters of the varying SED "\
            "parameters.")

    parser.add_argument('--prefix', dest='prefix', help="Identifying prefix of"\
            " the simulation. Should start with 'FAKE'. [FAKE]")

    parser.add_argument('-n', '--n_bins', type=tools.type_power_of_two,
            default=8192, dest='n_bins', help="Number of time bins in each "\
            "Fourier segment. Must be a power of 2, positive, integer. [8192]")

    parser.add_argument('--dt', type=tools.type_positive_float,
            default=0.0078125,
            dest='dt', help="Timestep between bins, in seconds. [0.0078125]")

    parser.add_argument('-g', '--n_seg', type=tools.type_positive_int,
            default=198, dest='n_seg', help="Number of segments to compute. "\
            "[198]")

    parser.add_argument('-s', '--n_spec', type=tools.type_positive_int,
            default=24, dest='n_spectra', help="Number of SED per QPO phase. "\
            "[24]")

    parser.add_argument('--n_par', type=tools.type_positive_int, default=9,
            dest='n_params', help="Total number of spectral parameters per "\
            "SED. [9]")

    parser.add_argument('-c', '--chan', type=tools.type_positive_int,
            default=64, dest='n_chans', help="Number of detector energy "\
            "channels for the data mode used for real data. [64]")

    parser.add_argument('-e', '--epoch', type=tools.type_positive_int,
            default=5,
            choices={1,2,3,4,5}, dest='epoch', help="RXTE observation epoch. "\
            "[5]")

    parser.add_argument('-x', '--exposure', type=tools.type_positive_float,
            default=13224.0, \
            dest='exposure', help="Exposure time of the simulated observation,"\
            " in seconds. [13224.0]")

    parser.add_argument('--rsp', default="PCU2.rsp", dest='rsp_matrix',
            help="The file name (local path) of the XSPEC response matrix. "\
            "[./PCU2.rsp]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
            dest='test', help="Int flag: 0 if computing all segments, 1 if "\
            "computing only one segment for testing. [0]")

    parser.add_argument('--power', action='store_true', default=False, \
            dest='psd_flag', help="If present, saves and plots the "\
            "power spectrum of the fake QPO. [False]")

    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True

    main(args.outroot, args.funcfitfile, prefix=args.prefix, n_bins=args.n_bins,
            dt=args.dt, n_seg=args.n_seg, n_chans=args.n_chans,
            exposure=args.exposure, n_spectra=args.n_spectra,
            n_params=args.n_params, epoch=args.epoch,
            rsp_matrix=args.rsp_matrix, test=test, psd_flag=args.psd_flag)
