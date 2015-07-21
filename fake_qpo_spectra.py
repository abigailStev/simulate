#!/usr/bin/env

import numpy as np
import datetime
import subprocess
import os.path
import warnings
import simulate_lightcurves as simlc
import powerspec as psd
import rebin_powerspec as rb_psd
import ccf as xcf
import plot_ccf as p_xcf
import get_lags as lags

__author__ = 'Abigail Stevens'

################################################################################
def run_fakeit(day, prefix, obs_time, total_spec):
    """
    Makes fake energy spectra for each part of the QPO phase with the HEAsoft
    FTOOL "fakeit".

    Parameters
    ----------
    day : string
        Today's date, in YYMMDD format, for writing to file names.

    :param prefix:
    :param obs_time:
    :param total_spec:

    Returns
    -------
    list of strings
        The file names of the fake spectra created.
    """
    fs_par = [0.07788201, -1.09793882, 0.22469418]  ## Parameters for the sine
            ## wave that best fits the FracSctr variations with QPO phase
    g_par = [0.18449737, -0.56668321, 2.52189912]  ## Parameters for the sine
            ## wave that best fits the Gamma variations with QPO phase

    spec_num = np.arange(total_spec)

    tinybins = np.arange(0, 1.0, 0.01)
    fs_smooth_sine = fs_par[0] * np.sin(2.0 * np.pi * tinybins + fs_par[1]) + fs_par[2]
    g_smooth_sine = g_par[0] * np.sin(2.0 * np.pi * tinybins + g_par[1]) + g_par[2]
    fs_sine_split = np.array_split(fs_smooth_sine, total_spec)
    g_sine_split = np.array_split(g_smooth_sine, total_spec)
    fs_sine = np.asarray([])
    g_sine = np.asarray([])

    for (fs_element, g_element) in zip(fs_sine_split, g_sine_split):
        fs_sine = np.append(fs_sine, np.mean(fs_element))
        g_sine = np.append(g_sine, np.mean(g_element))

    nh_val = 0.6
    simpler_upsc = 1.0
    bb1_kt = 0.719494
    bb1_norm = 3705.67
    bb2_kt = 0.401327
    bb2_norm = 41716.1
    gauss_linee = 6.34799
    gauss_sigma = 1.00000
    gauss_norm = 0.0203765

    rsp_matrix = "/Users/abigailstevens/Reduced_data/GX339-BQPO/PCU2.rsp"
    # bkgd_spec = "/Users/abigailstevens/Reduced_data/GX339-BQPO/evt_bkgd_rebinned.pha"
    model = "phabs*(simpler*bbodyrad+bbodyrad+gaussian)"
    fake_spec_list = []

    for (simpler_gamma, simpler_fs, num_spec) in zip(g_sine, fs_sine, spec_num):
        filename = "/Users/abigailstevens/Dropbox/Research/simulate/spectra/" \
                   "%s_%s_nopoiss_%d" % (prefix, day, num_spec)
        fakeit_script = filename+"_fakeit.xcm"
        fake_spec=filename+".fak"

        mod_par = "& %.1f & %.6f & %.6f & %.1f & %.6f & %.2f & %.6f & %.2f & "\
                "%.6f & %.1f & %.7f" % (nh_val, simpler_gamma, simpler_fs,
                simpler_upsc, bb1_kt, bb1_norm, bb2_kt, bb2_norm, gauss_linee,
                gauss_sigma, gauss_norm)
        # print mod_par

        with open(fakeit_script, mode='w') as out:
            out.write("lmod simpler /Users/abigailstevens/Dropbox/Research/"
                    "xspecmods\n")
            out.write("abund wilm\n")
            out.write("xsect vern\n")
            out.write("mod "+model+"  "+mod_par+"\n")
            out.write("fakeit 1 & %s & & n & & %s & %.1f, 1.0, 1.0 \n" \
                      % (rsp_matrix, fake_spec, obs_time))

        # print fakeit_script
        # print fake_spec
        # subprocess.call(["open", fakeit_script])

        # subprocess.Popen("xspec %s" % (fakeit_script), shell=True)
                ## Setting shell=True allows you to run non- standard shell
                ## commands, and Popen lets us redirect the output

        subprocess.Popen("xspec %s > %s" % (fakeit_script, "dump.txt"), \
                shell=True)  ## Setting shell=True allows you to run non-
                ## standard shell commands, and Popen lets us redirect the
                ## output

        if os.path.exists(fake_spec):
            fake_spec_list.append(fake_spec)
        else:
            warnings.warn("Fake spectrum #%d doesn't exist." % num_spec, \
                    RuntimeWarning)

    return fake_spec_list

################################################################################
def make_powerspectrum(power_array, mean_rate_array, meta_dict, prefix, \
        filename):

    mean_rate = np.mean(mean_rate_array)
    power = np.mean(power_array, axis=1)

    out_file = filename+"_psd.fits"
    plot_file = filename+"_psd_rb.eps"

    freq, power, leahy_power, fracrms_power, fracrms_err = psd.normalize(power,\
            meta_dict, mean_rate, False)

    psd.fits_out(out_file, prefix, meta_dict, mean_rate, freq, fracrms_power, \
            fracrms_err, leahy_power)

    rebin_const = 1.01

    rb_freq, rb_power, rb_err, freq_min, freq_max = rb_psd.geometric_rebinning(\
            freq, fracrms_power, fracrms_err, rebin_const)

    rb_psd.plot_rb(plot_file, rebin_const, prefix, rb_freq, rb_freq*rb_power, \
            rb_freq*rb_err)

    subprocess.call(["open", plot_file])


################################################################################
def make_crosscorrelation(cross_spec_array, ci, ref, meta_dict, prefix, \
        filename):

    ref.mean_rate = np.mean(ref.mean_rate_array[1:])
    ref.power = np.mean(ref.power_array, axis=1)
    ci.power = np.mean(ci.power_array, axis=2)
    print np.shape(ci.power)
    ci.mean_rate = np.mean(ci.mean_rate_array, axis=1)
    print np.shape(ci.mean_rate)
    cross_spec = np.mean(cross_spec_array, axis=2)

    out_file = filename+".fits"
    xcf.save_for_lags(out_file, prefix, meta_dict, ci.mean_rate,
        ref.mean_rate, cross_spec, ci.power, ref.power)

    out_file = filename+"_ccf.fits"
    plot_file = filename+"_ccf.eps"
    t = np.arange(0, meta_dict['n_bins'])

    ccf = xcf.UNFILT_cs_to_ccf(cross_spec, meta_dict, ref, False)
    ccf_error = xcf.standard_ccf_err(cross_spec_array, meta_dict, ref, False)

    xcf.fits_out(out_file, prefix, "None", meta_dict, ci.mean_rate, \
        ref.mean_rate, t, ccf, ccf_error, False)

    time_bins = np.arange(meta_dict['n_bins'])
    pos_time_bins = time_bins[0:meta_dict['n_bins']/2]
    neg_time_bins = time_bins[meta_dict['n_bins']/2:] - meta_dict['n_bins']
    time_bins = np.append(neg_time_bins, pos_time_bins)

    print np.shape(ccf)
    ccf = ccf[:,15]
    pos_time_ccf = ccf[0:meta_dict['n_bins']/2]
    neg_time_ccf = ccf[meta_dict['n_bins']/2:]
    ccf = np.append(neg_time_ccf, pos_time_ccf)

    ccf_error = ccf_error[:,15]
    pos_time_ccf_err = ccf_error[0:meta_dict['n_bins']/2]
    # neg_time_ccf_err = pos_time_ccf_err[::-1]
    neg_time_ccf_err = ccf_error[meta_dict['n_bins']/2:]
    ccf_err = np.append(neg_time_ccf_err, pos_time_ccf_err)

    p_xcf.make_plot(time_bins, ccf, ccf_err, meta_dict['n_bins'], prefix, \
            plot_file, 15, 128)

    subprocess.call(["open", plot_file])


################################################################################
def make_lagspectrum(filename, prefix):
    """
    Makes lag-energy and lag-frequency spectra of the simulated data.

    Parameters
    ---------

    filename : string
        The root of the filename, including full directory path.

    prefix : string
        The identifying prefix of the simulated data.

    Returns
    -------
    nothing

    """
    in_file = filename+"_cs.fits"
    out_file = filename+"_lag.fits"
    lags.main(in_file, out_file, filename, prefix, "eps", 4.0, 7.0, 3, 20)


################################################################################
def main():

    day = datetime.datetime.today().strftime("%y%m%d")
    prefix = "FAKE-GX339B"

    filename = "/Users/abigailstevens/Dropbox/Research/simulate/out_sim/" \
            "%s_%s_nopoiss" % (prefix, day)

    num_seconds = 64
    detchans = 64
    dt_mult = 64
    t_res = 1.0/8192.0
    obs_time = 13056
    num_seg = 204
    # num_seg = 50
    # num_seg = 2
    adjust_seg = 0
    epoch = 5

    dt = dt_mult * t_res
    n_bins = num_seconds * int(1.0 / dt)
    nyquist_freq = 1.0 / (2.0 * dt)
    df = 1.0 / float(num_seconds)

    meta_dict = {'dt': dt, 't_res': t_res, 'num_seconds': num_seconds, \
                 'df': df, 'nyquist': nyquist_freq, 'n_bins': n_bins, \
                 'detchans': detchans, 'adjust_seg': adjust_seg, \
                 'num_seg': num_seg, 'obs_epoch': epoch}

    ref = xcf.Lightcurve()
    ci = xcf.Lightcurve()
    mean_rate_array_1D = 0
    power_array = np.zeros((meta_dict['n_bins'], 1))
    cross_spec_array = np.zeros((meta_dict['n_bins'], meta_dict['detchans'], \
            1), dtype=np.complex128)
    ref.power_array = np.zeros((meta_dict['n_bins'], 1))
    ref.mean_rate_array = 0
    ci.power_array = np.zeros((meta_dict['n_bins'], meta_dict['detchans'], 1))
    ci.mean_rate_array = np.zeros((meta_dict['detchans'], 1))

    # total_spec = 25
    mean = 5.42089871724  ## Hz
    std_dev = 0.352722903769
    frequencies = np.random.normal(loc=mean, scale=std_dev, size=num_seg)
    total_spectra = 1 / frequencies / dt
    total_spectra = total_spectra.astype(int)
    print total_spectra

    ################################
    ## Looping through the segments
    ################################

    for (segment, total_spec) in zip(np.arange(num_seg), total_spectra):

        #############################################
        ## Making fake energy spectra for one period
        #############################################

        fake_spec_list = run_fakeit(day, prefix, obs_time, total_spec)
        spectra = np.zeros((total_spec, meta_dict['detchans']))

        for (spec_file, num_spec) in zip(fake_spec_list, np.arange(total_spec)):
            single_spectrum = simlc.read_fakeit_spectra(spec_file)
            assert len(single_spectrum) == detchans, "ERROR: Channels of fake "\
                    "energy spectrum do not match the expected detector energy"\
                    " channels."
            spectra[num_spec,:] = single_spectrum

        num_repeat = np.ceil(float(meta_dict['n_bins']+num_seg) / \
                float(total_spec))
        spectra = np.tile(spectra, (num_repeat, 1))

        ## Change lightcurve into count RATE units in here by dividing by the
        ## exposure time.
        ## This part I can move up per segment to change it slightly.
        lightcurve_2D = spectra[segment:meta_dict['n_bins']+segment,:]/ obs_time
        lightcurve_1D = np.sum(lightcurve_2D, axis=1)
        lightcurve_ref = xcf.stack_reference_band(lightcurve_2D, \
                meta_dict['obs_epoch'])

        power_segment, mean_rate_segment = psd.make_ps(lightcurve_1D)

        cs_seg, ci_seg, ref_seg = xcf.make_cs(lightcurve_2D, lightcurve_ref, \
                meta_dict)

        ## Adding this segment to the total arrays
        power_array = np.hstack((power_array, np.reshape(power_segment,
                (meta_dict['n_bins'], 1)) ))
        mean_rate_array_1D = np.vstack((mean_rate_array_1D, mean_rate_segment))
        print "Mean count rate per segment:", mean_rate_segment
        cross_spec_array = np.dstack((cross_spec_array, cs_seg))
        ci.power_array = np.dstack((ci.power_array, np.reshape(ci_seg.power, \
                (meta_dict['n_bins'], meta_dict['detchans'], 1)) ))
        ref.power_array = np.hstack((ref.power_array, np.reshape(ref_seg.power,
                (meta_dict['n_bins'], 1)) ))
        ci.mean_rate_array = np.hstack((ci.mean_rate_array,
                np.reshape(ci_seg.mean_rate, (meta_dict['detchans'], 1)) ))
        ref.mean_rate_array = np.vstack((ref.mean_rate_array, \
                ref_seg.mean_rate))

    ## Cutting off the initializing zeros
    mean_rate_array_1D = mean_rate_array_1D[1:]
    power_array = power_array[:,1:]
    ref.power_array = ref.power_array[:,1:]
    ref.mean_rate_array = ref.mean_rate_array[1:]
    ci.power_array = ci.power_array[:,:,1:]
    ci.mean_rate_array = ci.mean_rate_array[:,1:]
    cross_spec_array = cross_spec_array[:,:,1:]

    ###########################
    ## Making a power spectrum
    ###########################

    make_powerspectrum(power_array, mean_rate_array_1D, meta_dict, prefix, \
          filename)

    #######################################
    ## Making a cross correlation function
    #######################################

    make_crosscorrelation(cross_spec_array, ci, ref, meta_dict, prefix, \
            filename)

    ##################################################
    ## Making a lag-energy and lag_frequency spectrum
    ##################################################

    make_lagspectrum(filename, prefix)



################################################################################
if __name__ == "__main__":

    main()
