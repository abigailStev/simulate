import numpy as np
import subprocess
import os
import warnings
import argparse
import simulate_lightcurves as simlc
import powerspec as psd
import rebin_powerspec as rb_psd
import ccf as xcf
import plot_ccf as p_xcf
import get_lags as lags
import tools

__author__ = 'Abigail Stevens <A.L.Stevens at uva.nl>'

"""
Makes legit fake QPO and does spectral timing from a set of energy spectral
parameters (either tied value of parameters of best-fit sine wave to QPO-phase
variations), from energy_spectra/multifit_plots.py.
"""
class Parameter(object):
    """
    This Parameter class is slightly different from the one in multfit_plots.py.
    This one doesn't have any strings for model name or plot label, and requires
     num_spectra.
    """
    def __init__(self, num_spectra):
        self.value = np.zeros(num_spectra)
        self.error = 0
        self.lo_v = 0
        self.hi_v = 0
        self.pos_err = 0
        self.neg_err = 0
        self.varying = False
        self.sinefit = np.zeros(num_spectra)
        self.best_fit = np.zeros(3)  ## three parameters used to define a sine
        self.phase = None
        self.phase_err = None


################################################################################
def average_of_parameter(temp_stack):
    # print np.shape(temp_stack)
    temp_stack = np.average(temp_stack[1:,:], axis=0)
    # print np.shape(temp_stack)
    # print temp_stack

    return temp_stack

################################################################################
def read_parameters(sinefit_file, num_spectra, num_params):

    tinybins = np.arange(0, 1.0, 0.01)
    model = []
    parameters = np.array([Parameter(num_spectra) for i in range(num_params)])
    # print np.shape(parameters)

    with open(sinefit_file, 'r') as f:

        for line in f:
            line_element = line.strip().split("    ")
            model.append(line_element[0])
            assert line_element[0] == model[0]

            for (element, i) in zip(line_element[1:], range(num_params)):
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
                            np.repeat(temp, num_spectra) ))

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
            smooth_sine = component.best_fit[0] * np.sin(2.0 * np.pi * \
                    tinybins + component.best_fit[1]) + component.best_fit[2]
            sine_split = np.array_split(smooth_sine, num_spectra)
            component.value = np.array([np.mean(element) for element in \
                    sine_split])
            print component.best_fit
        else:
            component.value = average_of_parameter(component.value)
            print component.value

    print("Model: %s" % model[0])
    return parameters, model[0]

################################################################################
def run_fakeit(parameters, model, exposure, num_spectra, out_root):
    """
    Makes fake energy spectra for each part of the QPO phase with the HEAsoft
    FTOOL "fakeit".

    Parameters
    ----------
    parameters : list
        List of energy spectral parameters.

    model : str
        The full energy spectral model to fake.

    exposure : float
        Total exposure time of the observation, in seconds.

    num_spectra : int
        Total number of spectra per QPO phase to make/use.

    out_root : str
        Root of name of the output file (i.e., without extension).

    Returns
    -------
    list of strings
        The file names of the fake spectra created.
    """

    rsp_matrix = "GX339-BQPO_PCU2.rsp"
    fake_spec_list = []
    current_dir = os.getcwd()
    # print "Current dir:", current_dir
    out_dir = os.path.dirname(out_root)
    # print out_dir
    os.chdir(out_dir)
    out_file = os.path.basename(out_root)

    for spec_num in np.arange(num_spectra):
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
            warnings.warn("Fake spectrum #%d doesn't exist. Exiting." % \
                    spec_num, RuntimeWarning)
            exit()

    os.chdir(current_dir)

    return fake_spec_list


################################################################################
def make_powerspectrum(power_array, mean_rate_array, meta_dict, prefix, \
        out_root):
    """
    Makes a power spectrum of the simulated data and plots it.

    Parameters
    ----------
    power_array : np.array of floats (n_bins x num_seg)
        2-D array of the unnormalized power spectrum.

    mean_rate_array : np.array of floats (num_seg)
        1-D array of the mean count rate.

    meta_dict : dict
        Dictionary of important parameters.

    prefix : string
        Identifying prefix of the simulated data (should start with 'FAKE').

    out_root : string
        Root of the output file name, including full directory path.

    Returns
    -------
    nothing

    """
    print np.shape(power_array)
    mean_rate = np.mean(mean_rate_array)
    power = xcf.seg_average(power_array)

    out_file = out_root+"_psd.fits"
    plot_file = out_root+"_psd_rb.eps"

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
        out_root):
    """
    Makes the ccf of the simulated data and plots it.

    Parameters
    ---------
    cross_spec_array : np.array of complex numbers (n_bins x detchans x num_seg)
        3-D array of the cross spectrum.

    ci : Lightcurve object
        Channel of interest.

    ref : Lightcurve object
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

    """
    ref.mean_rate = xcf.seg_average(ref.mean_rate_array)[0]
    ref.power = xcf.seg_average(ref.power_array)
    ci.power = xcf.seg_average(ci.power_array)
    # print np.shape(ci.power)
    ci.mean_rate = xcf.seg_average(ci.mean_rate_array)
    # print np.shape(ci.mean_rate)
    cross_spec = xcf.seg_average(cross_spec_array)

    out_file = out_root+".fits"
    xcf.save_for_lags(out_file, prefix, meta_dict, ci.mean_rate,
        ref.mean_rate, cross_spec, ci.power, ref.power)

    out_file = out_root+"_ccf.fits"
    plot_file = out_root+"_ccf.eps"
    t = np.arange(0, meta_dict['n_bins'])

    meta_dict['dt'] = np.repeat(meta_dict['dt'], meta_dict['num_seg'])

    ccf = xcf.UNFILT_cs_to_ccf(cross_spec, meta_dict, ref, False)
    ccf_error = xcf.standard_ccf_err(cross_spec_array, meta_dict, ref, False)

    xcf.fits_out(out_file, prefix, "None", meta_dict, ci.mean_rate, \
        ref.mean_rate, t, ccf, ccf_error, False, -1, -1)

    time_bins = np.arange(meta_dict['n_bins'])
    pos_time_bins = time_bins[0:meta_dict['n_bins']/2]
    neg_time_bins = time_bins[meta_dict['n_bins']/2:] - meta_dict['n_bins']
    time_bins = np.append(neg_time_bins, pos_time_bins)

    # print np.shape(ccf)
    ccf = ccf[:,15]
    pos_time_ccf = ccf[0:meta_dict['n_bins']/2]
    neg_time_ccf = ccf[meta_dict['n_bins']/2:]
    ccf = np.append(neg_time_ccf, pos_time_ccf)

    ccf_error = ccf_error[:,15]
    pos_time_ccf_err = ccf_error[0:meta_dict['n_bins']/2]
    neg_time_ccf_err = ccf_error[meta_dict['n_bins']/2:]
    ccf_err = np.append(neg_time_ccf_err, pos_time_ccf_err)

    p_xcf.make_plot(time_bins, ccf, ccf_err, meta_dict['n_bins'], prefix, \
            plot_file, 15, 128)

    subprocess.call(["open", plot_file])


################################################################################
def make_lagspectrum(out_root, prefix):
    """
    Makes lag-energy and lag-frequency spectra of the simulated data.

    Parameters
    ---------

    out_root : string
        The root of the out_root, including full directory path.

    prefix : string
        Identifying prefix of the simulated data (should start with 'FAKE').

    Returns
    -------
    nothing

    """
    in_file = out_root+"_cs.fits"
    out_file = out_root+"_lag.fits"
    lags.main(in_file, out_file, out_root, prefix, "eps", 4.0, 7.0, 3, 20)
    subprocess.call(["open", out_root + "_lag-energy.eps"])


################################################################################
def main(out_root, sinefit_file, prefix, n_bins, dt, num_seg, detchans,\
        exposure, num_spectra, num_params, epoch, test, psd_flag):
    """
    Main of fake_qpo_spectra.py

    Parameters
    ----------
    out_root : str
        Description.

    sinefit_file : str
        Description.

    prefix : str
        The identifying prefix of the data. Should have 'FAKE' in it for
        simulated data.

    num_seconds : float
        Number of seconds per segment of lightcurve.

    dt : float
        Timestep between time bins.

    num_seg : int
        Number of segments to compute (i.e., length of lightcurve).

    exposure : float
        The exposure time of the fake observation.

    num_spectra : int
        Number of phase-resolved energy spectra per QPO phase.

    num_params : int
        Total number of spectral parameters in one energy spectrum.

    epoch : int
        RXTE observation epoch.

    test : bool
        If true, runs one segment for testing.

    psd_flag : bool
        If true, computes and prints the power density spectrum of the fake QPO.

    Returns
    -------
    nothing
    """

    if test:
        num_seg = 1

    adjust_seg = 0
    num_seconds = n_bins * dt
    nyquist_freq = 1.0 / (2.0 * dt)
    df = 1.0 / float(num_seconds)

    meta_dict = {'dt': dt, 'num_seconds': num_seconds, 'df': df, \
            'nyquist': nyquist_freq, 'n_bins': n_bins, 'detchans': detchans, \
            'adjust_seg': adjust_seg, 'num_seg': num_seg, 'obs_epoch': epoch,
            'exposure': exposure}
    print_seg = 20
    # print(meta_dict)

    ############################################################
    ## Reading in energy spectral parameters from the text file
    ############################################################
    parameters, model = read_parameters(sinefit_file, num_spectra, num_params)

    #############################################
    ## Making fake energy spectra for one period
    #############################################

    fake_spec_list = run_fakeit(parameters, model, meta_dict['exposure'], \
            num_spectra, out_root)

    spectra = np.zeros((num_spectra, meta_dict['detchans']))

    for (spec_file, num_spec) in zip(fake_spec_list, np.arange(num_spectra)):

        single_spectrum = simlc.read_fakeit_spectra(spec_file)
        assert len(single_spectrum) == detchans, "ERROR: Channels of fake "\
                "energy spectrum do not match the expected detector energy"\
                " channels."
        spectra[num_spec,:] = single_spectrum

    num_repeat = np.ceil(float(meta_dict['n_bins']+meta_dict['num_seg']) / \
            float(num_spectra))
    spectra = np.tile(spectra, (num_repeat, 1))

    #############################################
    ## Looping through the segments to do timing
    #############################################
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

    # mean = 5.42089871724  ## Hz
    # std_dev = 0.352722903769
    # frequencies = np.random.normal(loc=mean, scale=std_dev, size=num_seg)
    # num_spectra = 1 / frequencies / dt
    # num_spectra = num_spectra.astype(int)
    # print num_spectra

    print("Segments computed:")
    # for (segment, num_spectra) in zip(np.arange(num_seg), num_spectra):
    for segment in np.arange(num_seg):

        if num_seg % print_seg == 0:
            print("\t%d" % num_seg)

        ## Change lightcurve into count RATE units in here by dividing by the
        ## exposure time.
        ## This part I can move up per segment to change it slightly.
        lightcurve_2D = spectra[segment:meta_dict['n_bins']+segment,:] / \
                meta_dict['exposure']
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
        # print "Mean count rate per segment:", mean_rate_segment
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

    # print mean_rate_array_1D[1:4]
    # print ref.mean_rate_array[1:4]

    ###########################
    ## Making a power spectrum
    ###########################
    if psd_flag:
        make_powerspectrum(power_array, mean_rate_array_1D, meta_dict, prefix, \
                out_root)

    #######################################
    ## Making a cross correlation function
    #######################################

    make_crosscorrelation(cross_spec_array, ci, ref, meta_dict, prefix, \
            out_root)

    print "place 2"

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
            " sinefitfile [OTHER ARGUMENTS]", \
            description="Computes a fake QPO power spectrum, ccf, and lag-"\
            "energy spectrum given the parameters of the best-fit sine wave to "\
            "the energy spectral parameter variations with QPO phase.", \
            epilog="For optional arguments, default values are given in "\
            "brackets at end of description.")

    parser.add_argument('outroot', help="The root of the file names to write "\
            "the different outputs to.")

    parser.add_argument('sinefitfile', help="Full path of the text file with "\
            "the best-fit sine wave paramters of the varying energy spectral "\
            "parameters.")

    parser.add_argument('--prefix', dest='prefix', help="Identifying prefix of"\
            " the simulation. Should start with 'FAKE'. [FAKE]")

    parser.add_argument('-n', '--n_bins', type=tools.type_power_of_two,
            default=1, dest='n_bins', help="Number of time bins in each "\
            "Fourier segment. Must be a power of 2, positive, integer. [1]")

    parser.add_argument('--dt', type=tools.type_positive_float, default=1, dest='dt',
            help="Timestep between bins, in seconds. [1]")

    parser.add_argument('-g', '--n_seg', type=tools.type_positive_int,
            default=1, dest='num_seg', help="Number of segments to compute. [1]")

    parser.add_argument('-s', '--n_spec', type=tools.type_positive_int, default=1,
            dest='num_spectra', help="Number of energy spectra per QPO phase. "\
            "[1]")

    parser.add_argument('--n_par', type=tools.type_positive_int, default=1,
            dest='num_params', help="Total number of spectral parameters per "\
            "energy spectrum.")

    parser.add_argument('-c', '--chan', type=tools.type_positive_int,
            default=64, dest='detchans', help="Number of detector energy "\
            "channels for the data mode used for real data. [64]")

    parser.add_argument('-e', '--epoch', type=tools.type_positive_int, default=5,
            choices={1,2,3,4,5}, dest='epoch', help="RXTE observation epoch. "\
            "[5]")

    parser.add_argument('-x', '--exposure', type=tools.type_positive_float, default=1, \
            dest='exposure', help="Exposure time of the simulated observation,"\
            " in seconds. [1]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
            dest='test', help="Int flag: 0 if computing all segments, 1 if "\
            "computing only one segment for testing. [0]")

    parser.add_argument('--power', action='store_true', default=False, \
            dest='psd_flag', help="If present, saves and plots the "\
            "power spectrum of the fake QPO. [False]")

    args = parser.parse_args()

    main(args.outroot, args.sinefitfile, args.prefix, args.n_bins, args.dt, \
            args.num_seg, args.detchans, args.exposure, args.num_spectra, \
            args.num_params, args.epoch, args.test, args.psd_flag)
