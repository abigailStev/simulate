#!/usr/bin/env python
"""
Makes fake lag-energy spectra given the energy spectral parameters from
multifit_plots.py in sed_fit_bootstrap.sh. Calls fake_qpo_spectra.py.
Used in ccf_bootstrap.sh.

Be sure that the directories point correctly for your setup.

Example call:  python sim_qpo_bootstrap.py --prefix FAKE-cygx1

"""
__author__ = 'Abigail Stevens <A.L.Stevens at uva.nl>'
__year__ = "2015"

import argparse
import numpy as np
import os.path
import tools  ## in https://github.com/abigailStev/whizzy_scripts
import subprocess
import fake_qpo_spectra
from datetime import datetime
import multifit_plots
import overplot_lags
from fake_qpo_spectra import chisquared_lagspectrum as lag_chisq
import sed_pars
from astropy.io import fits
from fake_qpo_spectra import fit_function
from fake_qpo_spectra import average_of_parameter

HOME_DIR = os.path.expanduser("~")
CCF_BASEDIR = HOME_DIR + "/Dropbox/Research/cross_correlation"
ES_BASEDIR = HOME_DIR + "/Dropbox/Research/energy_spectra"
LAG_BASEDIR = HOME_DIR + "/Dropbox/Research/lags"
EXE_DIR = HOME_DIR + "/Dropbox/Research/simulate"
DATA_BASEDIR = HOME_DIR + "/Reduced_data"
DUMP_FILE = "dump.txt"  # Name of dumping file for intermediary steps

N_PARAMS = int(os.getenv('n_params', 9))
N_SPECTRA = int(os.getenv('n_spectra', 24))
FIT_SPECIFIER = os.getenv('fit_specifier', "1BB-FS-G-Tin-fzs-fzNbb")

DYLD_LIB_PATH = os.environ.get('DYLD_LIBRARY_PATH')
HEADAS_DIR = os.environ.get('HEADAS')
HEAINIT = ". " + HEADAS_DIR + "/headas-init.sh"

PLOT_EXT = "eps"

def plot_final_multifit(var_pars, out_rootname="out", \
        parfit_file="funcfit.fits"):
    """
    Plots the SED parameter variations of the original run with the best-fit
    function averaged over all the bootstrapped runs.

    Parameters
    ----------
    out_rootname : str
        Full path root name (directories plus base) of the output file. Used for
        the data file and the plot file.

    parfit_file : str
        Full path name of the

    Returns
    -------
    nothing
    """
    data_file = out_rootname + ".fits"
    plot_file = out_rootname + "." + PLOT_EXT
    varpar_fits_file = parfit_file.replace("_funcfit.txt", "_varpars.txt")


    print plot_file, parfit_file, data_file

    multifit_plots.make_var_plots(plot_file, N_SPECTRA, var_pars, quiet=False)


def plot_final_lagenergy(lag_data_file, out_rootname="out"):
    """
    Plots the lag-energy of this simulation with the data.

    TODO: gives the chisquared fit (from fake_qpo_spectra.chisquared_lagspectrum).

    Parameters
    ----------
    lag_data_file : str
        The lag-energy spectrum from the original data that we are fitting to.

    out_rootname : str
        The dir+basename of the output file that the simulated lag-energy
        spectrum was written to, and to write the lag-energy plot to.

    Returns
    -------
    nothing
    """

    in_file_list = [out_rootname + "_lag.fits", lag_data_file]
    labels = ["Simulation", "Data"]
    overplot_lags.make_plot(in_file_list, labels, out_rootname)


def avg_boot_var_pars(var_pars, boot_num=2):
    """

    :param var_pars:
    :param boot_num:
    :return:
    """
    n_varpar = np.shape(var_pars)[1]
    # print n_varpar

    avg_var_pars = np.asarray([])

    for par in range(n_varpar):
        parameter_values = np.zeros(N_SPECTRA)
        for boot in range(boot_num):
            # print var_pars[boot, par].value
            parameter_values = np.vstack((parameter_values, \
                    var_pars[boot, par].value))

        # print np.shape(parameter_values)
        parameter_values = parameter_values[1:]
        # print np.shape(parameter_values)
        variance = np.var(parameter_values, axis=0, ddof=1)
        err = np.sqrt(variance)

        avg_par = sed_pars.Parameter(mod_name=var_pars[0,par].mod_name,
                label=var_pars[0,par].label, par_name=var_pars[0,par].par_name)
        avg_par.value = np.mean(parameter_values, axis=0)
        avg_par.pos_err = err
        avg_par.neg_err = err
        avg_par.varying = True
        avg_var_pars = np.append(avg_var_pars, avg_par)

    return avg_var_pars



def read_all_log_files(es_dir, out_rel_name, boot_num):
    """

    :param es_dir:
    :param out_rel_name:
    :param boot_num:
    :return:
    """
    sed_parameters = np.tile(sed_pars.Parameter(), N_PARAMS)
    var_pars = []

    for i in range(1, boot_num + 1):
        log_file = es_dir + "/" + out_rel_name + "_b-" + str(i) + "_xspec.log"
        print log_file

        boot_sed_parameters, n_spec_log = multifit_plots.read_log_file(log_file,
                quiet=False)
        assert n_spec_log == N_SPECTRA
        sed_parameters = np.vstack((sed_parameters, boot_sed_parameters))
        boot_var_par = multifit_plots.determine_varying_parameters(boot_sed_parameters, \
            n_spectra=N_SPECTRA, quiet=True)
        var_pars.append(boot_var_par)

    sed_parameters = sed_parameters[1:]
    # print sed_parameters.__str__()
    # print var_pars
    var_pars = np.asarray(var_pars)
    # print np.shape(var_pars)
    # print var_pars.__str__()
    return var_pars


def avg_fit_func(var_pars, parfit_file):
    """
    Gets the average of the fit function parameters from each bootstrap
    iteration.

    Parameters
    ----------
    var_pars : np.array of sed_pars.Parameter objects
        1-D array of the SED parameters that vary with QPO phase.

    parfit_file : str
        Full path of the "_funcfit.txt" file that was made in
        energy_spectra/multifit_plots.py and read in by ./fake_qpo_spectra.py.

    Returns
    -------
    var_pars : np.array of sed_pars.Parameter objects
        Same as input, but with best_fit, func_fit, phase, and phase_err
        assigned.

    """
    tinybins = np.arange(0, 1.0, 0.01)

    for par in var_pars:
        par.best_fit = np.zeros(5)

    with open(parfit_file, 'r') as f:
        for line in f:
            i = 0
            line_element = line.strip().split("    ")
            for element in line_element[1:]:
                if "[" in element or "]" in element or "," in element:
                    # print "Best fit par before:", var_pars[i].best_fit
                    temp = element.replace('[', '').replace(']', '').split(',')
                    # print "Temp:", temp
                    var_pars[i].best_fit = np.vstack((var_pars[i].best_fit,
                            np.array(temp, dtype=np.float)))
                    var_pars[i].phase = np.append(var_pars[i].phase,
                            var_pars[i].best_fit[1] / (2.0 * np.pi))
                    # print "Best fit par after:", var_pars[i].best_fit
                    i += 1

    for par in var_pars:
        print "\t", par.mod_name, par.par_name
        # print np.shape(par.phase)
        par.best_fit = np.mean(par.best_fit[1:,:], axis=0)
        par.funcfit = fake_qpo_spectra.fit_function(np.arange(-0.02, 1.02,
                0.01), par.best_fit)
        print "Predicted phase:", par.best_fit[1] / (2.0 * np.pi)
        par.phase_err = np.sqrt(np.var(par.phase[1:], ddof=1))
        par.phase = np.mean(par.phase[1:])
        print "Average phase:", par.phase
        print "Err on phase:", par.phase_err

    return var_pars


def main(prefix="--", dt_mult=64, n_seconds=64, boot_num=2, testing=False,
        day=datetime.now().strftime("%y%m%d")):
    """
    Main of sim_qpo_bootstrap.py python script.

    Parameters
    ----------
    prefix : str
        The identifying prefix of the data (object nickname or data ID). [--]

    dt_mult : float
        Multiple of dt (dt is from data file) for timestep between bins.
        Must be a power of 2, positive. [64]

    n_seconds : int
        Number of seconds in each Fourier segment. Must be a power of 2,
        positive. [64]

    boot_num : int
        Number of bootstrap realizations to do/done. Entering 0 means don't
        bootstrap. [2]

    testing : bool
        If present, sets to True. Used for testing script so it'll only run one
        segment. [False]

    day : str
        Date of the analysis, in YYMMDD format. [default is today's date]

    Returns
    -------
    nothing
    """
    ccf_dir = CCF_BASEDIR + "/out_ccf/" + prefix
    es_dir = ES_BASEDIR + "/out_es/" + prefix + "/bootstrapped"
    lag_dir = LAG_BASEDIR + "/out_lags/" + prefix
    out_dir = EXE_DIR + "/out_sim/" + prefix + "/bootstrapped"
    data_dir = DATA_BASEDIR + "/" + prefix
    rsp_matrix = prefix + "_PCU2.rsp"
    parfit_file = es_dir + "/" + prefix + "_" + day + "_" + FIT_SPECIFIER + \
            "_funcfit.txt"
    ccf_data_file = ccf_dir + "/" + prefix + "_" + day + "_t" + str(dt_mult) + \
            "_" + str(n_seconds) + "sec_adj.fits"
    lag_data_file = lag_dir + "/" + prefix + "_" + day + "_t" + str(dt_mult) + \
            "_" +str(n_seconds) + "sec_adj_lag.fits"
    epoch = 5

    out_rel_name = prefix + "_" + day + "_" + FIT_SPECIFIER
    out_rootname = out_dir + "/" + out_rel_name  # with no extension

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if os.path.isfile(data_dir + "/PCU2.rsp"):
        subprocess.call(["cp", data_dir + "/PCU2.rsp", \
                out_dir + "/" + rsp_matrix])
    else:
        print("ERROR: Response matrix doesn't exist in the reduced data "\
                "directory.")
        exit()

    if "heasoft" not in DYLD_LIB_PATH:
        subprocess.Popen(HEAINIT, shell=True)

    try:
        hdulist = fits.open(ccf_data_file)
    except IOError:
        print "\tERROR: File does not exist: %s" % ccf_data_file
        exit()

    obs_time = hdulist[0].header['EXPOSURE']
    dt = hdulist[0].header['DT']
    n_seg = hdulist[0].header['SEGMENTS']
    detchans = hdulist[0].header['DETCHANS']
    n_bins = hdulist[0].header['N_BINS']
    hdulist.close()

    if os.path.isfile(parfit_file):
        fake_qpo_spectra.main(out_rootname, parfit_file, dt=dt, \
                rsp_matrix=rsp_matrix, exposure=obs_time, test=testing)
    else:
        print("\t ERROR: Parameter fit file does not exist: %s" % parfit_file)
        exit()

    varying_parameters = read_all_log_files(es_dir, out_rel_name, boot_num)
    avg_var_pars = avg_boot_var_pars(varying_parameters, boot_num=boot_num)

    avg_var_pars = avg_fit_func(avg_var_pars, parfit_file)

    plot_final_multifit(avg_var_pars, out_rootname=out_rootname, \
        parfit_file=parfit_file)

    plot_final_lagenergy(lag_data_file, out_rootname=out_rootname)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__, epilog="For optional"\
            " arguments, default values are given in brackets at end of "\
            "description.")

    parser.add_argument("--prefix", default="GX339-BQPO", dest="prefix",
            help="Identifying prefix of data (object nickname or data ID). "\
            "[GX339-BQPO]")

    parser.add_argument("--dt_mult", default=64, dest="dt_mult",
            type=tools.type_positive_int, help="Multiple of dt (dt is time "\
            "resolution from data file) for timestep between bins. Must be a "\
            "power of 2, positive, int. [64]")

    parser.add_argument("--nsec", default=64, dest="n_seconds",
            type=tools.type_positive_int, help="Number of seconds in each "\
            "Fourier segment. Must be a power of 2, positive, int. [64]")

    parser.add_argument("--testing", type=int, default=0, choices={0,1},
            dest='test', help="Int flag: 0 if computing all segments, 1 if "\
            "computing only one segment for testing. [0]")

    parser.add_argument("--day", default=datetime.now().strftime("%y%m%d"),
            dest="day", help="Date of the analysis, in YYMMDD format. "\
            "[%s]" % datetime.now().strftime("%y%m%d"))

    parser.add_argument("--boot", default=2, dest="boot_num", type=int,
            help="Number of bootstrap realizations done. [2]")

    args = parser.parse_args()

    main(prefix=args.prefix, dt_mult=args.dt_mult, n_seconds=args.n_seconds,
            day=args.day, boot_num=args.boot_num, testing=args.test)