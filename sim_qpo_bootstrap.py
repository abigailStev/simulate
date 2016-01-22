#!/usr/bin/env python
"""
Makes fake lag-energy spectra given the energy spectral parameters from
multifit_plots.py in sed_fit_bootstrap.sh. Calls fake_qpo_spectra.py.
Used in ccf_bootstrap.sh.

Be sure that the directories point correctly for your setup.

Notes: HEASOFT 6.14.* and Python 2.7.* (with supporting libraries) must be
installed in order to run this script.

Example call:  python sim_qpo_bootstrap.py --prefix cygx1

"""
__author__ = 'Abigail Stevens <A.L.Stevens at uva.nl>'
__year__ = "2015"

import argparse
import numpy as np
import os.path
import subprocess
from datetime import datetime
from astropy.table import Table
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

## These are things I've written.
## Their directories are in my PYTHONPATH bash environment variable.
import fake_qpo_spectra
import multifit_plots  ## in energy_spectra
import overplot_lags  ## in lag_spectra
from fake_qpo_spectra import chisquared_lagspectrum as lag_chisq  ## in lag_spectra
import sed_pars  ## in energy_spectra
import tools  ## in whizzy_scripts

## Global variables -- change these to suit your setup
HOME_DIR = os.path.expanduser("~")
CCF_BASEDIR = HOME_DIR + "/Dropbox/Research/cross_correlation"
ES_BASEDIR = HOME_DIR + "/Dropbox/Research/energy_spectra"
LAG_BASEDIR = HOME_DIR + "/Dropbox/Research/lag_spectra"
EXE_DIR = HOME_DIR + "/Dropbox/Research/simulate"
DATA_BASEDIR = HOME_DIR + "/Reduced_data"
DUMP_FILE = "dump.txt"  # Name of dumping file for intermediary steps

## Gets these from 'source'/'export' in sed_fit_bootstrap.sh; if nothing,
## uses the default value (the second argument).
N_SPECTRA = int(os.getenv('n_spectra', 24))
# N_PARAMS = int(os.getenv('n_params', 9))
# FIT_SPECIFIER = os.getenv('fit_specifier', "1BB-FS-G-Tin-fzs-fzNbb")
N_PARAMS = int(os.getenv('n_params', 11))
FIT_SPECIFIER = os.getenv('fit_specifier', "2BB-FS-G-kT-fzs-fzNbb8857")

print "Fit specifier: %s" % FIT_SPECIFIER
print "Number of spectra: %d" % N_SPECTRA
print "Number of parameters: %d" % N_PARAMS

## Gets HEASOFT directory name and script name
DYLD_LIB_PATH = os.environ.get('DYLD_LIBRARY_PATH')
HEADAS_DIR = os.environ.get('HEADAS')
HEAINIT = ". " + HEADAS_DIR + "/headas-init.sh"

## The file type you want the figures to be
## 'eps' is recommended for journal-quality (but not smartphone-friendly) plots
## 'png' is recommended for lower-quality (but universally open-able) plots
PLOT_EXT = "eps"


################################################################################
def plot_final_multifit(var_pars, out_rootname="out", title=" "):
    """
    Plot the SED parameter variations averaged over all the bootstrapped runs.
    The error on the parameters and the function come from the bootstrap-to-
    bootstrap variations.

    Parameters
    ----------
    var_pars : np.array of sed_pars.Parameter
        1-D array of the varying SED parameters, averaged

    out_rootname : str
        Full path root name (directories plus base) of the output file. Used for
        the data file and the plot file.

    Returns
    -------
    Nothing, but makes a plot.
    """
    print "Var par shape:", np.shape(var_pars)
    plot_file = out_rootname + "." + PLOT_EXT

    multifit_plots.make_var_plots(plot_file, N_SPECTRA, var_pars, quiet=False,
            title=title)
    subprocess.call(['cp', plot_file,
          "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])


################################################################################
def plot_final_lagenergy(lag_data_file, out_rootname="out", mod=" "):
    """
    Plot the lag-energy of this simulation (generated from the average of the
    bootstraps) with the data, and gives the chisquared/dof fit.

    Caveat: These degrees of freedom should be taken with a grain of salt, since
    technically the simulation is generated from the data. But this is the best
    estimate we can come up with for the DOF for now.

    Parameters
    ----------
    lag_data_file : str
        The lag-energy spectrum from the original data that we are fitting to.

    out_rootname : str
        The dir+basename of the output file that the simulated lag-energy
        spectrum was written to, and to write the lag-energy plot to.

    Returns
    -------
    Nothing, but makes a plot.

    Files created
    -------------
    *_lag.eps :
        The lag-energy plot of the data with this simulation (generated from the
        average of the bootstraps).

    """

    in_file_list = [lag_data_file, out_rootname + "_lag.fits"]
    labels = ["Data", r"%s" % mod]
    overplot_lags.make_plot(in_file_list, labels, out_rootname)
    chisquared, dof = lag_chisq(lag_data_file, out_rootname + "_lag.fits")
    print "\n\tLag-energy fit: %.5f / %d" % (chisquared, dof)


################################################################################
def plot_deltaphase_hist(phase_points, out_rootname="out", title=" "):
    """
    Plot histograms of the delta phase for param 0 and param 2 wrt param 1.

    Parameters
    ----------
    phase_points : np.array of floats
        2-D array of the phase of each parameter from each bootstrap.
        Size = (number of varying parameters, number of bootstraps).

    out_rootname : str
        The dirname+basename of the output file, to be appended to for saving
        the plot file.

    Returns
    -------
    Nothing, but writes to file *_dphase_hist.PLOT_EXT

    """
    xLocator = MultipleLocator(0.02)  ## loc of minor ticks on x-axis
    font_prop = font_manager.FontProperties(size=20)

    plot_file = out_rootname + "_dphase_hist." + PLOT_EXT

    delta_phase_1 = np.abs(phase_points[1,:] - phase_points[0,:])
    print np.mean(delta_phase_1)
    delta_phase_2 = np.abs(phase_points[1,:] - phase_points[2,:])
    print np.mean(delta_phase_2)

    # hist_1, edges_1 = np.histogram(delta_phase_1, bins=50, range=(0, 1))
    # hist_2, edges_2 = np.histogram(delta_phase_2, bins=50, range=(0, 1))

    fig, ax = plt.subplots(1, 1, figsize=(10, 7.5), dpi=300, tight_layout=True)
    bins_1, edges_1, patches_1 = ax.hist(delta_phase_1, bins=100, range=[0, 1],
            facecolor='red')
    bins_2, edges_2, patches_2 = ax.hist(delta_phase_2, bins=100, range=[0, 1],
            facecolor='blue')
    ax.set_xlabel(r'$\Delta$ Normalized parameter phase',
            fontproperties=font_prop)
    ax.set_ylabel('Bootstrap iterations', fontproperties=font_prop)
    ax.set_xlim(0, 0.41)
    ax.set_ylim(0, 3500)
    ax.set_xticks(np.arange(0, 0.41, 0.1))
    ax.xaxis.set_minor_locator(xLocator)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.set_title(r'%s' % title, fontproperties=font_prop)

    plt.savefig(plot_file)
    plt.close()
    print("Plot saved to: %s" % plot_file)
    subprocess.call(['open', plot_file])
    subprocess.call(['cp', plot_file,
            "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])


################################################################################
def plot_phase_contours(phase_points, avg_var_pars, out_rootname="out", title=" "):
    """
    Make a contour plot (currently just a scatter plot) of the phases of two
    parameters from each bootstrap.

    Parameters
    ----------
    phase_points : np.array of floats
        2-D array of the phase of each parameter from each bootstrap.
        Size = (number of varying parameters, number of bootstraps).

    avg_var_pars : np.array or list of str
        1-D array of the names of the varying parameters.

    out_rootname : str
        The dirname+basename of the output file, to be appended to for saving
        the plot file.

    """
    xLocator = MultipleLocator(0.02)  ## loc of minor ticks on x-axis
    font_prop = font_manager.FontProperties(size=20)

    plot_file = out_rootname + "_phasecontour1." + PLOT_EXT
    fig, ax = plt.subplots(1, 1, figsize=(10, 7.5), dpi=300, tight_layout=True)
    x = phase_points[0,:] + 1
    y = phase_points[1,:] + 1
    ax.scatter(x, y, color='black')

    # cov_xy = np.cov(x, y, ddof=1)
    # # print cov_xy
    # # print np.shape(cov_xy)
    # # sigma_xy = cov_xy[0, 1] * 100
    # # sigma_x = cov_xy[0, 0] * 100
    # # sigma_y = cov_xy[1, 1] * 100
    # mu_x = np.mean(x)
    # sigma_x = np.std(x, ddof=1)
    # mu_y = np.mean(y)
    # sigma_y = np.std(y, ddof=1)
    # print sigma_x, sigma_y

    # g = models.Gaussian2D(amplitude=1.0, x_mean=mu_x, y_mean=mu_y,
    #         cov_matrix=cov_xy)
    # g = models.Gaussian2D(amplitude=1.0, x_mean=mu_x, y_mean=mu_y,
    #         x_stddev=sigma_x, y_stddev=sigma_y, theta=0.785398)
    # print "Is this working?"
    # fit_g = fitting.LevMarLSQFitter()
    # print "How about now?"
    # with warnings.catch_warnings():
    #     warnings.simplefilter('ignore')
    #     g = fit_g(g_init, np.arange(0.75, 0.85, 0.01),
    #             np.arange(0.75, 0.85, 0.01), data)

    # plt.imshow(g(x,y), origin='lower', interpolation='nearest')
    # plt.imshow(data, origin='lower', interpolation='nearest')

    # X,Y = np.meshgrid(np.arange(0, 1, 0.01), np.arange(0, 1, 0.01))
    # X,Y = np.meshgrid(np.arange(0.75, 0.85, 0.01), np.arange(0.75, 0.85, 0.01))
    # Z = mlab.bivariate_normal(X, Y, sigmax=sigma_x, sigmay=sigma_y, mux=mu_x,
    #         muy=mu_y)
    # CS = plt.contour(X,Y,Z)
    # ax.clabel(CS)
    ax.plot(np.mean(x), np.mean(y), '*',
            mfc='white', mew=1, mec='yellow', ms=10)
    ax.set_xlabel(r"Normalized phase: %s" % avg_var_pars[0],
            fontproperties=font_prop)
    ax.set_ylabel(r"Normalized phase: %s" % avg_var_pars[1],
            fontproperties=font_prop)
    # ax.set_xticks(np.arange(0, 1.05, 0.1))
    ax.xaxis.set_minor_locator(xLocator)
    # ax.set_yticks(np.arange(0, 1.05, 0.1))
    ax.yaxis.set_minor_locator(xLocator)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.set_title(r'%s' % title, fontproperties=font_prop)
    plt.savefig(plot_file)
    plt.close()
    print("Plot saved to: %s" % plot_file)
    subprocess.call(['open', plot_file])
    subprocess.call(['cp', plot_file,
          "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])

    plot_file = out_rootname+"_phasecontour2." + PLOT_EXT
    x = phase_points[0,:] + 1
    y = phase_points[2,:]
    fig, ax = plt.subplots(1, 1, figsize=(10, 7.5), dpi=300, tight_layout=True)
    ax.scatter(x, y, color='black')
    ax.plot(np.mean(x), np.mean(y), '*', mfc='white', mew=1, mec='yellow',
            ms=10)
    ax.set_xlabel(r"Normalized phase: %s" % avg_var_pars[0],
            fontproperties=font_prop)
    ax.set_ylabel(r"Normalized phase: %s" % avg_var_pars[2],
            fontproperties=font_prop)
    # ax.set_xticks(np.arange(0, 1.05, 0.1))
    ax.xaxis.set_minor_locator(xLocator)
    # ax.set_yticks(np.arange(0, 1.05, 0.1))
    ax.yaxis.set_minor_locator(xLocator)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.set_title(r'%s' % title, fontproperties=font_prop)
    plt.savefig(plot_file)
    plt.close()
    print("Plot saved to: %s" % plot_file)
    subprocess.call(['open', plot_file])
    subprocess.call(['cp', plot_file,
          "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])

    plot_file = out_rootname+"_phasecontour3." + PLOT_EXT
    x = phase_points[1,:] + 1
    y = phase_points[2,:]
    fig, ax = plt.subplots(1, 1, figsize=(10, 7.5), dpi=300, tight_layout=True)
    ax.scatter(x, y, color='black')
    ax.plot(np.mean(x), np.mean(y), '*', mfc='white', mew=1, mec='yellow',
            ms=10)
    ax.set_xlabel(r"Normalized phase: %s" % avg_var_pars[1],
            fontproperties=font_prop)
    ax.set_ylabel(r"Normalized phase: %s" % avg_var_pars[2],
            fontproperties=font_prop)
    # ax.set_xticks(np.arange(0, 1.05, 0.1))
    ax.xaxis.set_minor_locator(xLocator)
    # ax.set_yticks(np.arange(0, 1.05, 0.1))
    ax.yaxis.set_minor_locator(xLocator)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.set_title(r'%s' % title, fontproperties=font_prop)
    plt.savefig(plot_file)
    plt.close()
    print("Plot saved to: %s" % plot_file)
    subprocess.call(['open', plot_file])
    subprocess.call(['cp', plot_file,
            "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])


################################################################################
def avg_boot_var_pars(var_pars, boot_num=2):
    """
    Computes the average value of the varying parameters across the bootstrap
    iterations. The value at each QPO phase is the mean over all bootstraps,
    and the error is the square root of the variance of the values over all
    bootstraps.

    Parameters
    ----------
    var_pars : np.array of sed_pars.Parameter()
        2-D array of the SED parameters that vary with QPO phase, over
        all the bootstrap iterations.

    boot_num : int
        Number of bootstrap iterations implemented.

    Returns
    -------
    avg_var_pars : np.array of sed_pars.Parameter()
        1-D array; similar to input var_pars but averaged over all bootstrap
        iterations.

    bad_boot : int
        The number of bootstrap iterations that gave a 'nan' fit. Should not be
        more than ~1/1000 of boot_num (for large boot_nums).

    """
    n_varpar = np.shape(var_pars)[1]
    n_boot = np.shape(var_pars)[0]  ## Number of bootstraps we keep, since a few
                                    ## might have 'nan' covariance matrices from
                                    ## a screwy fit function.
    bad_boot = boot_num - n_boot
    print "Number of varying parameters:", n_varpar
    print "Number of good bootstraps:", n_boot
    print "Number of bad bootstraps:", boot_num - n_boot

    avg_var_pars = np.asarray([])
    for par in range(n_varpar):
        parameter_values = np.zeros(N_SPECTRA)

        for boot in range(n_boot):
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

    return avg_var_pars, bad_boot


################################################################################
def read_xspec_log_files(es_dir, out_rel_name, boot_num=2):
    """
    Read in all XSPEC log files (with chatter set to 4) that were generated in
    sed_fit_bootstrap.sh, and append each bootstrap iteration's values to its
    sed_pars.Parameter.

    Parameters
    ----------
    es_dir : str
        The directory with all the energy spectroscopy files from
        sed_fit_bootstrap.sh.

    out_rel_name : str
        The relative (i.e. local) name for the output files.

    boot_num : int
        Number of bootstrap iterations implemented.

    Returns
    -------
    var_pars : np.array of sed_pars.Parameter()
        2-D array of the SED parameters that vary with QPO phase, over
        all the bootstrap iterations.
    """
    all_parameters = np.tile(sed_pars.Parameter(), N_PARAMS)
    untied = []
    for i in range(1, boot_num + 1):
    # for i in range(1, 10):
        log_file = es_dir + "/" + out_rel_name + "_b-" + str(i) + "_xspec.log"
        # print log_file
        if os.path.isfile(log_file):
            boot_sed_parameters, n_spec_log = multifit_plots.read_log_file(log_file,
                    quiet=True)
            assert n_spec_log == N_SPECTRA
            boot_var_par = []
            for component in boot_sed_parameters:
                # print component.mod_name, component.par_name
                if len(component.value) > 1:
                    if component.value[0] != component.value[1] or \
                            component.value[7] != component.value[0]:
                        component.varying = True
                        boot_var_par.append(component)

            all_parameters = np.vstack((all_parameters, boot_sed_parameters))
            untied.append(boot_var_par)
        else:
            pass

    all_parameters = all_parameters[1:]
    print np.shape(all_parameters)
    untied = np.asarray(untied)
    # print untied
    print np.shape(untied)
    n_untied = np.shape(untied)[1]
    # print var_pars.__str__()
    return untied, n_untied


################################################################################
def avg_best_fit(avg_untied, parfit_file, n_untied=3, boot_num=2):
    """
    Gets the average of the best-fit function parameters from each bootstrap
    iteration.

    Parameters
    ----------
    avg_untied : np.array of sed_pars.Parameter objects
        1-D array of the untied spectral parameters, averaged over bootstrap
        iteration.

    parfit_file : str
        Full path of the "_funcfit.txt" file that was made in
        energy_spectra/multifit_plots.py and read in by ./fake_qpo_spectra.py.

    n_untied : int
        The number of untied spectral parameters in the model. [3]

    boot_num : int
        Number of bootstrap iterations implemented. [2]

    Returns
    -------
    avg_untied : np.array of sed_pars.Parameter objects
        Same as input, but with best_fit, func_fit, phase, and phase_err
        assigned.

    phase_points : np.array of floats
        2-D array of the phases from each bootstrap,
        Size = (good bootstraps, number of untied parameters)

    model_name : str
        The name of the spectral model, used for printing on top of plots.

    """
    all_bestfits = np.zeros((n_untied, 5, 1))

    with open(parfit_file, 'r') as f:
        j = 0
        for line in f:
            i = 0
            boot_bestfit = np.zeros(5)
            line_element = line.strip().split("    ")
            model_name = line_element[0][1:-1]
            for element in line_element[1:]:
                if "[" in element or "]" in element or "," in element:
                    parameter_bestfit = np.array(element.replace('[',
                            '').replace(']', '').split(','), dtype=np.float)
                    boot_bestfit = np.vstack((boot_bestfit, parameter_bestfit))
                    i += 1
            if np.isnan(boot_bestfit).any():
                print "It's nan:, %d" % j
            else:
                all_bestfits = np.dstack((all_bestfits, boot_bestfit[1:,]))
            j += 1
    all_bestfits = all_bestfits[:,:,1:]
    print np.shape(all_bestfits)

    phase_points = np.zeros(np.shape(all_bestfits)[-1])
    for par,i in zip(avg_untied, range(n_untied)):
        par.best_fit = np.mean(all_bestfits[i,:,:], axis=-1)
        par.funcfit = fake_qpo_spectra.fit_function(np.arange(-0.02, 1.02,
                0.01), par.best_fit)
        phase = all_bestfits[i,1,:] / (2.0 * np.pi)
        phase_points = np.vstack((phase_points, phase))
        par.phase = np.mean(phase)
        par.phase_err = np.sqrt(np.var(phase, ddof=1))
        print par.mod_name, par.par_name
        print "\tAverage phase:", par.phase
        print "\tErr on phase:", par.phase_err

    phase_points = phase_points[1:,]

    return avg_untied, phase_points, model_name


################################################################################
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
    Nothing
    """

    ## Initializations of file names and directories
    ccf_dir = CCF_BASEDIR + "/out_ccf/" + prefix
    es_dir = ES_BASEDIR + "/out_es/" + prefix + "/bootstrapped"
    # es_dir = ES_BASEDIR + "/out_es/" + prefix
    lag_dir = LAG_BASEDIR + "/out_lags/" + prefix
    out_dir = EXE_DIR + "/out_sim/" + prefix + "/bootstrapped"
    # out_dir = EXE_DIR + "/out_sim/" + prefix
    data_dir = DATA_BASEDIR + "/" + prefix
    rsp_matrix = prefix + "_PCU2.rsp"
    parfit_file = es_dir + "/" + prefix + "_" + day + "_" + FIT_SPECIFIER + \
            "_funcfit.txt"
    # parfit_file = es_dir + "/" + prefix + "_" + day + "_" + FIT_SPECIFIER + \
    #         "_funcfit_test.txt"
    ccf_data_file = ccf_dir + "/" + prefix + "_" + day + "_t" + str(dt_mult) + \
            "_" + str(n_seconds) + "sec_adj.fits"
    lag_data_file = lag_dir + "/" + prefix + "_" + day + "_t" + str(dt_mult) + \
            "_" +str(n_seconds) + "sec_adj_lag.fits"
    print lag_data_file

    epoch = 5
    out_rel_name = prefix + "_" + day + "_" + FIT_SPECIFIER
    out_rootname = out_dir + "/" + out_rel_name  # with no extension

    ## Make sure the output directory exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ## Copy the response matrix into the reduced data directory, if it doesn't
    ## exist
    if os.path.isfile(data_dir + "/PCU2.rsp"):
        subprocess.call(["cp", data_dir + "/PCU2.rsp", \
                out_dir + "/" + rsp_matrix])
    else:
        print("ERROR: Response matrix doesn't exist in the reduced data "\
                "directory.")
        exit()

    ## Make sure heasoft is running
    if "heasoft" not in DYLD_LIB_PATH:
        subprocess.Popen(HEAINIT, shell=True)

    ## Get parameters from the ccf data file header
    try:
        ccf_table = Table.read(ccf_data_file)
    except IOError:
        print("\tERROR: File does not exist: %s" % ccf_data_file)
        exit()

    obs_time = ccf_table.meta['EXPOSURE']
    n_seg = ccf_table.meta['SEGMENTS']
    detchans = ccf_table.meta['DETCHANS']
    n_bins = ccf_table.meta['N_BINS']
    n_seconds = ccf_table.meta['SEC_SEG']
    dt = ccf_table.meta['DT']

    # ## Run fake_qpo_spectra.py
    # if os.path.isfile(parfit_file):
    #     fake_qpo_spectra.main(out_rootname, parfit_file, dt=dt, \
    #             rsp_matrix=rsp_matrix, exposure=obs_time, test=testing,
    #             n_params=N_PARAMS)
    # else:
    #     print("\t ERROR: Parameter fit file does not exist: %s" % parfit_file)
    #     exit()
    # print "Finished fake_qpo_spectra"

    ## Read in values of all parameters from the XSPEC log files
    untied_parameters, n_untied = read_xspec_log_files(es_dir, out_rel_name,
            boot_num=boot_num)
    avg_untied, bad_boot = avg_boot_var_pars(untied_parameters,
            boot_num=boot_num)


    ## Compute the average fit function to the untied parameter variations
    ## Using boot_num-bad_boot as boot_num, because we're not counting the bad
    ## bootstrap runs.
    avg_untied, phase_points, model_name = avg_best_fit(avg_untied, parfit_file,
            n_untied, boot_num=boot_num-bad_boot)

    # print model_name
    model_name = model_name.replace("phabs*", "phabs$\\times\\,$")
    # print model_name
    np.savetxt("%s_phasepoints.txt" % FIT_SPECIFIER, phase_points)

    ## Plot the "multi-fit" parameter variations to the average values
    plot_final_multifit(avg_untied, out_rootname=out_rootname, title=model_name)

    ## Make contour plots of the phases
    plot_phase_contours(phase_points, avg_untied, out_rootname=out_rootname,
            title=model_name)
    plot_deltaphase_hist(phase_points, out_rootname=out_rootname,
            title=model_name)

    ## Plot the simulated lag-energy vs the data
    plot_final_lagenergy(lag_data_file, out_rootname=out_rootname,
            mod=model_name)
    print "Done!"


################################################################################
if __name__ == "__main__":

    #########################################
    ## Parse input arguments and call 'main'
    #########################################

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
