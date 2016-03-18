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
__year__ = "2015-2016"

import argparse
import numpy as np
import os.path
import subprocess
from scipy.optimize import leastsq
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
# N_PARAMS = int(os.getenv('n_params', 10))
# FIT_SPECIFIER = os.getenv('fit_specifier', "pBB-FS-G-p-fzs-fzNbb")
N_PARAMS = int(os.getenv('n_params', 11))
FIT_SPECIFIER = os.getenv('fit_specifier', "2BB-FS-G-kT-fzs-fzNbb8857-2")

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
def check_equal(iterator):
    return len(set(iterator)) <= 1


################################################################################
def plot_final_multifit(var_pars, out_rootname="out", title=" "):
    """
    Plot the spectral parameter variations and the fit function.
    The error on the parameters and the function are the rms of the
    bootstrap-to-bootstrap variations.

    Parameters
    ----------
    var_pars : np.array of sed_pars.Parameter
        1-D array of the varying SED parameters, averaged

    phi_max_err : np.array of floats
        1-D array of the error on the phase of the maximum value of the fit
        function, taken as the standard error on the bootstrapping.

    out_rootname : str
        Full path root name (directories plus base) of the output file. Used for
        the data file and the plot file.

    title : str
        The title to write on the plot, typically the XSPEC model name.

    Returns
    -------
    Nothing, but makes a plot and saves it to out_rootname.PLOT_EXT
    """
    print "Var par shape:", np.shape(var_pars)
    print var_pars
    plot_file = out_rootname + "." + PLOT_EXT

    multifit_plots.make_var_plots(plot_file, N_SPECTRA, var_pars, quiet=False,
            title=title)
    subprocess.call(['cp', plot_file,
          "/Users/abigailstevens/Dropbox/Research/CCF_paper1/images/"])


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

    delta_phase_1 = phase_points[0,:] - phase_points[1,:]
    print np.mean(delta_phase_1)
    delta_phase_2 = phase_points[2,:] - phase_points[1,:]
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
    ax.set_ylim(0, 5000)
    ax.set_xticks(np.arange(0, 0.41, 0.1))
    ax.xaxis.set_minor_locator(xLocator)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.set_title(r'%s' % title, fontproperties=font_prop)

    plt.savefig(plot_file)
    plt.close()
    print("Plot saved to: %s" % plot_file)
    subprocess.call(['open', plot_file])


################################################################################
def plot_phase_contours(phase_points, avg_var_pars, out_rootname="out",
        title=" "):
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


################################################################################
def write_final_besfit_pars(write_func, mod_string, mod_parameters):
    """
    Writes the final fit function best fit parameters to the file write_func.
    Only writes for the varying parameters.

    Parameters
    ----------
    write_func : str
        Filename where the best fit parameters will be written. Overwrites any
        previous file.

    mod_string : string
        The full xspec model name.

    mod_parameters : list of sed_pars.Parameter objects
        1-D list of the model parameters.

    Returns
    -------
    Nothing, but overwrites the file write_func.

    """
    if write_func != "":
        # if not quiet:
        #     print("Writing function parameters to: %s" % write_func)
        with open(write_func, 'w') as out:
            out.write("%s    " % (mod_string))
            for parameter in mod_parameters:
                if parameter.varying:
                    # print parameter.best_fit
                    try:
                        out.write("[%.4e,%.4e,%.4e,%.4e,%.4e]    " % \
                                (parameter.best_fit[0], parameter.best_fit[1],
                                 parameter.best_fit[2], parameter.best_fit[3],
                                 parameter.best_fit[4]))
                    except TypeError:
                        out.write("[nan,nan,nan,nan,nan]    ")
                else:
                    out.write("%.4e    " % parameter.value[0])
            out.write("\n")


################################################################################
def fit_function(t, p):
    """
    Computing a function to fit to the SED parameter variations.

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
def function_residuals(p, data, data_err, t):
    """
    Getting the residual of the data with the current fit function. Dividing by
    error bar to weight it appropriately like in weighted least squares, e.g.
    S. Vaughan 2013 eqn 6.12 (modified -- not squaring because according to
    scipy.optimize.leastsq documentation, it will square for me to compute the
    real residual and covariance matrix (which will also make it look exactly
    like eqn 6.12))

    Parameters
    ------
    p : np.array of floats
        1-D array of the function parameters.

    data : np.array of floats
        1-D array of the data we want to fit to; in this case, the list of SED
        fit parameters over QPO phase.

    data_err : np.array of floats
        1-D array of the error on the data.

    t : np.array of floats
        1-D array of the time steps for the fitting function.

    Returns
    -------
    np.array of floats
        1-D array of a modified weighted least squared residual of the current
        function fit with the data. From S. Vaughan 2013 eqn 6.12 and
        scipy.optimize.leastsq documentation.
    """
    residual = np.abs(data - fit_function(t, p)) / data_err
    return residual


################################################################################
def get_phase(parameter, num_spectra, quiet):
    """
    Fitting a function to an energy spectrum fit parameter to determine the
    phase of the parameter changes.

    Parameters
    ----------
    parameter : Parameter object
        The spectral energy distribution parameter.

    num_spectra : int
        The number of energy spectra in use (the number of energy spectra per
        QPO phase).

    quiet : bool
        If True, suppresses printing to the screen.

    Returns
    -------
    Parameter object
        The energy spectrum parameter, with funcfit, phase, and phase_err
        assigned.

    """
    t = np.arange(num_spectra) / 23.5
    p = [1., 0., 1., 0., np.mean((np.min(parameter.value), \
            np.max(parameter.value)))]  ## Amplitude, phase shift, mean

    parameter.error = np.mean((parameter.pos_err, parameter.neg_err), axis=0)

    p_best = leastsq(function_residuals, p, args=(parameter.value, \
            parameter.error, t), full_output=1)
    # print "P best:", p_best
    best_fit = p_best[0]

    # if not quiet:
    #     print("\tBest fit: %s" % str(best_fit))

    # plt.errorbar(t, parameter.value, xerr=None, yerr=parameter.error)
    # plt.plot(t, fit_function(t, best_fit))
    # plt.xlim(0,1)
    # plt.show()

    ## Error on phase from S. Vaughan 2013 p 168
    bonus_matrix = p_best[1]  ## A Jacobian approximation to the Hessian of the
            ## least squares objective function.
    resid_var = np.var(function_residuals(best_fit, parameter.value, \
            parameter.error, t), ddof=1)
    ## As outlined in the scipy.optimize.leastsq documentation, multiply the
    ## bonus matrix by the variance of the residuals to get the covariance
    ## matrix.

    try:
        cov_matrix = bonus_matrix * resid_var
    except TypeError:
        # print("\t %s" % str(resid_var))
        # print("\t %s" % str(bonus_matrix))
        parameter.best_fit = np.nan
        parameter.funcfit = np.nan
        parameter.phase = np.nan
        parameter.phase_err = np.nan
        return parameter

    parameter.best_fit = best_fit
    if num_spectra == 24:
        parameter.funcfit = fit_function(np.arange(-0.02, 1.02, 0.001), best_fit)
    elif num_spectra == 47:
        parameter.funcfit = fit_function(np.arange(-0.02, 2.02, 0.001), best_fit)
    else:
        parameter.funcfit = fit_function(np.arange(-1.02, 2.02, 0.001), best_fit)
    parameter.phase = best_fit[1] / (2.0 * np.pi)
    parameter.phase_err = np.sqrt(cov_matrix[1][1]) / (2.0 * np.pi)

    return parameter


################################################################################
def get_boot_varpar_err(var_par_vals, xspec_log_data_file):
    """
    Computes the average value of the varying parameters across the bootstrap
    iterations. The value at each QPO phase is the mean over all bootstraps,
    and the error is the square root of the variance of the values over all
    bootstraps.

    Parameters
    ----------
    var_par_vals : np.array of floats
        2-D array of the values of the spectral parameters that vary with QPO
        phase, over all the bootstrap iterations.

    xspec_log_data_file : str
        File name of the XSPEC log of the fit to the real data, with chatter
        set to 4.

    boot_num : int
        Number of bootstrap iterations implemented.

    Returns
    -------
    data_var_par : np.array of sed_pars.Parameter()
        1-D array; the varying parameters, with values from the data XSPEC log
        and errors from the rms of the bootstrap parameter values.

    """
    print "Shape of varpar vals:", np.shape(var_par_vals)
    n_varpar = np.shape(var_par_vals)[1]
    # print "Number of varying parameters:", n_varpar

    data_var_par = read_data_xspec_log(xspec_log_data_file)
    assert len(data_var_par) == n_varpar

    if np.shape(var_par_vals)[-1] == 1:
        err = np.ones((N_SPECTRA, N_PARAMS))
    else:
        err = np.sqrt(np.var(var_par_vals, axis=-1, ddof=1))
    print err

    for i in range(n_varpar):
        assert len(err[:,i]) == N_SPECTRA, "Untied parameter value errors "\
                "should have length N_SPECTRA"
        data_var_par[i].pos_err = err[:,i]
        data_var_par[i].neg_err = err[:,i]

    return data_var_par


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
    all_par_vals = np.zeros((N_SPECTRA, N_PARAMS))

    for i in range(1, boot_num + 1):
    # for i in range(1, 10):
        boot_log_file = es_dir + "/" + out_rel_name + "_b-" + str(i) + \
                        "_xspec.log"
        # print log_file
        if os.path.isfile(boot_log_file):
            par_vals = xspec_boot_log_to_array(boot_log_file)
            # print "Shape of par vals:", np.shape(par_vals)
            all_par_vals = np.dstack((all_par_vals, par_vals))
        else:
            pass

    # print "All par vals:", np.shape(all_par_vals)
    all_par_vals = all_par_vals[:,:,1:]
    # print "All par vals:", np.shape(all_par_vals)

    good_boots = np.shape(all_par_vals)[-1]
    # print "Good boots:", good_boots
    n_varpars = 0
    delete_index = []

    for j in range(N_PARAMS):
        if not check_equal(all_par_vals[:,j,0].flatten()):
            n_varpars += 1
        else:
            delete_index.append(j)
    # print "Shape all par vals:", np.shape(all_par_vals)
    untied_varpar_vals = np.delete(all_par_vals, delete_index, axis=1)
    # print "Shape untied par vals:", np.shape(untied_varpar_vals)
    # print untied_varpar_vals

    return untied_varpar_vals, n_varpars, good_boots


################################################################################
def xspec_boot_log_to_array(log_file):
    if not os.path.isfile(log_file) or os.path.getsize(log_file) == 0:
        print log_file
        raise Exception("Log file does not exist or is empty.")

    par_names = np.array(['nH', 'Tin', ' p ', 'Gamma', 'FracSctr', 'UpScOnly',
            'norm', 'LineE', 'Sigma', 'Redshift', 'kT', 'PhoIndex'])
    comp_names = np.array(['phabs', 'wabs', 'diskpbb', 'diskbb', 'bbodyrad', \
            'simpl', 'simpler', 'nthcomp', 'gaussian', 'diskline', 'cutoffpl'])

    par_values_1D = []

    with open(log_file, 'r') as f:
        for line in f:
            if any(elt in line for elt in par_names) and \
                    any(elt in line for elt in comp_names):
                if "frozen" in line:
                    par_values_1D = np.append(par_values_1D,
                            float(line.split()[-2]))
                else:
                    par_values_1D = np.append(par_values_1D,
                                float(line.split()[-3]))

    try:
        par_values = np.reshape(par_values_1D, (N_SPECTRA, N_PARAMS))
    except ValueError:
        print "N_SPECTRA and/or N_PARAMS are incorrect for this log file."
        par_values = np.asarray([])

    return par_values


################################################################################
def read_data_xspec_log(log_file):
    """
    Reads the XSPEC log file of the fit to the real data and assigns parameters
    to Parameter objects, with the expectation that chatter is set to 4.

    Parameters
    ----------
    log_file : str
        Full path of the XSPEC log file, with chatter set to 4. Assuming the
        spectral models listed below are the only ones used (or the only
        interesting ones).

    Returns
    -------
    varying_parameters : list of Parameter objects
        A 1-D list of untied parameters for the model.
    """
    if not os.path.isfile(log_file) or os.path.getsize(log_file) == 0:
        print log_file
        raise Exception("Log file does not exist or is empty.")

    mod_parameters = [sed_pars.Phabs().nH, sed_pars.Simpler().Gamma,
            sed_pars.Simpler().FracSctr, sed_pars.Simpler().UpScOnly,
            sed_pars.Simpl().Gamma, sed_pars.Simpl().FracSctr,
            sed_pars.Simpl().UpScOnly, sed_pars.Diskbb().Tin,
            sed_pars.Diskbb().norm, sed_pars.Diskpbb().Tin,
            sed_pars.Diskpbb().p, sed_pars.Diskpbb().norm,
            sed_pars.Bbodyrad().kT, sed_pars.Bbodyrad().norm,
            sed_pars.Gaussian().LineE, sed_pars.Gaussian().Sigma,
            sed_pars.Gaussian().norm]

    #################################################
    ## Reading in parameter values from the log file
    #################################################

    with open(log_file, 'r') as f:
        for line in f:
            for parameter in mod_parameters:
                if parameter.mod_name in line and parameter.par_name in line:
                    if "frozen" in line:
                        parameter.value = np.append(parameter.value,
                            float(line.split()[-2]))
                        parameter.par_num = np.append(parameter.par_num,
                            int(line.split()[1]))
                    else:
                        parameter.value = np.append(parameter.value,
                                float(line.split()[-3]))
                        parameter.par_num = np.append(parameter.par_num,
                                int(line.split()[1]))

    #############################################################
    ## Delete parameters if they're not used/present
    ## Determine if parameter varies across QPO phase or is tied
    ## Assign zero error to parameters
    #############################################################

    # num_spectra = np.amax([len(nth.Gamma.value), len(simpler.Gamma.value)])
    unused_parameters = []

    for parameter in mod_parameters:
        # print parameter.mod_name, parameter.par_name
        if len(parameter.value) > 1:
            # print parameter.value[0], parameter.value[1], parameter.value[3]
            if parameter.value[0] != parameter.value[1] or \
                    parameter.value[7] != parameter.value[0]:
                parameter.varying = True
                num_spectra = len(parameter.value)
        elif len(parameter.value) == 0:
            unused_parameters.append(parameter)

    for elt in unused_parameters:
        mod_parameters.remove(elt)

    var_pars = 0
    varying_parameters = []
    for parameter in mod_parameters:
        if parameter.varying:
            var_pars += 1
            varying_parameters.append(parameter)
            # print parameter.mod_name, parameter.par_name

    assert len(mod_parameters) == N_PARAMS

    return varying_parameters


################################################################################
def get_best_fit(avg_untied, boot_parfit_file, n_untied=3):
    """
    Gets the average of the best-fit function parameters from each bootstrap
    iteration.

    Parameters
    ----------
    avg_untied : np.array of sed_pars.Parameter objects
        1-D array of the untied spectral parameters, averaged over bootstrap
        iteration.

    boot_parfit_file : str
        Full path of the "_funcfit.txt" file that was made in
        energy_spectra/multifit_plots.py that contains the fit function
        parameters for all bootstraps, one bootstrap per line.

    data_parfit_file : str
        Full path of the "_funcfit.txt" file of the fit to the data that was
        made in energy_spectra/multifit_plots.py and read in by
        ./fake_qpo_spectra.py.

    n_untied : int
        The number of untied spectral parameters in the model. [3]

    Returns
    -------
    avg_untied : np.array of sed_pars.Parameter objects
        Same as input, but with best_fit, func_fit, phase, and phase_err
        assigned.

    [[Currently commented out]] phase_points : np.array of floats
        2-D array of the phases from each bootstrap,
        Size = (good bootstraps, number of untied parameters)

    model_name : str
        The name of the spectral model, used for printing on top of plots.

    """

    all_boot_bestfits = np.zeros((n_untied, 5, 1))
    boot_model_name = " "
    t_bins = np.arange(0,1,0.001)
    all_boot_bestfits = np.zeros((n_untied, 5, 1))
    all_max_phase_boot = np.zeros((n_untied,1))
    all_min_phase_boot = np.zeros((n_untied,1))

    with open(boot_parfit_file, 'r') as f:
        j = 0
        for line in f:
            i = 0
            boot_bestfit = np.zeros(5)
            line_element = line.strip().split("    ")
            boot_model_name = line_element[0][1:-1]
            for element in line_element[1:]:
                if "[" in element or "]" in element or "," in element:
                    parameter_bestfit = np.array(element.replace('[',
                            '').replace(']', '').split(','), dtype=np.float)
                    boot_bestfit = np.vstack((boot_bestfit, parameter_bestfit))
                    i += 1
            if np.isnan(boot_bestfit).any():
                print "It's nan:, %d" % j
            else:
                all_boot_bestfits = np.dstack((all_boot_bestfits,
                        boot_bestfit[1:,]))

                max_phase_boot = []
                min_phase_boot = []
                for k in range(0,3):
                    this_boot = np.array([boot_bestfit[k+1,0],
                                     boot_bestfit[k+1,1],
                                     boot_bestfit[k+1,2],
                                     boot_bestfit[k+1,3],
                                     boot_bestfit[k+1,4]])
                    func = fit_function(t_bins, this_boot)
                    max_phase_boot.append(t_bins[np.argmax(func)])
                    min_phase_boot.append(t_bins[np.argmin(func)])

                all_max_phase_boot = np.hstack((all_max_phase_boot,
                        np.reshape(np.array(max_phase_boot), (3,1))))
                all_min_phase_boot = np.hstack((all_min_phase_boot,
                        np.reshape(np.array(min_phase_boot), (3,1))))

            j += 1
    all_boot_bestfits = all_boot_bestfits[:,:,1:]
    all_max_phase_boot = all_max_phase_boot[:,1:]
    all_min_phase_boot = all_min_phase_boot[:,1:]
    print "Shape all boot bestfits:", np.shape(all_boot_bestfits)

    # phase_points = np.zeros(np.shape(all_boot_bestfits)[-1])
    for par,i in zip(avg_untied, range(n_untied)):

        par = get_phase(par, N_SPECTRA, False)
        # par.best_fit = np.mean(all_bestfits[i,:,:], axis=-1)
        # par.best_fit = data_bestfit[i]
        # par.funcfit = fake_qpo_spectra.fit_function(np.arange(-0.02, 1.02,
        #         0.01), par.best_fit)
        boot_phase = all_boot_bestfits[i,1,:] / (2.0 * np.pi)
        # phase_points = np.vstack((phase_points, boot_phase))
        par.phase = par.best_fit[1] / (2.0 * np.pi)
        par.phase_err = np.sqrt(np.var(boot_phase, ddof=1))
        print "\n\t", par.mod_name, par.par_name
        print "Phase:", par.phase
        print "Err on phase:", par.phase_err
        print "Bestfit:", par.best_fit
        print "Funcfit:", par.funcfit
        print "Phi max:", t_bins[np.argmax(par.funcfit)]
        par.phi_max = t_bins[np.argmax(par.funcfit)]
        par.phi_min = t_bins[np.argmin(par.funcfit)]
        par.phi_max_err = np.sqrt(np.var(all_max_phase_boot[i,:], ddof=1))
        par.phi_min_err = np.sqrt(np.var(all_min_phase_boot[i,:], ddof=1))
        print "Phi max err:", par.phi_max_err


    # phase_points = phase_points[1:,]
    # return avg_untied, phase_points, model_name

    return avg_untied, boot_model_name


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
    boot_parfit_file = es_dir + "/" + prefix + "_" + day + "_" + \
            FIT_SPECIFIER + "_funcfit.txt"
    ccf_data_file = ccf_dir + "/" + prefix + "_" + day + "_t" + str(dt_mult) + \
            "_" + str(n_seconds) + "sec_adj.fits"
    lag_data_file = lag_dir + "/" + prefix + "_" + day + "_t" + str(dt_mult) + \
            "_" +str(n_seconds) + "sec_adj_lag.fits"
    data_parfit_file = ES_BASEDIR + "/out_es/" + prefix + "/" + prefix + "_" + \
            day + "_" + FIT_SPECIFIER + "-wMCMC_funcfit.txt"
    xspec_log_data_file = ES_BASEDIR + "/out_es/" + prefix + "/" + prefix + \
            "_" + day + "_" + FIT_SPECIFIER + "_xspec.log"
    final_bestfit_file = ES_BASEDIR + "/out_es/" + prefix + "/" + prefix + "_" + \
            day + "_" + FIT_SPECIFIER + "_final_bestfit.txt"

    epoch = 5
    out_rel_name = prefix + "_" + day + "_" + FIT_SPECIFIER
    out_rootname = out_dir + "/" + out_rel_name  # with no extension

    ## Make sure the output directory exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # ## Copy the response matrix into the reduced data directory, if it doesn't
    # ## exist
    # if os.path.isfile(data_dir + "/PCU2.rsp"):
    #     subprocess.call(["cp", data_dir + "/PCU2.rsp", \
    #             out_dir + "/" + rsp_matrix])
    # else:
    #     print("ERROR: Response matrix doesn't exist in the reduced data "\
    #             "directory.")
    #     exit()
    #
    # ## Make sure heasoft is running
    # try:
    #     if "heasoft" not in DYLD_LIB_PATH:
    #         subprocess.Popen(HEAINIT, shell=True)
    # except TypeError:  # This triggers if 'DYLD_LIB_PATH' is 'NoneType'
    #     subprocess.Popen(HEAINIT, shell=True)
    #
    # assert os.path.isfile(boot_parfit_file)
    # assert os.path.isfile(data_parfit_file)
    # assert os.path.isfile(ccf_data_file)
    # assert os.path.isfile(lag_data_file)
    # assert os.path.isfile(xspec_log_data_file)
    #
    # ## Get parameters from the ccf data file header
    # try:
    #     ccf_table = Table.read(ccf_data_file)
    # except IOError:
    #     print("\tERROR: File does not exist or could not be read: %s" \
    #           % ccf_data_file)
    #     exit()
    #
    # obs_time = ccf_table.meta['EXPOSURE']
    # n_seg = ccf_table.meta['SEGMENTS']
    # detchans = ccf_table.meta['DETCHANS']
    # n_bins = ccf_table.meta['N_BINS']
    # n_seconds = ccf_table.meta['SEC_SEG']
    # dt = ccf_table.meta['DT']
    #
    # ## Run fake_qpo_spectra.py
    # if os.path.isfile(data_parfit_file):
    #     fake_qpo_spectra.main(out_rootname, data_parfit_file, dt=dt,
    #             rsp_matrix=rsp_matrix, exposure=obs_time, test=testing,
    #             n_params=N_PARAMS)
    #     print "Finished fake_qpo_spectra"
    #
    # else:
    #     print("\t ERROR: Parameter fit file does not exist: %s" %
    #           data_parfit_file)
    #     exit()

    ## Read in values of all parameters from the XSPEC log files
    untied_par_values, n_untied, good_boots = read_xspec_log_files(es_dir,
            out_rel_name, boot_num=boot_num)
    avg_untied = get_boot_varpar_err(untied_par_values,
            xspec_log_data_file)

    ## Compute the average fit function to the untied parameter variations
    ## Using boot_num-bad_boot as boot_num, because we're not counting the bad
    ## bootstrap runs.
    avg_untied, model_name = get_best_fit(avg_untied,
            boot_parfit_file, n_untied)

    write_final_besfit_pars(final_bestfit_file, model_name, avg_untied)
    print final_bestfit_file

    model_name = model_name.upper()
    model_name = model_name.replace("PHABS*", "PHABS$\\times\\,$")
    model_name = model_name.replace("*", " * ")
    model_name = model_name.replace("+", " + ")

    # print model_name
    # np.savetxt("%s_phasepoints.txt" % FIT_SPECIFIER, phase_points)

    ## Plot the "multi-fit" parameter variations to the average values
    plot_final_multifit(avg_untied, out_rootname=out_rootname,
            title=model_name)

    ## Make contour plots of the phases
    # plot_phase_contours(phase_points, avg_untied, out_rootname=out_rootname,
    #         title=model_name)
    # plot_deltaphase_hist(phase_points, out_rootname=out_rootname,
    #         title=model_name)

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
