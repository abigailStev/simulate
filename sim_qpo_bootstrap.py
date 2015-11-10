#!/usr/bin/env python
"""
Makes fake lag-energy spectra given the energy spectral parameters from
multifit_plots.py in sed_fit_bootstrap.sh. Calls fake_qpo_spectra.py.
Used in ccf_bootstrap.sh.
"""
__author__ = 'Abigail Stevens <A.L.Stevens at uva.nl>'

import argparse
import os.path
import tools
import subprocess
import fake_qpo_spectra
from datetime import datetime

HOME_DIR = os.path.expanduser("~")
CCF_BASEDIR = HOME_DIR + "/Dropbox/Research/cross_correlation"
ES_BASEDIR = HOME_DIR + "/Dropbox/Research/energy_spectra"
EXE_DIR = HOME_DIR + "/Dropbox/Research/simulate"
DATA_BASEDIR = HOME_DIR + "/Reduced_data"
DUMP_FILE = "dump.txt"  # Name of dumping file for intermediary steps

N_PARAMS = os.getenv('n_params', 9)
N_SPECTRA = os.getenv('n_spectra', 24)
FIT_SPECIFIER = os.getenv('fit_specifier', "1BB-FS-G-Tin-fzs-fzNbb")

DYLD_LIB_PATH = os.environ.get('DYLD_LIBRARY_PATH')
HEADAS_DIR = os.environ.get('HEADAS')
HEAINIT = ". " + HEADAS_DIR + "/headas-init.sh"

TAB_EXT = "fits"
PLOT_EXT = "eps"

def main(prefix="--", dt_mult=64, n_seconds=64, boot_num=0, testing=False,
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
        Number of bootstrap realizations to do/done. [0]

    testing : bool
        If present, sets to True. Used for testing script so it'll only run one
        segment. [False]

    day : str
        Date of the analysis, in YYMMDD format. [default is today's date]

    Returns
    -------
    nothing
    """
    ccf_dir = CCF_BASEDIR + "/out_ccf/" + prefix + "/bootstrapped"
    es_dir = ES_BASEDIR + "/out_es/" + prefix + "/bootstrapped"
    out_dir = EXE_DIR + "/out_sim/" + prefix + "/bootstrapped"
    data_dir = DATA_BASEDIR + "/" + prefix
    rsp_matrix = prefix + "_PCU2.rsp"
    parfit_file = es_dir + "/" + prefix + "_" + day + "_" + FIT_SPECIFIER + \
                  "_funcfit.txt"
    ccf_file = ccf_dir + "/" + prefix + "_" + day + "_t" + str(dt_mult) + "_" +\
               str(n_seconds) + "sec_adj.fits"

    # print parfit_file
    # print ccf_file

    out_rel_name = prefix + "_" + day + "_" + FIT_SPECIFIER
    out_rootname = out_dir + "/" + out_rel_name  # with no extension

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if os.path.isfile(data_dir+"/PCU2.rsp"):
        subprocess.call(["cp", data_dir+"/PCU2.rsp", out_dir+"/"+rsp_matrix])
    else:
        print("ERROR: Response matrix doesn't exist in the reduced data directory.")
        exit()

    if "heasoft" not in DYLD_LIB_PATH:
        subprocess.Popen(HEAINIT, shell=True)

    obs_time = tools.get_key_val(ccf_file, 0, 'EXPOSURE')
    dt = tools.get_key_val(ccf_file, 0, 'DT')
    # n_seg = tools.get_key_val(ccf_file, 0, 'SEGMENTS')
    # detchans = tools.get_key_val(ccf_file, 0, 'DETCHANS')
    # n_bins = tools.get_key_val(ccf_file, 0, 'N_BINS')

    # obs_time = 13224.3984375
    # dt = 0.008153062878233013
    # n_seg = 198
    # n_bins = 8192
    # epoch = 5

    if os.path.isfile(parfit_file):
        fake_qpo_spectra.main(out_rootname, parfit_file, dt=dt, \
                rsp_matrix=rsp_matrix, exposure=obs_time, test=testing)
    else:
        print("\t ERROR: Parameter fit file does not exist: %s" % parfit_file)
        exit()


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

    parser.add_argument("--testing", default=False, dest="test",
            action="store_true", help="If present, sets to True. Used for "\
            "testing script so it'll only run one segment. [False]")

    parser.add_argument("--day", default=datetime.now().strftime("%y%m%d"),
            dest="day", help="Date of the analysis, in YYMMDD format. "\
            "[%s]" % datetime.now().strftime("%y%m%d"))

    parser.add_argument("--boot", default=10, dest="boot_num", type=int,
            help="Number of bootstrap realizations done. [10]")

    args = parser.parse_args()

    main(prefix=args.prefix, dt_mult=args.dt_mult, n_seconds=args.n_seconds,
            day=args.day, boot_num=args.boot_num, testing=args.test)