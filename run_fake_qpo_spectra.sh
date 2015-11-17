#!/usr/bin/env bash

################################################################################
##
## Bash script for making fake lag-energy spectra given the energy spectral
## parameters from multifit_plots.py in sed_fitting.sh
##
## Example call: ./run_fake_qpo_spectra.sh GX339-BQPO 64 64 0 150131
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: HEASOFT 6.11.*, bash 3.*, and Python 2.7.* (with supporting libraries)
## 		  must be installed in order to run this script.
##
## Written by Abigail Stevens, A.L.Stevens at uva.nl, 2015
##
################################################################################

## Checking the number of input arguments
if (( $# != 5 )); then
    echo -e "\tUsage: ./run_fake_qpo_spectra.sh <prefix> <dt multiple> <num "\
        "seconds> <testing> <date>\n"
    exit
fi

prefix=$1
dt_mult=$2
numsec=$3
testing=$4
day=$5

################################################################################

## If heainit isn't running, start it
if (( $(echo $DYLD_LIBRARY_PATH | grep heasoft | wc -l) < 1 )); then
	. $HEADAS/headas-init.sh
fi

home_dir=$(ls -d ~)
ccf_dir="$home_dir/Dropbox/Research/cross_correlation/out_ccf/${prefix}"
es_dir="$home_dir/Dropbox/Research/energy_spectra/out_es/${prefix}"
exe_dir="$home_dir/Dropbox/Research/simulate"
out_dir="$exe_dir/out_sim/${prefix}"
data_dir="$home_dir/Reduced_data/$prefix"
rsp_matrix="${prefix}_PCU2.rsp"
dump_file=dump.txt # Name of dumping file for intermediary steps

#fit_specifier="1BB-Nicer"
#fit_specifier="1BB-FS-G-Tin-fzs-fzNbb"
#fit_specifier="2BB-FS-G-kT"
echo Fit specifier: "$fit_specifier"
#n_spectra=24
#n_params=9
#echo "$n_spectra"
#echo "$n_params"

parfit_file="$es_dir/${prefix}_${day}_${fit_specifier}_sines.txt"
ccf_file="$ccf_dir/${prefix}_${day}_t${dt_mult}_${numsec}sec_adj.fits"

#prefix="FAKE-NICER"
#prefix="FAKE-${prefix}-Nicer"
out_name="${prefix}_${day}_${fit_specifier}"
out_root="${out_dir}/${out_name}"

plot_ext="eps"

################################################################################
################################################################################

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi
if [ -e "$data_dir/PCU2.rsp" ]; then
    cp "$data_dir/PCU2.rsp" "$out_dir/$rsp_matrix"
else
    echo "ERROR: Response matrix doesn't exist in the reduced data directory."
    exit
fi

#if [ ! -e "$ccf_file" ]; then echo -e "\tCCF file does not exist."; exit; fi
#obs_time=$(python -c "from tools import get_key_val; print get_key_val('$ccf_file', 0, 'EXPOSURE')")
#dt=$(python -c "from tools import get_key_val; print get_key_val('$ccf_file', 0, 'DT')")
#n_seg=$(python -c "from tools import get_key_val; print get_key_val('$ccf_file', 0, 'SEGMENTS')")
#detchans=$(python -c "from tools import get_key_val; print get_key_val('$ccf_file', 0, 'DETCHANS')")
#n_bins=$(python -c "from tools import get_key_val; print get_key_val('$ccf_file', 0, 'N_BINS')")

obs_time=13224.3984375
#obs_time=640
#dt=0.0078125
dt=0.008153062878233013
n_seg=198
#n_seg=10
detchans=64
#detchans=238
n_bins=8192
epoch=5

if [ -e "${parfit_file}" ]; then
    python "$exe_dir"/fake_qpo_spectra.py "${out_root}" "${parfit_file}" \
            --prefix "${prefix}" --n_bins "${n_bins}" --dt "${dt}" \
            --n_seg "${n_seg}" --n_spec "${n_spectra}" --n_par "${n_params}" \
            --chan "${detchans}" --epoch "${epoch}" --exposure "${obs_time}" \
            --test "${testing}"  --rsp "${rsp_matrix}"
else
    echo -e "\t ERROR: Parameter fit file does not exist."
fi
