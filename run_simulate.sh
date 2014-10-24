#!/bin/bash

home_dir=$(ls -d ~)  # the -d flag is extremely important here
day=$(date +%y%m%d)  # make the date a string and assign it to 'day'
exe_dir="$home_dir/Dropbox/Research/simulate"
out_dir="$exe_dir/out_sim"

freq=401.0
# bb_spec="spectra/100000s_mean.fak"
# pl_spec="spectra/100000s_mean.fak"
bb_spec="spectra/100000s_bb_nopoiss.fak"
pl_spec="spectra/100000s_pl_nopoiss.fak"
amp_ci=0.065
amp_ref=0.065
exposure=2400.0
phase_spec=0.0
num_sec=4

if [ ${bb_spec: -4} == ".fak" ] && [ ${pl_spec: -4} == ".fak" ]; then
	bb_exposure=$(python -c "from tools import get_key_val; print get_key_val('$bb_spec', 1, 'EXPOSURE')")
	pl_exposure=$(python -c "from tools import get_key_val; print get_key_val('$pl_spec', 1, 'EXPOSURE')")
	if [ $bb_exposure == $pl_exposure ]; then
		exposure=$bb_exposure
	fi
fi

echo "$exposure"

# python "$exe_dir"/simulate_lightcurves.py -h

time python "$exe_dir"/simulate_lightcurves.py --freq "$freq" --bb "$bb_spec" --pl "$pl_spec" --amp_ci "$amp_ci" --amp_ref "$amp_ref" --num_seconds "$num_sec" --exposure "$exposure" --phase_spec "$phase_spec" #--test

# if [ -e "$exe_dir/plot.png" ]; then
# 	open -a ImageJ "$exe_dir/plot.png"
# fi

if [ -e "$exe_dir/sim_power.png" ]; then
	open -a ImageJ "$exe_dir/sim_power.png"
fi