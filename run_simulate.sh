#!/bin/bash

home_dir=$(ls -d ~)  # the -d flag is extremely important here
day=$(date +%y%m%d)  # make the date a string and assign it to 'day'
exe_dir="$home_dir/Dropbox/Research/simulate"
out_dir="$exe_dir/out_sim"

freq=401.0
bb_spec="spectra/fakeit_mean.fak"
pl_spec="spectra/fakeit_mean.fak"
amp_ci=0.09
amp_ref=0.09
exposure=1.0
phase_spec=-1.0
num_sec=1

if [ ${bb_spec: -4} == ".fak" ] && [ ${pl_spec: -4} == ".fak" ]; then
	bb_exposure=$(python -c "from tools import get_key_val; print get_key_val('$bb_spec', 1, 'EXPOSURE')")
	pl_exposure=$(python -c "from tools import get_key_val; print get_key_val('$pl_spec', 1, 'EXPOSURE')")
	if [ $bb_exposure == $pl_exposure ]; then
		exposure=$bb_exposure
	fi
fi

# python "$exe_dir"/simulate_lightcurves.py -h

python "$exe_dir"/simulate_lightcurves.py --freq "$freq" --bb "$bb_spec" --pl "$pl_spec" --amp_ci "$amp_ci" --amp_ref "$amp_ref" --num_seconds "$num_sec" --exposure "$exposure" --phase_spec "$phase_spec"

if [ -e "$exe_dir/plot.png" ]; then
	open -a ImageJ "$exe_dir/plot.png"
fi

if [ -e "$exe_dir/powerspectrum.png" ]; then
	open -a ImageJ "$exe_dir/powerspectrum.png"
fi