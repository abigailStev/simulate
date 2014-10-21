#!/bin/bash

home_dir=$(ls -d ~)  # the -d flag is extremely important here
day=$(date +%y%m%d)  # make the date a string and assign it to 'day'
exe_dir="$home_dir/Dropbox/Research/simulate"
out_dir="$exe_dir/out_sim"

freq=401.0
bb_spec="spectra/fakeit_mean.fak"
pl_spec="spectra/fakeit_mean.fak"
exposure=2
phase_ci=-1.0  # Negative means ci lags ref, positive means ref lags ci



python "$exe_dir"/simulate_lightcurves.py -h

python "$exe_dir"/simulate_lightcurves.py --freq "$freq" --bb "$bb_spec" --pl "$pl_spec" --exposure "$exposure" --phase_ci "$phase_ci"