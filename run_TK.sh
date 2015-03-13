#!/bin/bash

home_dir=$(ls -d ~)
sim_dir="$home_dir/Dropbox/Research/simulate"
pow_dir="$home_dir/Dropbox/Research/power_spectra"

prefix="FAKE-TK-GX339B"
day=$(date +%y%m%d)  # make the date a string and assign it to 'day'

nbins=8192
dt=0.0078125
rebin_const=1.01
meanrate=2331.0
numsegments=267
variance=4.61214030908e-07

cd "$sim_dir"

python "$sim_dir"/TimmerKoenig.py "$nbins" "$dt" "$meanrate" "$numsegments" "$variance"

# python "$pow_dir"/plot_powerspec.py "$sim_dir"/TK_power.fits \
# 		-o "$sim_dir"/TK_psd.png -p "$prefix"

python "$pow_dir"/rebin_powerspec.py "$sim_dir"/TK_power.fits \
		"$sim_dir"/TK_power_rb.fits -o "$sim_dir"/TK_psd_rb.png -p "$prefix" \
		-c "$rebin_const"

if [ -e "$sim_dir/TK_psd_rb.png" ]; then open "$sim_dir/TK_psd_rb.png"; fi