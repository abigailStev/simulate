#!/bin/bash

home_dir=$(ls -d ~)
sim_dir="$home_dir/Dropbox/Research/simulate"
pow_dir="$home_dir/Dropbox/Research/power_spectra"
ccf_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$sim_dir/out_sim"

prefix="FAKE-TK-GX339B"
# day=$(date +%y%m%d)  # make the date a string and assign it to 'day'

nbins=8192
dt=0.0078125
rebin_const=1.01
meanrate=1000.0  # for one pcu
# meanrate=0
# numsegments=267
numsegments=10
# variance=4.61214030908e-07
# variance=2.1273510648e-07
variance=1e-6
# variance=1
pl_scale=3.29617761e-04
qpo_scale=2.06908556e-02
# pl_scale=0.00033744
# qpo_scale=0.02769814
fake_e_spec="$out_dir/${prefix}.fak"


cd "$sim_dir"

for (( i=0; i<1; i++ )); do
	pow_out="$out_dir/TK_power_${i}.fits"
	ccf_out="$out_dir/TK_ccf_${i}.fits"
	pow_rb_out="$out_dir/TK_power_rb_${i}.fits"
	pow_plot="$out_dir/TK_psd_${i}.png"
	pow_rb_plot="$out_dir/TK_psd_rb_${i}.png"
	ccf_plot="$out_dir/TK_ccf_${i}"
	ccf_2D_plot="$out_dir/TK_ccf_2D_${i}.png"
	
# 	echo "$i"
python "$sim_dir"/TimmerKoenig.py "$nbins" "$dt" "$meanrate" \
		"$numsegments" "$variance" "$pl_scale" "$qpo_scale" "$fake_e_spec" \
		"$pow_out" "$ccf_out"

	python "$pow_dir"/plot_powerspec.py "$pow_out" -o "$pow_plot" -p "$prefix"
	# if [ -e "$pow_plot" ]; then open "$pow_plot"; fi

	python "$pow_dir"/rebin_powerspec.py "$pow_out" "$pow_rb_out" -o "$pow_rb_plot"\
		-p "$prefix" -c "$rebin_const"
	if [ -e "$pow_rb_plot" ]; then open "$pow_rb_plot"; fi

	if [ -e "$ccf_out" ]; then
		python "$ccf_dir"/plot_ccf.py "$ccf_out" -o "$ccf_plot" -p "$prefix" 
	# 	if [ -e "$ccf_plot"_chan_06.png ]; then open "$ccf_plot"_chan_06.png; fi
		python "$ccf_dir"/plot_2d.py "$ccf_out" -o "$ccf_2D_plot" -p "${prefix}"
		if [ -e "$ccf_2D_plot" ]; then open "$ccf_2D_plot"; fi
	fi
done