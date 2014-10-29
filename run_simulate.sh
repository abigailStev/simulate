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
numsec=4

if [ ${bb_spec: -4} == ".fak" ] && [ ${pl_spec: -4} == ".fak" ]; then
	bb_exposure=$(python -c "from tools import get_key_val; print get_key_val('$bb_spec', 1, 'EXPOSURE')")
	pl_exposure=$(python -c "from tools import get_key_val; print get_key_val('$pl_spec', 1, 'EXPOSURE')")
	if [ $bb_exposure == $pl_exposure ]; then
		exposure=$bb_exposure
	fi
fi

# echo "$exposure"

# python "$exe_dir"/simulate_lightcurves.py -h

# time python "$exe_dir"/simulate_lightcurves.py --freq "$freq" --bb "$bb_spec" --pl "$pl_spec" --amp_ci "$amp_ci" --amp_ref "$amp_ref" --num_seconds "$numsec" --exposure "$exposure" --phase_spec "$phase_spec" #--test
# 
# # if [ -e "$exe_dir/plot.png" ]; then
# # 	open -a ImageJ "$exe_dir/plot.png"
# # fi
# 
# if [ -e "$exe_dir/sim_power.png" ]; then
# 	open -a ImageJ "$exe_dir/sim_power.png"
# fi




time python ccf_simulation.py

ccf_file="$out_dir/ccf_out.dat"
ensp_rsp_matrix="$home_dir/Dropbox/Research/energy_spectra/out_es/P70080_141029_70080-01-01-02_PCU2.rsp"
rsp_matrix="$out_dir/P70080_141029_70080-01-01-02_PCU2.rsp"
cp "$ensp_rsp_matrix" "$rsp_matrix"
propID="P70080"
day="141029"
dt=1
spec_type=1
tab_ext="dat"
ensp_exe_dir="$home_dir/Dropbox/Research/energy_spectra"
ensp_out_dir="out_es"
dump_file="dum.dat"


if [ -e "$ccf_file" ]; then
	obs_time=$(python -c "from tools import read_obs_time; print read_obs_time('$ccf_file')")
	echo "PIPELINE EXPOSURE TIME =" $obs_time "s"
# 	echo `echo "$obs_time / 16.0" | bc -l`
else
	obs_time = 0
	echo -e "\n\t Couldn't get observation time from header, set obs_time to zero.\n"
fi



## Generating energy spectra at each phase bin
# for tbin in {25..40..5}; do  ## should work in bash 4.*, but i have 3.2.*
for (( tbin=25; tbin<=40; tbin+=5 )); do
# # 	echo "$tbin"
	out_end="${day}_t${dt}_${numsec}sec_pbin_${tbin}"
	out_end_mean="${day}_t${dt}_${numsec}sec_mean"
	out_file="$out_dir/$out_end"
	echo "$out_end"
	echo "$out_end_mean"
	echo "$out_file"
	
	if [ -e "${ccf_file}" ]; then
		python "$ensp_exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_file}.${tab_ext}" -b "$tbin" -s "$spec_type"
	else
		echo -e "\n\t ${ccf_file} does not exist, energyspec.py was NOT run.\n"
	fi
	
	if [ -e "$rsp_matrix" ] && [ -e "${out_file}.${tab_ext}" ]; then
		
		cd "$out_dir"
		ascii2pha infile="${out_end}.${tab_ext}" \
			outfile="${out_end}.pha" \
			chanpres=yes \
			dtype=2 \
			qerror=yes \
			rows=- \
			tlmin=0 \
			detchans=64 \
			pois=no \
			telescope=RXTE \
			instrume=PCA \
			detnam=PCU2 \
			filter=NONE \
			exposure=$obs_time \
			clobber=yes \
			respfile="$rsp_matrix" > $dump_file
		echo "XSPEC data: ${out_file}.pha"
		echo -e "XSPEC resp: $rsp_matrix\n"
	else
		echo -e "\n\t${rsp_matrix} and/or ${out_file}.${tab_ext} do NOT exist, ascii2pha was NOT run.\n"
	fi
# 	if [ ! -e "${out_end}.pha" ]; then
# 		echo -e "\n\tERROR: ASCII2pha didn't run, ${out_end}.pha does not exist.\n"
# 	fi
done