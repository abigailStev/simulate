#!/bin/bash

home_dir=$(ls -d ~)  # the -d flag is extremely important here
day=$(date +%y%m%d)  # make the date a string and assign it to 'day'
# day=141103
exe_dir="$home_dir/Dropbox/Research/simulate"
out_dir="$exe_dir/out_sim"
ccf_dir="$home_dir/Dropbox/Research/cross_correlation"


freq=401.0
# bb_spec="spectra/100000s_mean.fak"
# pl_spec="spectra/100000s_mean.fak"
bb_spec="$exe_dir/spectra/100000s_bb_nopoiss.fak"
pl_spec="$exe_dir/spectra/100000s_pl_nopoiss.fak"
amp_ci=0.055
amp_ref=0.055
exposure=100000.0
numsec=4
obs_time=0  # this is set down below
dt=1

phase=45  # degrees
# phase=180    # degrees
# phase=0
out_root="FAKE_${day}_t${dt}_${numsec}sec_phase${phase}"
ccf_file="$out_dir/${out_root}_ccf.fits"
plot_root="$out_dir/${out_root}_ccf"
ccfs_plot="$out_dir/${out_root}_manyccfs.png"


#################################
## Running simulate_lightcurves
#################################

# python "$exe_dir"/simulate_lightcurves.py -h

# time python "$exe_dir"/simulate_lightcurves.py --freq "$freq" --bb "$bb_spec" --pl "$pl_spec" --amp_ci "$amp_ci" --amp_ref "$amp_ref" --num_seconds "$numsec" --exposure "$exposure" --phase "$phase" #--test
# 
# if [ -e "$exe_dir/plot.png" ]; then
# 	open -a ImageJ "$exe_dir/plot.png"
# fi
# if [ -e "$exe_dir/sim_power.png" ]; then
# 	open -a ImageJ "$exe_dir/sim_power.png"
# fi


###########################
## Running ccf_simulation
###########################

time python "$exe_dir/"ccf_simulation.py "$ccf_file" --bb "$bb_spec" --pl "$pl_spec" --freq "$freq" --amp_ci "$amp_ci" --amp_ref "$amp_ref" --num_seconds "$numsec" --phase "$phase" --test # --noisy 

if [ -e "$ccf_file" ]; then
	python "$ccf_dir/"plot_ccf.py "$ccf_file" -o "$plot_root" -p "FAKE"
	open -a ImageJ "${plot_root}_chan_06.png"
	python "$ccf_dir"/plot_multi.py "$ccf_file" "$ccfs_plot" "${numsec}"
	open -a ImageJ "$ccfs_plot"
fi

## If heainit isn't running, start it
if (( $(echo $DYLD_LIBRARY_PATH | grep heasoft | wc -l) < 1 )); then
	echo -e "\nStarting heainit"
	. $HEADAS/headas-init.sh
fi

ensp_rsp_matrix="$home_dir/Dropbox/Research/energy_spectra/out_es/P70080_141029_70080-01-01-02_PCU2.rsp"
rsp_matrix="$out_dir/P70080_141029_70080-01-01-02_PCU2.rsp"
cp "$ensp_rsp_matrix" "$rsp_matrix"
propID="FAKE"
tab_ext="dat"
ensp_exe_dir="$home_dir/Dropbox/Research/energy_spectra"
ensp_out_dir="out_es"
dump_file="dum.dat"



if [ -e "$ccf_file" ]; then
	if [ "${ccf_file##*.}" == "dat" ]; then
		obs_time=$(python -c "from tools import read_obs_time; print read_obs_time('$ccf_file')")
	elif [ "${ccf_file##*.}" == "fits" ]; then
		obs_time=$(python -c "from tools import get_key_val;  print get_key_val('$ccf_file', 0, 'EXPOSURE')")
	fi
	echo "EXPOSURE TIME =" $obs_time "s"
else
	obs_time=0
	echo -e "\tERROR: Couldn't get observation time from header, set obs_time to zero."
	exit
fi

if [ $(echo " $obs_time == 0.0" | bc) -eq 1 ]; then
	echo -e "\tERROR: Exposure time is zero. Exiting script."
	exit
fi


## Making ccf only spectra and plotting them

spec_type=1  # 0 for mean+ccf, 1 for ccf, 2 for mean
xspec_script="$out_dir/${out_root}_xspec.xcm"
spectrum_plot="${out_root}_spectrum_ccf"

if [ -e "$xspec_script" ]; then
	rm "$xspec_script"
fi
touch "$xspec_script"
i=1
cd "$out_dir"
## Generating energy spectra at each phase bin
# for tbin in {25..40..5}; do  ## should work in bash 4.*, but i have 3.2.*
for (( tbin=25; tbin<=40; tbin+=5 )); do
# # 	echo "$tbin"
	out_end="${out_root}_ccf_${tbin}bin"
	out_end_mean="${out_root}_mean"
# 	out_file="$out_dir/$out_end"
# 	echo "$out_end"
# 	echo "$out_end_mean"
# 	echo "$out_file"
	
	if [ -e "${ccf_file}" ]; then
		python "$ensp_exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_end}.${tab_ext}" -b "$tbin" -s "$spec_type"
	else
		echo -e "\t${ccf_file} does NOT exist, energyspec.py was NOT run."
	fi
	
	if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ]; then
		
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
# 		echo "XSPEC data: ${out_file}.pha"
# 		echo -e "XSPEC resp: $rsp_matrix\n"
	else
		echo -e "\t${rsp_matrix} and/or ${out_end}.${tab_ext} do NOT exist, ascii2pha was NOT run."
	fi
# 	if [ ! -e "${out_end}.pha" ]; then
# 		echo -e "\tASCII2pha didn't run, so ${out_end}.pha does not exist."
# 	fi
	echo "data $i:$i $out_end" >> $xspec_script
	((i+=1))
done

echo "ignore 1-4: **-3 11 27-**" >> $xspec_script
echo "notice 1-4: 3 27" >> $xspec_script
echo "cpd /xw" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "mod pow & 0" >> $xspec_script
echo "iplot eeufspec" >> $xspec_script
echo "@ccf_nomean.pco $spectrum_plot" >> $xspec_script
echo "exit" >> $xspec_script
cd "$out_dir"
xspec < "$xspec_script" > "$dump_file"
# echo "$spectrum_plot.eps"
if [ -e "$spectrum_plot.eps" ]; then
	awk '/%%EndComments/{print "%% phase = "'$phase'" degrees, of PL wrt BB"}1' $spectrum_plot.eps > temp.eps
	mv -f temp.eps "$spectrum_plot.eps"
# 	open -a "TextWrangler" "$spectrum_plot.eps"
	open "$spectrum_plot.eps"
	echo "Spectrum: $spectrum_plot.eps"
fi


## Making ccf+mean spectra and plotting them
spec_type=0  # 0 for mean+ccf, 1 for ccf, 2 for mean
xspec_script="$out_dir/${out_root}_xspec.xcm"
spectrum_plot="${out_root}_spectrum_ccfwmean"

if [ -e "$xspec_script" ]; then
	rm "$xspec_script"
fi
touch "$xspec_script"
i=1
cd "$out_dir"
## Generating energy spectra at each phase bin
# for tbin in {25..40..5}; do  ## should work in bash 4.*, but i have 3.2.*
for (( tbin=25; tbin<=40; tbin+=5 )); do
# # 	echo "$tbin"
	out_end="${out_root}_ccfwmean_${tbin}bin"
	out_end_mean="${out_root}_mean"
# 	out_file="$out_dir/$out_end"
# 	echo "$out_end"
# 	echo "$out_end_mean"
# 	echo "$out_file"

	if [ -e "${ccf_file}" ]; then
		python "$ensp_exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_end}.${tab_ext}" -b "$tbin" -s "$spec_type"
	else
		echo -e "\t${ccf_file} does not exist, energyspec.py was NOT run."
	fi
	
	if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ]; then
		
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
# 		echo "XSPEC data: ${out_end}.pha"
# 		echo -e "XSPEC resp: $rsp_matrix\n"
	else
		echo -e "\t${rsp_matrix} and/or ${out_end}.${tab_ext} do NOT exist, ascii2pha was NOT run."
	fi
# 	if [ ! -e "${out_end}.pha" ]; then
# 		echo -e "\n\tERROR: ASCII2pha didn't run, ${out_end}.pha does not exist.\n"
# 	fi
	echo "data $i:$i $out_end" >> $xspec_script
	((i+=1))
done

echo "ignore 1-4: **-3 11 27-**" >> $xspec_script
echo "notice 1-4: 3 27" >> $xspec_script
echo "cpd /xw" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "mod pow & 0" >> $xspec_script
echo "iplot eeufspec" >> $xspec_script
echo "@mean.pco $spectrum_plot" >> $xspec_script
echo "exit" >> $xspec_script
cd "$out_dir"
xspec < "$xspec_script" > "$dump_file"
if [ -e "$spectrum_plot.eps" ]; then
	awk '/%%EndComments/{print "%% phase = "'$phase'" degrees, of PL wrt BB"}1' $spectrum_plot.eps > temp.eps
	mv -f temp.eps "$spectrum_plot.eps"
# 	open -a "TextWrangler" "$spectrum_plot.eps"
	open "$spectrum_plot.eps"
	echo "Spectrum: $spectrum_plot.eps"
fi