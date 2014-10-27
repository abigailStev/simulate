import argparse
import ccf as crosscorr
import powerspec
import tools
import simulate_lightcurves as sim_lc
import numpy as np
import plot_ccf


###############################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Simulates a light curve and \
		runs it through ccf.py.')
	parser.add_argument('--test', action='store_true', dest='test', help='If \
		present, only does a short test run.')
	args = parser.parse_args()

	freq = 401.0
	bb_spec = "spectra/100000s_bb_nopoiss.fak"
	pl_spec = "spectra/100000s_pl_nopoiss.fak"
	amp_ci = 0.065
	amp_ref = 0.065
	phase_spec = 0.0
	num_seconds = 4
	dt_mult = 1
	exposure = 100000.0
	mean_ci = 1.0
	mean_ref = 1.0
	phase_ci = 0.0
	phase_spec = 0.0

	t_res = 1.0 / 8192.0  # The time resolution of the data, in seconds.
	dt = dt_mult * t_res
	n_bins = num_seconds * int(1.0 / dt)
	assert tools.power_of_two(n_bins)  # n_bins must be a power of 2 for the FFT
	
	spec_bb = sim_lc.read_fakeit_spectra(bb_spec)
	spec_pl = sim_lc.read_fakeit_spectra(pl_spec)
	mean_curve_ci = 0
	mean_curve_ref = 0
	mean_ps_ref = 0
	mean_rate_ref = 0
	
	mean_rate_whole_ci = np.zeros(64, dtype=np.float64)
	mean_rate_whole_ref = 0
	sum_rate_whole_ci = np.zeros(64, dtype=np.float64)
	sum_rate_whole_ref = 0
	sum_power_ci = np.zeros((n_bins, 64), dtype=np.float64)
	sum_power_ref = np.zeros(n_bins, dtype=np.float64)
	ccf_filtered = np.zeros((n_bins, 64))
	cs_sum = np.zeros((n_bins, 64), dtype=np.complex128)

    
    
	## Making the signals
	sine_ci, sine_ref = sim_lc.generate_sines(dt, n_bins, freq, amp_ci, \
		amp_ref, mean_ci, mean_ref, phase_ci, phase_spec)
	curve_ci_bb, curve_ref_bb = sim_lc.make_lightcurves(spec_bb, sine_ci, \
		sine_ref, n_bins)
	curve_ci_pl, curve_ref_pl = sim_lc.make_lightcurves(spec_pl, sine_ci, \
		sine_ref, n_bins)
	
		
	############################
	## Looping through segments
	############################
	for num_segments in xrange(1, 601): # 'i' tracks the number of segments
		curve_ci, curve_ref = sim_lc.add_lightcurves(curve_ci_bb, curve_ref_bb,\
			curve_ci_pl, curve_ref_pl, dt, exposure)
		mean_curve_ci += curve_ci
		mean_curve_ref += curve_ref
		
		ps_ref, rate_ref = powerspec.make_ps(curve_ref)
		mean_ps_ref += ps_ref
		mean_rate_ref += rate_ref
		
# 		print "Shape mean curve ci:", np.shape(mean_curve_ci)
# 		print "Shape mean curve ref:", np.shape(mean_curve_ref)
		
		mean_rate_segment_ci = np.zeros(64, dtype=np.float64)
		mean_rate_segment_ref = np.zeros(64, dtype=np.float64)
		cs_segment = np.zeros((n_bins, 64), dtype=np.complex128)

		cs_segment, mean_rate_segment_ci, mean_rate_segment_ref, power_ci, power_ref = crosscorr.make_cs(curve_ci, curve_ref, n_bins, dt)
		
		sum_rate_whole_ci += mean_rate_segment_ci
		sum_rate_whole_ref += mean_rate_segment_ref
		sum_power_ci += power_ci
		sum_power_ref += power_ref
		cs_sum += cs_segment  # This adds indices
		
# 		sum_rate_ci += np.mean(curve.ci)
		
		if num_segments == 1 and args.test == True: 
			break
		
		if num_segments % 10 == 0: print "\t", num_segments
	## End of for-loop
	
	## Averages
	mean_ps_ref /= float(num_segments)
	mean_rate_ref /= float(num_segments)
	mean_curve_ci /= float(num_segments)
	mean_curve_ref /= float(num_segments)
	mean_rate_whole_ci = sum_rate_whole_ci / float(num_segments)
	mean_rate_whole_ref = sum_rate_whole_ref / float(num_segments)
	mean_power_ci = sum_power_ci / float(num_segments)
	mean_power_ref = sum_power_ref / float(num_segments)
	cs_avg = cs_sum / float(num_segments)
	
	print "Mean count rate in reference band:", np.mean(mean_curve_ref)

	ccf_filtered, ccf_error = crosscorr.cs_to_ccf_w_err(cs_avg, dt, n_bins, \
		num_seconds, num_segments, mean_rate_whole_ci, mean_rate_whole_ref, \
		mean_power_ci, mean_power_ref)
		
	t = np.arange(0, n_bins)
	in_file = "Faked: " + bb_spec + " + " + pl_spec
	out_file = "ccf_out.dat"
	crosscorr.output(out_file, in_file, dt, n_bins, num_seconds, num_segments,
        mean_rate_whole_ci, mean_rate_whole_ref, t, ccf_filtered, ccf_error)
	
# 	sim_lc.plot_curves(n_bins, mean_curve_ci[:,6], mean_curve_ref, "plot.png")

	sim_lc.power_spectra_things(mean_ps_ref, dt, n_bins, num_seconds, \
		num_segments, mean_rate_ref)
	
	
	plot_root = "./ccf"
	propID = "FAKE"
	plot_ccf.main(out_file, plot_root, propID)
	
