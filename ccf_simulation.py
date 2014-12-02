import argparse
import numpy as np
from datetime import datetime
from astropy.io import fits
import os

import ccf as crosscorr
import powerspec
import tools
import simulate_lightcurves as sim_lc
import plot_ccf

"""
		ccf_simulation.py

Simulates two light curves from fake spectra and computes their cross-
correlation function. 

Required arguments:
out_file - Name of (ASCII/txt/dat) output file which the table of
    cross-correlation function data will be written to.
    
Optional arguments:
test - If present, does a single run. If not, does a full run. 

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as python 2.7, can be downloaded
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

'tools', 'simulate_lightcurves', 'plot_ccf', 'ccf', and 'powerspec' are available on GitHub. 

"""
###############################################################################
def fits_output(out_file, bb_spec, pl_spec, freq, phase, dt, n_bins, \
	num_seconds, num_segments, amp_ci, amp_ref, mean_ci, mean_ref, \
	mean_rate_whole_ci, mean_rate_whole_ref, t, ccf_filtered, ccf_error, noisy):
	"""
			fits_output
	
	
	
	"""
	print "Output sent to %s" % out_file
	## Making header
	prihdr = fits.Header()
	prihdr.set('TYPE', "Cross-correlation of simulation")
	prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
	prihdr.set('BB_SPEC', bb_spec, "Blackbody (bb) spectrum file")
	prihdr.set('PL_SPEC', pl_spec, "Power law (pl) spectrum file")
	prihdr.set('SIG_FREQ', freq, "Hz; Frequency of simulated periodic signal")
	prihdr.set('PHASE', phase, "radians; phase shift of pl from bb")
	prihdr.set('AMP_CI', amp_ci, "Amplitude of signal in channels of interest")
	prihdr.set('AMP_REF', amp_ref, "Amplitude of signal in reference band")
	prihdr.set('MEAN_CI', mean_ci, "Mean of signal in channels of interest")
	prihdr.set('MEAN_REF', mean_ref, "Mean of signal in reference band")
	prihdr.set('DT', dt, "seconds")
	prihdr.set('N_BINS', n_bins, "time bins per segment")
	prihdr.set('SECONDS', num_seconds, "seconds per segment")
	prihdr.set('SEGMENTS', num_segments, "segments in the whole light curve")
	prihdr.set('EXPOSURE', num_segments * n_bins * dt, \
		"seconds, of light curve")
	prihdr.set('RATE_CI', str(mean_rate_whole_ci.tolist()), "counts/second")
	prihdr.set('RATE_REF', mean_rate_whole_ref, "counts/second")
	prihdr.set('NOISY', noisy, "True if adding Poisson noise to simulation")
	prihdu = fits.PrimaryHDU(header=prihdr)
	
	chan = np.arange(0,64)
	energy_channels = np.tile(chan, len(t))
	ccf_error = np.tile(ccf_error, len(t))
	time_bins = np.repeat(t, len(chan))
	assert len(energy_channels) == len(time_bins)
	## Making FITS table
	col1 = fits.Column(name='TIME_BIN', format='K', array=time_bins)
	col2 = fits.Column(name='CCF', unit='Counts/second', format='D', \
		array=ccf_filtered.real.flatten('C'))
	col3 = fits.Column(name='ERROR', unit='', format='D', array=ccf_error)
	col4 = fits.Column(name='CHANNEL', unit='', format='I', \
		array=energy_channels)
	cols = fits.ColDefs([col1, col2, col3, col4])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	## If the file already exists, remove it (still working on just updating it)
	assert out_file[-4:].lower() == "fits", 'ERROR: Output file must have extension ".fits".'
	if os.path.isfile(out_file):
		os.remove(out_file)
	## Writing to a FITS file
	thdulist = fits.HDUList([prihdu, tbhdu])
	thdulist.writeto(out_file)	
## End of function 'fits_output'

	
###############################################################################
def dat_output(out_file, bb_spec, pl_spec, freq, phase, dt, n_bins, \
	num_seconds, num_segments, amp_ci, amp_ref, mean_ci, mean_ref, \
	mean_rate_whole_ci, mean_rate_whole_ref, t, ccf_filtered, ccf_error, noisy):
    """
            dat_output

    Writes the simulation parameters and cross-correlation function to an output
    file.

    Passed: out_file - Name of output file.
    		bb_spec - 
    		pl_spec - 
    		freq - 
    		phase - 
            dt - Size of each time bin, in seconds.
            n_bins - Number of (time) bins per segment.
            num_seconds - Number of seconds in each Fourier segment.
            num_segments - Number of segments the light curve was split up into.
            amp_ci - 
            amp_ref - 
            mean_ci - 
            mean_ref - 
            mean_rate_whole_ci - Mean count rate of light curve 1, averaged over
                all segments.
            mean_rate_whole_ref - Mean count rate of light curve 2, averaged
                over all segments.
            t - Integer time bins to plot against the ccf.
            ccf_filtered - CCF amplitudes, filtered in frequency space.
            ccf_error - Error on the filtered CCF.

    Returns: nothing

    """
    print "Output sent to %s" % out_file

    with open(out_file, 'w') as out:
        out.write("#\t\tCross-correlation function of simulated data")
        out.write("\n# Date(YYYY-MM-DD localtime): %s" % str(datetime.now()))
        out.write("\n# Blackbody (bb) spectrum: %s" % bb_spec)
        out.write("\n# Power law (pl) spectrum: %s" % pl_spec)
        out.write("\n# Frequency of signal = %.2f Hz" % freq)
        out.write("\n# Phase shift of PL with respect to BB = %f" % phase)
        out.write("\n# Time bin size = %.21f seconds" % dt)
        out.write("\n# Number of bins per segment = %d" % n_bins)
        out.write("\n# Number of seconds per segment = %d" % num_seconds)
        out.write("\n# Number of segments per light curve = %d" % num_segments)
        out.write("\n# Exposure time = %d seconds" \
                  % (num_segments * num_seconds))
        out.write("\n# Fractional amplitude of signal in channels of interest = %f" % amp_ci)
        out.write("\n# Fractional amplitude of signal in reference band = %f" % amp_ref)
        out.write("\n# Fractional mean of signal in channels of interest = %f" % mean_ci)
        out.write("\n# Fractional mean of signal in reference band = %f" % mean_ref)
        out.write("\n# Mean count rate of ci = %s" \
        	% str(list(mean_rate_whole_ci)))
        out.write("\n# Mean count rate of ref band = %.8f" \
                  % np.mean(mean_rate_whole_ref))
        out.write("\n# Noisy = %s" % str(noisy))
        out.write("\n# ")
        out.write("\n# Column 1: Time bins")
        out.write("\n# Columns 2-65: Filtered ccf per energy channel, real \
            part [count rate]")
        out.write("\n# Columns 66-129: Error on filtered ccf per energy \
            channel, real part [count rate]")
        out.write("\n# ")
        for j in xrange(0, n_bins):
            out.write("\n%d" % t[j])
            for i in xrange(0, 64):
                out.write("\t%.6e" % ccf_filtered[j][i].real)
            for i in xrange(0, 64):
                out.write("\t%.6e" % ccf_error[i].real)
        ## End of for-loops
    ## End of with-block
    
## End of function 'dat_output'


###############################################################################
def main(out_file, bb_spec, pl_spec, freq, dt_mult, num_seconds, amp_ci, \
	amp_ref, mean_ci, mean_ref, phase, test, noisy):
	"""
			main
			
	Summary.
	
	Passed: 
	
	Returns:
	
	"""
	t_res = 1.0 / 8192.0  # The time resolution of the data, in seconds.
	dt = dt_mult * t_res
	n_bins = num_seconds * int(1.0 / dt)
	assert tools.power_of_two(n_bins)  # n_bins must be a power of 2 for the FFT
	
	exposure_bb = float(tools.get_key_val(bb_spec, 1, "EXPOSURE"))
	exposure_pl = float(tools.get_key_val(pl_spec, 1, "EXPOSURE"))
	assert exposure_bb == exposure_pl
	exposure = exposure_bb
	
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
	sine_ci_bb, sine_ref_bb = sim_lc.generate_sines(dt, n_bins, freq, amp_ci, \
		amp_ref, mean_ci, mean_ref, 0.0)
	sine_ci_pl, sine_ref_pl = sim_lc.generate_sines(dt, n_bins, freq, amp_ci, \
		amp_ref, mean_ci, mean_ref, phase)
	curve_ci_bb, curve_ref_bb = sim_lc.make_lightcurves(spec_bb, sine_ci_bb, \
		sine_ref_bb, n_bins)
	curve_ci_pl, curve_ref_pl = sim_lc.make_lightcurves(spec_pl, sine_ci_pl, \
		sine_ref_pl, n_bins)
	
		
	############################
	## Looping through segments
	############################
# 	for num_segments in xrange(1, 41948): # tracks the number of segments
# 	for num_segments in xrange(1, 101):  # tracks the number of segments
	for num_segments in xrange(1, 15000): # tracks the number of segments


		curve_ci, curve_ref = sim_lc.add_lightcurves(curve_ci_bb, curve_ref_bb,\
			curve_ci_pl, curve_ref_pl, dt, exposure, noisy)
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

		cs_segment, mean_rate_segment_ci, mean_rate_segment_ref, power_ci, \
			power_ref = crosscorr.make_cs(curve_ci, curve_ref, n_bins, dt)
		
		sum_rate_whole_ci += mean_rate_segment_ci
		sum_rate_whole_ref += mean_rate_segment_ref
		sum_power_ci += power_ci
		sum_power_ref += power_ref
		cs_sum += cs_segment  # This adds indices
		
# 		sum_rate_ci += np.mean(curve.ci)
		if num_segments % 10 == 0: print "\t", num_segments
		if num_segments == 1 and args.test == True: 
			break
		
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
		mean_power_ci, mean_power_ref, noisy)

	ccf_filtered[np.where(np.isnan(ccf_filtered))] = 0.0
	ccf_error[np.where(np.isnan(ccf_error))] = 0.0
	
	t = np.arange(0, n_bins)
	if out_file[-3:].lower() == "dat":
		dat_output(out_file, bb_spec, pl_spec, freq, phase, dt, n_bins, \
			num_seconds, num_segments, amp_ci, amp_ref, mean_ci, mean_ref, \
			mean_rate_whole_ci, mean_rate_whole_ref, t, ccf_filtered, \
			ccf_error, noisy)
	elif out_file[-4:].lower() == "fits":
		fits_output(out_file, bb_spec, pl_spec, freq, phase, dt, n_bins, \
			num_seconds, num_segments, amp_ci, amp_ref, mean_ci, mean_ref, \
			mean_rate_whole_ci, mean_rate_whole_ref, t, ccf_filtered, \
			ccf_error, noisy)
	else:
		raise Exception("ERROR: Output file must have extension .dat or .fits.")
# 	sim_lc.plot_curves(n_bins, mean_curve_ci[:,6], mean_curve_ref, "plot.png")

	sim_lc.power_spectra_things(mean_ps_ref, dt, n_bins, num_seconds, \
		num_segments, mean_rate_ref, noisy)
	

	
###############################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Simulates a light curve and \
		runs it through ccf.py.')
	parser.add_argument('out_file', help='Output file name.')
	parser.add_argument('--freq', type=float, required=True, dest='freq', \
		help='Frequency of the periodic pulsation, in Hz.')
	parser.add_argument('--bb', required=True, dest='bb_spec', help='Name of \
		.fak spectrum file for the blackbody component of the energy spectrum.')
	parser.add_argument('--pl', required=True, dest='pl_spec', help='Name of \
		.fak spectrum file for the power law component of the energy spectrum.')
	parser.add_argument('--num_seconds', type=tools.type_power_of_two, \
		default=1, dest='num_seconds', help='Number of seconds in each segment.\
		Must be a power of 2. [1]')
	parser.add_argument('--dt_mult', type=tools.type_power_of_two, default=1, \
		dest='dt_mult', help='Multiple of 1/8192 seconds for timestep between \
		bins. Must be a power of 2. [1]')
	parser.add_argument('--mean_ci', type=tools.type_positive_float, \
		default=1.0, dest='mean_ci', help='Mean value of the signal for the \
		channels of interest. [1.0]')
	parser.add_argument('--mean_ref', type=tools.type_positive_float, \
		default=1.0, dest='mean_ref', help='Mean value of the signal for the \
		reference band. [1.0]')
	parser.add_argument('--amp_ci', type=tools.type_positive_float, \
		default=0.2, dest='amp_ci', help='Fractional amplitude of the signal \
		for the channels of interest. [0.2]')
	parser.add_argument('--amp_ref', type=tools.type_positive_float, \
		default=0.2, dest='amp_ref', help='Fractional amplitude of the signal \
		for the reference band. [0.2]')
	parser.add_argument('--phase', type=float, default=0.0, \
		dest='phase', help='Phase difference of the power law variability \
		to the blackbody variability in the energy spectrum, in degrees. [0.0]')
	parser.add_argument('--test', action='store_true', dest='test', help='If \
		present, only does a short test run.')
	parser.add_argument('--noisy', action='store_true', dest='noisy', help='If \
		present, adds Poisson noise to the fake data.')
	args = parser.parse_args()
	
	main(args.out_file, args.bb_spec, args.pl_spec, args.freq, args.dt_mult, \
		args.num_seconds, args.amp_ci, args.amp_ref, args.mean_ci, \
		args.mean_ref, np.deg2rad(args.phase), args.test, args.noisy)
	
