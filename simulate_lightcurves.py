import argparse
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from astropy.io import fits
import tools
from ccf import stack_reference_band
import powerspec
from scipy import fftpack

"""
		simulate_lightcurves

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as Python 2.7, can be downloaded 
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

'tools' can be found in the github repo 'whizzy_scripts'
'stack_reference_band' can be found in the github repo 'cross_correlation'

"""
###############################################################################
def plot_curves(n_bins, curve_ci, curve_ref, plot_file):
	"""
			plot_curves
	
	Plots two light curves of length n_bins.
	
	Passed: n_bins - Length of light curves. Must be a power of 2.
			curve_1 - Amplitudes of ci light curve.
			curve_2 - Amplitudes of ref light curve.
			plot_file - Name of file to save plot to.
			
	Returns: nothing
	
	"""
	
	bins = np.arange(0,n_bins) # Bins to plot against

	fig, ax = plt.subplots()
	ax.plot(bins, curve_ci, linewidth=1.5, label="Curve 'ci'")
	ax.plot(bins, curve_ref, linewidth=1.5, label="Curve 'ref'")
	plt.xlim(0, 30)
	plt.xlabel('Arbitrary time bins')
	plt.ylabel('Count rate')

	## The following legend code was found on stack overflow I think, 
	## or a pyplot tutorial
	legend = ax.legend(loc='upper right')

	## Set the fontsize
	for label in legend.get_texts():
		label.set_fontsize('small')

	for label in legend.get_lines():
		label.set_linewidth(2)  # the legend line width

	plt.savefig(plot_file)

# 	plt.show()
	
## End function 'plot_curves'
	
	
###############################################################################
def power_spectra_things(mean_ps_ref, dt, n_bins, num_seconds, num_segments, \
	mean_rate_ref, noisy):
	"""
			power_spectra_things
	
	Computes power spectrum things and plots the power spectrum.
	
	"""
	ps_freq, mean_ps_ref, leahy_power, rms2_power, rms2_err_power = \
		powerspec.normalize(mean_ps_ref, n_bins, dt, num_seconds, num_segments,\
		mean_rate_ref, noisy)
	
	np.savetxt( "sim_power.dat", mean_ps_ref)
	fig, ax = plt.subplots()
	ax.plot(ps_freq, rms2_power, linewidth=2)
# 	ax.plot(ps_freq, leahy_power-2, linewidth=2)
	plt.xlabel(r'$\nu$ [Hz]')
	plt.ylabel(r'Noise-subtracted fractional rms$^2$ power')
	plt.xlim(0,800)
	plt.ylim(0,)
	# plt.xscale('symlog') # this works much better than 'log'
	# plt.yscale('symlog')
	plt.title("Power spectrum")
	plt.savefig("sim_power.png", dpi=120)
# 	plt.show()
	plt.close()
## End of function 'power_spectra_things'
	
	
###############################################################################
def determine_extra_bins(bpp):
	"""
	
	"""
	
	i = 1
	while np.abs(bpp*i - np.rint(bpp*i)) > 0.1 and \
		np.abs(bpp*i - np.rint(bpp*i)) < 0.9 and i < 100:
		i+= 1
	
	print np.rint(bpp*i)
	return np.rint(bpp*i)
## End of function 'determine_extra_bins'
	

###############################################################################
def generate_sines(dt, n_bins, freq, amp_ci, amp_ref, mean_ci, mean_ref, \
	phase):
	"""
			generate_sines
		
	Creates two synthetic light curves with Poisson noise on sine wave signals 
	of length n_bins. For the blackbody spectrum, 'phase_diff' is zero. For the 
	power law spectrum, 'phase_diff' is non-zero.
		
	Passed: dt - The timestep or amount of time per bin, in seconds.
			n_bins - Number of bins, the integer length of the light curve. 
				Must be a power of two.
			freq - Frequency of the sine wave signals, in Hz. Works well when 
				this is a power of 2 (otherwise we get aliasing). We assume that 
				the signals of both light curves have the same frequency.
			amp_ci - Amplitude of the sine wave signal for ci. To simulate a 
				noise-only process, set amp_ci = 0.
			amp_ref - Amplitude of the sine wave signal for ref. To simulate a 
				noise-only process, set amp_ref = 0.
			mean_ci - Mean of the sine wave signal for ci.
			mean_ref - Mean of the sine wave signal for ref.
			phase - The phase shift for the power law variability, in 
				radians?
	
	Returns: sine_ci - Relative sine wave signal for ci of length n_bins.
			 sine_ref - Relative sine wave signal for ref of length n_bins.
			 
	"""
	
	## Initializing new variables
	period = 1.0/freq  # Period of sine waves, in seconds
	bins_per_period = period / dt  # Number of bins per period of sine wave 
								   # signal
	smooth_sine_ci = []  # A smooth, frequently sampled sine wave (with amp_ci 
						 # and mean_ci)
	smooth_sine_ref = []  # A smooth, frequently sampled sine wave (with amp_ref 
						  # and mean_ref)
	sine_ci = []  # Averaged smooth_sine_ci over every 10 tiny_bins; this is 
				  # what gets Poisson noise
	sine_ref = []  # Averaged smooth_sine_ref over every 10 tiny_bins; this is 
				   # what gets Poisson noise
	
	extra_bins = determine_extra_bins(bins_per_period)
	
	## Making two sine waves that are sampled over tiny_bins, so they're very 
	## smooth.
	tiny_bins = np.arange(0, n_bins+extra_bins, 0.1)

	smooth_sine_ci = amp_ci * np.sin(2.0 * np.pi * tiny_bins / bins_per_period + phase) + mean_ci
	smooth_sine_ref = amp_ref * np.sin(2.0 * np.pi * tiny_bins / bins_per_period + phase) + mean_ref
	
	
	
	
	
	sine_ci = np.mean(np.array_split(smooth_sine_ci, n_bins+extra_bins), axis=1)
	sine_ref = np.mean(np.array_split(smooth_sine_ref, n_bins+extra_bins), axis=1)
	
# 	print np.shape(sine_ci)
# 	print np.shape(sine_ref)
	
	return sine_ci, sine_ref, extra_bins
		
## End of function 'generate_sines'


###############################################################################
def read_fakeit_spectra(spec_file):
	"""
			read_fakeit_spectra
	
	Reads spectra (created by the heasoft tool 'fakeit') from a text file.
	
	Passed: spec_file - Filename of energy spectrum.
	
	Returns: spectrum - 
	
	"""
	
# 	print spec_file[-3:]
	print "Energy spectrum file: %s" % spec_file
	if spec_file[-3:] == "dat":
		table = np.loadtxt(spec_file, dtype=float)
		spectrum = table[:,1]	
	
	elif spec_file[-3:] == "fak":
		file_hdu = fits.open(spec_file)
		spectrum = file_hdu[1].data.field('COUNTS')
	else:
		raise Exception("ERROR: Spectrum format not recognized. Must be a fits file of type '.fak' or a text file of type '.dat'.")
		
	assert len(spectrum) == 64
# 	print np.shape(spectrum)
# 	print spectrum

	return spectrum
## End of function 'read_fakeit_spectra'


###############################################################################
def make_lightcurves(spec_xx, sine_ci, sine_ref, n_bins):
	"""
			make_lightcurves
	
	Multiplies a spectrum by a sinusoid to create a fake light curve per energy 
	channel, for ci and reference.
	
	Passed: spec_xx - 
			sine_ci - 
			sine_ref - 
			n_bins - 
	
	Returns: curve_ci_xx - 
			 curve_ref_xx - 
	
	"""
	curve_ci_xx = np.multiply(np.reshape(sine_ci, (n_bins,1)), \
		spec_xx[np.newaxis])
	curve_ref_xx = np.multiply(np.reshape(sine_ref, (n_bins,1)), \
		spec_xx[np.newaxis])
	curve_ref_xx = stack_reference_band(curve_ref_xx, 5)

	## Multiplying by a normalization factor so that we get ~200 counts/sec in 
	## the reference band
	
	curve_ci_xx *= 15.0
	curve_ref_xx *= 15.0
	
	return curve_ci_xx, curve_ref_xx
## End of function 'make_lightcurves'


###############################################################################
def add_lightcurves(curve_ci_bb, curve_ref_bb, curve_ci_pl, curve_ref_pl, dt, \
	exposure, noisy):
	"""
			add_lightcurves
	
	Adds together the light curves from the blackbody and the power law, and 
	adds Poisson noise to the light curve.
	
	Passed: curve_ci_bb - 
			curve_ref_bb - 
			curve_ci_pl - 
			curve_ref_pl - 
			dt - 
			exposure -
			noisy - 
	
	Returns: noisy_curve_ci - 
			 noisy_curve_ref - 
	
	"""
	
	curve_ci = curve_ci_bb + curve_ci_pl
	curve_ref = curve_ref_bb + curve_ref_pl
	
	if noisy:
		## Adding Poisson noise to curve_ci and curve_ref, changing to count rate units
		noisy_curve_ci = np.random.poisson(curve_ci * dt / exposure) / dt
		noisy_curve_ref = np.random.poisson(curve_ref * dt / exposure) / dt
		return noisy_curve_ci, noisy_curve_ref
	else:
		## Not contributing poisson noise, but still changing to count rate units
		rate_curve_ci = curve_ci / exposure
		rate_curve_ref = curve_ref / exposure
		return rate_curve_ci, rate_curve_ref

## End of function 'add_lightcurves'


###############################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Simulates the light curve of \
		a periodic pulsation with the blackbody component varying out of phase \
		with the power law component of the energy spectrum.', epilog='For \
		optional arguments, default values are given in brackets at end of \
		description.')
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
		default=1.0, dest='mean_ci', \
		help='Mean value of the signal for the channels of interest. [1.0]')
	parser.add_argument('--mean_ref', type=tools.type_positive_float, \
		default=1.0, dest='mean_ref', help='Mean value of the signal for the \
		reference band. [1.0]')
	parser.add_argument('--amp_ci', type=tools.type_positive_float, \
		default=0.5, dest='amp_ci', help='Fractional amplitude of the signal \
		for the channels of interest. [0.5]')
	parser.add_argument('--amp_ref', type=tools.type_positive_float, \
		default=0.5, dest='amp_ref', help='Fractional amplitude of the signal \
		for the reference band. [0.5]')
	parser.add_argument('--phase', type=float, default=0.0, \
		dest='phase', help='Phase difference of the power law variability \
		to the blackbody variability in the energy spectrum. [0.0]')
	parser.add_argument('--exposure', type=tools.type_positive_float, \
		default=1000.0, dest='exposure', help='Exposure time of the \
		observation, in seconds. [1000.0]')
	parser.add_argument('--test', action='store_true', dest='test', help='If \
		present, only does a short test run.')
	parser.add_argument('--noisy', action='store_true', dest='noisy', help='If \
		present, adds Poisson noise to the fake data.')
	args = parser.parse_args()


	
	t_res = 1.0 / 8192.0  # The time resolution of the data, in seconds.
	dt = args.dt_mult * t_res
	n_bins = args.num_seconds * int(1.0 / dt)
	assert tools.power_of_two(n_bins)  # n_bins must be a power of 2 for the FFT
	
	spec_bb = read_fakeit_spectra(args.bb_spec)
	spec_pl = read_fakeit_spectra(args.pl_spec)
	mean_curve_ci = 0
	mean_curve_ref = 0
	mean_ps_ref = 0
	mean_rate_ref = 0
	
	sine_ci, sine_ref, extra_bins = generate_sines(dt, n_bins, args.freq, \
		args.amp_ci, args.amp_ref, args.mean_ci, args.mean_ref, args.phase)
	curve_ci_bb, curve_ref_bb = make_lightcurves(spec_bb, sine_ci, sine_ref, \
		n_bins+extra_bins)
	curve_ci_pl, curve_ref_pl = make_lightcurves(spec_pl, sine_ci, sine_ref, \
		n_bins+extra_bins)
	
	extra_bin_choices = np.arange(int(extra_bins))
	
	for i in xrange(1, 101):  # 'i' tracks the number of segments
		start_bin = np.random.choice(extra_bin_choices)
		
		curve_ci, curve_ref = add_lightcurves(curve_ci_bb[start_bin:n_bins+start_bin], \
			curve_ref_bb[start_bin:n_bins+start_bin], \
			curve_ci_pl[start_bin:n_bins+start_bin], \
			curve_ref_pl[start_bin:n_bins+start_bin], dt, args.exposure, \
			args.noisy)
		mean_curve_ci += curve_ci
		mean_curve_ref += curve_ref
		ps_ref, rate_ref = powerspec.make_ps(curve_ref)
		mean_ps_ref += ps_ref
		mean_rate_ref += rate_ref
		
		if i == 1 and args.test == True: 
			break
		
		if i % 10 == 0: print "\t", i
	## End of for-loop
	
	mean_ps_ref /= i
	mean_rate_ref /= i
	print "Mean count rate in reference band:", mean_rate_ref
	
	plot_curves(n_bins, mean_curve_ci[:,6], mean_curve_ref, "plot.png")

	power_spectra_things(mean_ps_ref, dt, n_bins, args.num_seconds, i, \
		mean_rate_ref, args.noisy)
	
## End of main function

## End of 'simulate_lightcurves.py'