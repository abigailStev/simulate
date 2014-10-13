import argparse
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from astropy.io import fits

"""
		simulate_lightcurves

Creates two synthetic light curves by adding Poisson noise to a sine wave.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as Python 2.7, can be downloaded 
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

"""
###############################################################################
def plot_curves(n_bins, curve_1, curve_2, plot_file):
	"""
			plot_curves
	
	Plots two light curves of length n_bins.
	
	Passed: n_bins - Length of light curves. Must be a power of 2.
			curve_1 - Amplitudes of light curve 1
			curve_2 - Amplitudes of light curve 2
			plot_file - Name of plot file
			
	Returns: nothing
	
	"""
	
	bins = np.arange(0,n_bins) # Bins to plot against

	fig, ax = plt.subplots()
	ax.plot(bins, curve_1, linewidth=1.5, label="Curve 1")
	ax.plot(bins, curve_2, linewidth=1.5, label="Curve 2")
	plt.xlim(0, n_bins)
	plt.xlabel('Bins')
	plt.ylabel('Amplitude')

	## The following legend code was found on stack overflow I think, 
	## or a pyplot tutorial
	legend = ax.legend(loc='upper right')

	## Set the fontsize
	for label in legend.get_texts():
		label.set_fontsize('small')

	for label in legend.get_lines():
		label.set_linewidth(2)  # the legend line width

	plt.savefig(plot_file)

	plt.show()
	
	## End function 'plot_curves'
	
	
###############################################################################
def generate_sines(dt, n_bins, freq, amp_ci, amp_ref, mean_ci, mean_ref, \
	phase_ci, noisy):
	"""
			generate_sines
		
	Creates two synthetic light curves with Poisson noise on sine wave signals 
	of length n_bins. The ci curve is shifted by phase angle 'phase_ci'.
		
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
			mean_ci - Mean of the sine wave signal for ci, in count rate units.
			mean_ref - Mean of the sine wave signal for ref, in count rate 
				units.
			phase_ci - The phase shift for the ci sine wave signal, in radians?
			noisy - True for Poisson noise, False for no noise
	
	Returns: noisy_sine_ci - Noisified ci of length n_bins
			 noisy_sine_ref - Noisified ref of length n_bins
			 
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
	
	
	if noisy == False: ## For the 'no noise at all' case
		bins = np.arange(0,n_bins, 1)
		smooth_sine = [ (amp_ci * np.sin(2.0 * np.pi * x / bins_per_period + phase_ci) + mean_ci) for x in bins]
		another_smooth_sine = [ (amp_ref * np.sin(2.0 * np.pi * x / bins_per_period) + mean_ref) for x in bins]
		
		return smooth_sine, another_smooth_sine  ## for the 'no noise at all' case
	## End 'if not noisy'
	else: ## For the regular case with Poisson noise
		## Making two sine waves that are sampled over tiny_bins, so they're 
		## very smooth.
		tiny_bins = np.arange(0, n_bins, 0.1)
		smooth_sine_ci = [ (amp_ci * np.sin(2.0 * np.pi * x / bins_per_period + phase_ci) + mean_ci) for x in tiny_bins] # in units 'rate'
		smooth_sine_ref = [ (amp_ref * np.sin(2.0 * np.pi * x / bins_per_period) + mean_ref) for x in tiny_bins] # in units 'rate'

		## Taking the average amplitude of every 10 bins of smooth_sine as the 
		## value for sine
		i = 0
		j = 10
		while j <= len(tiny_bins):
			sine_ci.append(np.mean(smooth_sine_ci[i:j]))
			sine_ref.append(np.mean(smooth_sine_ref[i:j]))
			i = j
			j += 10
			## End of while loop
	
		## Adding Poisson noise to sine_ci and sine_ref
		## Dividing by dt puts values into 'count rate' units
		noisy_sine_ci = np.random.poisson(np.asarray(sine_ci) * dt) / dt
		noisy_sine_ref = np.random.poisson(np.asarray(sine_ref) * dt) / dt

# 		print "Mean noisy sine ci =", np.mean(noisy_sine_ci), ", mean noisy sine ref =", np.mean(noisy_sine_ref)
		
		return noisy_sine_ci, noisy_sine_ref
	## End 'else noisy'
		
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
	
	if spec_file[-3:] == "dat":
		table = np.loadtxt(spec_file, dtype=float)
		spectrum = table[:,1]	
	
	elif spec_file[-3:] == "fak":
		file_hdu = fits.open(spec_file)
		spectrum = file_hdu[1].data.field(0)
	else:
		print "\n\tERROR: Spectrum format not recognized. Exiting."
		exit()
	
	assert len(spectrum) == 64
# 	print np.shape(spectrum)

	return spectrum
## End of function 'read_fakeit_spectra'
	
	
###############################################################################
def make_lightcurves(spec_ci, spec_ref, sine_ci, sine_ref):
	"""
			make_lightcurves
	
	Multiplies a spectrum by a sinusoid to create a fake lightcurve per energy 
	channel, for ci and reference.
	
	Passed: spec_ci - 
			spec_ref - 
			sine_ci - 
			sine_ref - 
	
	Returns: curve_ci - 
			 curve_ref - 
	
	"""
	
	print "Shape of ci spectrum:", np.shape(spec_ci)
	print "Shape of sine ci:", np.shape(np.reshape(sine_ci, (n_bins,1)))
	curve_ci = np.multiply(np.reshape(sine_ci, (n_bins,1)), \
		spec_ci[np.newaxis])
	print "Shape of curve_ci:", np.shape(curve_ci)
	
	print "Shape of ref spectrum:", np.shape(spec_ref)
	print "Shape of sine ref:", np.shape(sine_ref)
	curve_ref = np.multiply(np.reshape(sine_ref, (n_bins,1)), \
		spec_ref[np.newaxis])
	print "Shape of curve_ref:", np.shape(curve_ref)
	
	return curve_ci, curve_ref
	
## End of function 'make_lightcurves'


###############################################################################
if __name__ == "__main__":
	"""
	"""
	dt = 1.0 / 8192.0  # The timestep or amount of time per bin, in seconds.
# 	n_bins = 128  # Number of bins in one light curve (assuming both curves are 
				  # of same length).
	n_bins = 8192
	freq = 401.0  # Frequency of sine wave signals, in Hz. Works well when this 	
				  # is a power of 2 (otherwise we get aliasing). Assuming that 
				  # the signals of both light curves have the same frequency.
	mean_ci = 5.5  # Mean count rate of ci.
	amp_ci = 1.0  # Amplitude of sine wave signal for ci
	mean_ref = 5.5  # Mean count rate of ref; when summed for the 24 energy 
	                 # channels used in reference band, should get 130 cts/sec
	amp_ref = 1.0  # Amplitude of sine wave signal for ref
	phase_ci = 0.0  # The phase shift of the ci sine wave signal
# 	amp_ci = 0.0	 # Amplitude of sine wave signal for ci; No signal
# 	amp_ref = 0.0  # Amplitude of sine wave signal for ref; No signal
	noisy = True  # Boolean flag: True gives a noisy light curve, False gives a 
				  # smooth one.
	
	exposure_time = 2 # ks  100  # ks
	
# 	parser = argparse.ArgumentParser()
# 	parser.add_argument('plot_file', help="Name of output plot file.")
# 	args = parser.parse_args()
	
# 	spectrum_file = "/Users/abigailstevens/Dropbox/Research/energy_spectra/out_es/P70080_140917_t1_4sec_pbin_35.dat"
	spectrum_file = "fakeit_mean.fak"
	
	spec_ci = read_fakeit_spectra(spectrum_file)
	spec_ref = read_fakeit_spectra(spectrum_file)
	
	sine_ci, sine_ref = generate_sines(dt, n_bins, freq, amp_ci, amp_ref, 
		mean_ci, mean_ref, phase_ci, noisy)
	
	curve_ci, curve_ref = make_lightcurves(spec_ci, spec_ref, sine_ci, sine_ref)
	
# 	plot_curves(n_bins, sine_ci, sine_ref, "plot.png")

## End of main function

## End of 'generate_sines.py'