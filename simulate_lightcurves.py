import argparse
import math
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

"""
		simulate_lightcurves

Creates two synthetic light curves by adding Poisson noise to a sine wave.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as Python 2.7, can be downloaded 
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

"""

###############################################################################
def generate_sines(dt, n_bins, freq, amp_ci, amp_ref, mean_ci, mean_ref, noisy):
	"""
			generate_sines
		
	Creates two synthetic light curves with Poisson noise on sine wave signals 
	of length n_bins.
		
	Passed: dt - float - The timestep or amount of time per bin, in seconds.
			n_bins - int - Number of bins, the integer length of the light 
				curve. Must be a power of two.
			freq - float - Frequency of the sine wave signals, in Hz. Works well
				when this is a power of 2 (otherwise we get aliasing). We assume 
				that the signals of both light curves have the same frequency.
			amp_ci - float - Amplitude of the sine wave signal for ci.
				To simulate a noise-only process, set amp_ci = 0.
			amp_ref - float - Amplitude of the sine wave signal for ref.
				To simulate a noise-only process, set amp_ref = 0.
			mean_ci - float - Mean of the sine wave signal for ci, in 
				count rate units.
			mean_ref - float - Mean of the sine wave signal for ref, in 
				count rate units.
			noisy - boolean - True for Poisson noise, False for no noise
	
	Returns: noisy_sine_ci - list of floats - Noisified ci of length n_bins
			 noisy_sine_ref - list of floats - Noisified ref of length n_bins
			 
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

#  	print "Bins per period:", bins_per_period
	
	if noisy == False: ## For the 'no noise at all' case
		
		bins = np.arange(0,n_bins, 1)
		smooth_sine = [ (amp_ci * math.sin(2.0 * math.pi * x / bins_per_period) + mean_ci) for x in bins]
		another_smooth_sine = [ (amp_ref * math.sin(2.0 * math.pi * x / bins_per_period) + mean_ref) for x in bins]
# 		print smooth_sine
		
		return smooth_sine, another_smooth_sine  ## for the 'no noise at all' case
	## End 'if not noisy'
	else: ## For the regular case with Poisson noise
		## Making two sine waves that are sampled over tiny_bins, so they're 
		## very smooth.
		tiny_bins = np.arange(0, n_bins, 0.1)
		smooth_sine_ci = [ (amp_ci * math.sin(2.0 * math.pi * x / bins_per_period) + mean_ci) for x in tiny_bins] # in units 'rate'
		smooth_sine_ref = [ (amp_ref * math.sin(2.0 * math.pi * x / bins_per_period) + mean_ref) for x in tiny_bins] # in units 'rate'

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
		noisy_sine_ci = [np.random.poisson(x*dt) / dt for x in sine_ci] 
		# puts values back into 'rate' units
		noisy_sine_ref = [np.random.poisson(x*dt) / dt for x in sine_ref] 
		# puts values back into 'rate' units
		
# 		print "Mean noisy sine ci =", np.mean(noisy_sine_ci), ", mean noisy sine ref =", np.mean(noisy_sine_ref)

		return noisy_sine_ci, noisy_sine_ref
	## End 'else noisy'
		
## End of function 'generate_sines'


###############################################################################
def read_fakeit_spectra(spec_file)
	"""
			read_fakeit_spectra
	
	Reads spectra (created by the heasoft tool 'fakeit') from a text file.
	
	Passed: spec_file - Filename of energy spectrum.
	
	Returns: spectrum - 
	
	"""
	
	spectrum = np.loadtxt(spec_file, dtype=float)
	
	print len(spectrum)
	assert len(spectrum) == 64
	
	return spectrum
	
	

###############################################################################
if __name__ == "__main__":
	"""
	"""
	dt = 1.0 / 8192.0  # The timestep or amount of time per bin, in seconds.
# 	n_bins = 128  # Number of bins in one light curve (assuming both curves are 
				  # of same length).
	n_bins = 32768
	freq = 401.0  # Frequency of sine wave signals, in Hz. Works well when this 	
				  # is a power of 2 (otherwise we get aliasing). Assuming that 
				  # the signals of both light curves have the same frequency.
	mean_ci = 10  # Mean count rate of ci.
	amp_ci = 1.0  # Amplitude of sine wave signal for ci
	mean_ref = 200  # Mean count rate of ref.
	amp_ref = 14.0  # Amplitude of sine wave signal for ref
# 	amp_ci = 0.0	 # Amplitude of sine wave signal for ci; No signal
# 	amp_ref = 0.0  # Amplitude of sine wave signal for ref; No signal
	noisy = True  # Boolean flag: True gives a noisy light curve, False gives a 
				  # smooth one.
	
	exposure_time = 2 # ks  100  # ks
	
# 	parser = argparse.ArgumentParser()
# 	parser.add_argument('plot_file', help="Name of output plot file.")
# 	args = parser.parse_args()
	
	curve_ci, curve_ref = generate_sines(dt, n_bins, freq, amp_ci, amp_ref, 
		mean_ci, mean_ref, noisy)
# 	print curve_ci

## End of 'generate_sines.py'