import argparse
from datetime import datetime
import os
import numpy as np
from scipy import fftpack
import matplotlib.pyplot as plt
import itertools
from astropy.io import fits

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens at uva.nl"
__year__ = "2015"
__description__ = "Simulates a power spectrum using the Timmer and Koenig \
method."

"""
		TimmerKoenig.py

Written in Python 2.7.

"""

################################################################################
def lc_out(dt, time, lightcurve):
	out_data = np.column_stack((time, lightcurve))
	np.savetxt('./TK_lightcurve.dat', out_data)
	
	detchans=64
	
	out_file = "./TK_lightcurve.lc"
	## Making header for standard power spectrum
	prihdr = fits.Header()
	prihdr.set('TYPE', "Fake light curve generated from T&K")
	prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
	prihdr.set('TIMEDEL', 0.0001220703125)
	prihdr.set('DETCHANS', detchans, "Number of detector energy channels")
	prihdu = fits.PrimaryHDU(header=prihdr)
	
	## Making FITS table for standard power spectrum
	col1 = fits.Column(name='TIME', unit='s', format='1D', array=time)
	col2 = fits.Column(name='RATE', unit='photon/s', format='1D', \
		array=lightcurve)
	cols = fits.ColDefs([col1, col2])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	
	## If the file already exists, remove it (still working on just updating it)
# 	assert out_file[-4:].lower() == "fits", \
# 		'ERROR: Standard output file must have extension ".fits".'
	if os.path.isfile(out_file):
# 		print "File previously existed. Removing and rewriting."
		os.remove(out_file)
		
	## Writing the standard power spectrum to a FITS file
	thdulist = fits.HDUList([prihdu, tbhdu])
	thdulist.writeto(out_file)	
## End of function 'lc_out'


################################################################################
def power_out(freq, power, dt, n_bins, nyquist, num_seg, mean_rate):
	out_file = "./TK_power.fits"
	detchans=64

	## Making header for standard power spectrum
	prihdr = fits.Header()
	prihdr.set('TYPE', "Fake power spectrum generated from T&K")
	prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
	prihdr.set('DETCHANS', detchans, "Number of detector energy channels")
	prihdr.set('DT', dt)
	prihdr.set('N_BINS', n_bins)
	prihdr.set('NYQUIST', nyquist)
	prihdr.set('SEGMENTS', num_seg)
	prihdr.set('MEANRATE', mean_rate)
	prihdu = fits.PrimaryHDU(header=prihdr)
	
	## Making FITS table for standard power spectrum
	col1 = fits.Column(name='FREQUENCY', unit='Hz', format='E', array=freq)
	col2 = fits.Column(name='POWER', unit='frac rms^2', format='E', \
		array=power)
	col3 = fits.Column(name='ERROR', unit='frac rms^2', format='E', \
		array=np.zeros(len(power)))
	cols = fits.ColDefs([col1, col2, col3])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	
	## If the file already exists, remove it (still working on just updating it)
	assert out_file[-4:].lower() == "fits", \
		'ERROR: Standard output file must have extension ".fits".'
	if os.path.isfile(out_file):
# 		print "File previously existed. Removing and rewriting."
		os.remove(out_file)
		
	## Writing the standard power spectrum to a FITS file
	thdulist = fits.HDUList([prihdu, tbhdu])
	thdulist.writeto(out_file)
## End of function 'power_out'

	
################################################################################
def powerlaw(w, beta):
    ## Gives a powerlaw of (1/w)^beta
    pl = np.where(w != 0, w ** (-beta), 0)  
    return pl

################################################################################    
def lorentzian(w, w_0, gamma):
    ## Gives a Lorentzian centered on w_0 with a FWHM of gamma
    numerator = gamma / (np.pi * 2.0)
    denominator = (w - w_0) ** 2 + (1.0/2.0 * gamma) ** 2
    L = numerator / denominator
    return L

################################################################################
def inv_frac_rms2_norm(amplitudes, dt, n_bins, mean_rate):
#     rms2_power = 2.0 * power / dt / float(n_bins) / (mean_rate ** 2)
    inv_rms2 = amplitudes * dt * n_bins * mean_rate ** 2 / 2.0
    return inv_rms2


################################################################################
def make_noise_seg(pos_freq, noise_psd_shape, dt, n_bins, noise_mean_rate):
	rand_r = np.random.standard_normal(len(pos_freq))
	if n_bins%2 == 0:
		rand_i = np.random.standard_normal(len(pos_freq)-1)
		rand_i = np.append(rand_i, 0.0) # because the nyquist frequency should only have a real value
	else:
		rand_i = np.random.standard_normal(len(pos_freq))

	## Creating the real and imaginary values from the lists of random numbers and the frequencies
	r_values = rand_r * np.sqrt(0.5 * noise_psd_shape)
	i_values = rand_i * np.sqrt(0.5 * noise_psd_shape)

	r_values[np.where(pos_freq == 0)] = 0
	i_values[np.where(pos_freq == 0)] = 0

	FT_pos = np.asarray([complex(r,i) for r,i in itertools.izip(r_values, i_values)])

	if pos_freq[0] == 0:
		FT_neg = np.conj(FT_pos[1:-1]) 
	else:
		FT_neg = np.conj(FT_pos)

	FT = np.append(FT_pos, FT_neg)
	FT = inv_frac_rms2_norm(FT, dt, n_bins, noise_mean_rate)
	
# 	return FT
	
	noise_unnorm_power = np.absolute(FT) ** 2
	noise_unnorm_power = noise_unnorm_power[0:len(pos_freq)]
	
	return noise_unnorm_power
	
# 	noise_unnorm_power[np.where(noise_unnorm_power < 0)] = 0
# 	noise_power = 2.0 * noise_unnorm_power * dt / float(n_bins) / (noise_mean_rate ** 2)
	
# 	return noise_power
	
	
################################################################################
def main(n_bins, dt, noise_mean_rate, num_seg, noise_psd_variance):
	
	###################
	## Initializations
	###################
	
	df = 1.0 / dt / float(n_bins)
	
	beta = 1.0  ## Negative slope of power law (negative is added later)
	## For QPOs, Q factor is w_0 / gamma
	w_1 = 5.4651495  ## Centroid frequency of QPO 1
	gamma_1 = 0.8835579  ## FWHM of QPO 1
	# w_2 = 300  ## Centroid frequency of QPO 2
	# gamma_2 = 30  ## FWHM of QPO 2
	
	##########################################
	## Making an array of Fourier frequencies
	##########################################
	
	frequencies = np.arange(float(-n_bins/2)+1, float(n_bins/2)+1)
	frequencies = frequencies * df # gives a freq resolution better than 1 Hz

	pos_freq = frequencies[np.where(frequencies >= 0)] # if n_bins is even, positive should have 2 more 
													   # than negative, because of the 0 freq and the 
													   # nyquist freq
	neg_freq = frequencies[np.where(frequencies < 0)]
	nyquist = pos_freq[-1]
	
# 	noise_psd_shape = noise_psd_variance * lorentzian(pos_freq, w_1, gamma_1)
	noise_psd_shape = noise_psd_variance * (lorentzian(pos_freq, w_1, gamma_1) + 0.1*powerlaw(pos_freq, beta))
	
	total_noise_lc = np.asarray([])
	total_noise_power = np.asarray([])
	
	##########################################################################
	## Making a new power spectrum from new random variables for each segment
	##########################################################################
	sum_power = np.zeros(len(pos_freq))
	for i in range(num_seg):
	
		noise_power = make_noise_seg(pos_freq, noise_psd_shape, dt, n_bins, noise_mean_rate)
		sum_power += noise_power
		
# 		noise_lc = fftpack.ifft(FT).real + noise_mean_rate
# 		total_noise_lc = np.append(total_noise_lc, noise_lc)
# 		total_noise_power = np.append(total_noise_power, noise_power)
	
	## End of for-loop through the segments
	
	power = sum_power / float(num_seg)
	
	power = 2.0 * power * dt / float(n_bins) / (noise_mean_rate ** 2)
	
# 	print np.mean(noise_lc)
# 
# 	time_bins = np.arange(len(total_noise_lc))
# 	print len(time_bins)
# 	time = time_bins * dt
	
	##########
	## Output
	##########
	power_out(pos_freq, power, dt, n_bins, nyquist, num_seg, noise_mean_rate)
	
# 	lc_out(dt, time, total_noise_lc)

## End of function 'main'


################################################################################
if __name__=='__main__':
	
	parser = argparse.ArgumentParser(description='This program makes a \
lightcurve from a Timmer and Koenig simulated power spectrum.')
	
	parser.add_argument('n_bins', type=int, help="Number of bins per segment.")
	parser.add_argument('dt', type=float, help="Time step between bins, in seconds.")
	parser.add_argument('noise_mean_rate', type=float, help="Mean count rate of the light curve.")
	parser.add_argument('num_seg', type=int, help="Number of segments to compute.")
	parser.add_argument('variance', type=float, help="Variance of the frac rms2 power spectrum.")
	
	args = parser.parse_args()

	main(args.n_bins, args.dt, args.noise_mean_rate, args.num_seg, args.variance)

################################################################################
