import argparse
from datetime import datetime
import os
import numpy as np
from scipy import fftpack
import matplotlib.pyplot as plt
import itertools
from astropy.io import fits
import ccf as xcf
import simulate_lightcurves as simlc

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
def lc_out(dt, lightcurve):
	"""
	
	"""
	time_bins = np.arange(len(lightcurve))
	time = time_bins * dt
	
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
def power_out(out_file, freq, power, dt, n_bins, nyquist, num_seg, mean_rate):
	"""
	
	"""
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
def ccf_out(out_file, dt, n_bins, detchans, num_segments, mean_rate_whole_ci, \
	mean_rate_whole_ref, t, ccf, ccf_error):
    """
            ccf_out

    Writes the cross-correlation function to a .fits output file.
    
    """
    num_seconds = dt * n_bins
    
    ## Getting data into a good output structure
    chan = np.arange(0, detchans)
    energy_channels = np.tile(chan, len(t))
    time_bins = np.repeat(t, detchans)
    assert len(energy_channels) == len(time_bins)
    
    print "\nOutput sent to: %s" % out_file
	
    ## Making FITS header (extension 0)
    prihdr = fits.Header()
    prihdr.set('TYPE', "Cross-correlation function generated from T&K")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', "Timmer and Koenig simulation")
    prihdr.set('BKGD', "None")
    prihdr.set('DT', dt, "seconds")
    prihdr.set('N_BINS', n_bins, "time bins per segment")
    prihdr.set('SEC_SEG', num_seconds, "seconds per segment")
    prihdr.set('SEGMENTS', num_segments, "segments in the whole light curve")
    prihdr.set('EXPOSURE', num_segments * num_seconds, \
    	"seconds, of light curve")
    prihdr.set('DETCHANS', detchans, "Number of detector energy channels")
    prihdr.set('RATE_CI', str(mean_rate_whole_ci.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_whole_ref, "counts/second")
    prihdr.set('FILTER', 'False')
    prihdu = fits.PrimaryHDU(header=prihdr)
    
    ## Making FITS table (extension 1)
    col1 = fits.Column(name='TIME_BIN', format='K', array=time_bins)
    col2 = fits.Column(name='CCF', unit='Counts/second', format='D', \
    	array=ccf.real.flatten('C'))
    col3 = fits.Column(name='ERROR', unit='', format='D', \
    	array=ccf_error.real.flatten('C'))
    col4 = fits.Column(name='CHANNEL', unit='', format='I', \
    	array=energy_channels)
    cols = fits.ColDefs([col1, col2, col3, col4])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    
    ## If the file already exists, remove it (still working on just updating it)
    assert out_file[-4:].lower() == "fits", \
    	'ERROR: Output file must have extension ".fits".'
    if os.path.isfile(out_file):
    	os.remove(out_file)
    	
    ## Writing to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file)	
## End of function 'ccf_out'


################################################################################
def sim_table_out(ccf):
	"""
			sim_table_out
		
	Writes the first 65 elements (not including 0) of the ccf to a table.
	
	"""
	out_file="./TK_simulation_table.dat"
	with open(out_file, 'a') as out:
		for element in ccf[0:100,2]:
			out.write("%.6e\t" % element.real)
		out.write("\n")
## End of function 'sim_table_out'

	
################################################################################
def powerlaw(w, beta):
	"""
    Gives a powerlaw of (1/w)^beta
    """
	pl = np.zeros(len(w))
	pl[1:] = w[1:] ** (beta)
	pl[0] = np.inf
	return pl


################################################################################    
def lorentzian(w, w_0, gamma):
	"""
    Gives a Lorentzian centered on w_0 with a FWHM of gamma
    """
	numerator = gamma / (np.pi * 2.0)
	denominator = (w - w_0) ** 2 + (1.0/2.0 * gamma) ** 2
	L = numerator / denominator
	return L


################################################################################
def gaussian(w, mean, std_dev):
	"""
	Gives a Gaussian with a mean of mean and a standard deviation of std_dev
	"""
	exp_numerator = -(w - mean)**2
	exp_denominator = 2 * std_dev**2
	G = np.exp(exp_numerator / exp_denominator)
	return G
	
	
################################################################################
def inv_frac_rms2_norm(amplitudes, dt, n_bins, mean_rate):
	"""
	rms2_power = 2.0 * power / dt / float(n_bins) / (mean_rate ** 2)
	"""
	inv_rms2 = amplitudes * dt * n_bins * mean_rate ** 2 / 2.0
	return inv_rms2


################################################################################
def inv_leahy_norm(amplitudes, dt, n_bins, mean_rate):
	"""
	leahy_power = 2.0 * power / dt / float(n_bins) / (mean_rate ** 2)
	"""
	inv_leahy = amplitudes * dt * n_bins * mean_rate / 2.0
	return inv_leahy
	
	
################################################################################
def make_noise_seg(pos_freq, noise_psd_shape, dt, n_bins, noise_mean_rate):
	"""
			make_noise_seg
	
	Makes a Fourier segment of the noise process of the defined shape.
	
	"""
	rand_r = np.random.standard_normal(len(pos_freq))
	rand_i = np.random.standard_normal(len(pos_freq)-1)
	rand_i = np.append(rand_i, 0.0)  ## because the nyquist frequency should 
									 ## only have a real value

	## Creating the real and imaginary values from the lists of random numbers
	r_values = rand_r * np.sqrt(0.5 * noise_psd_shape)
	i_values = rand_i * np.sqrt(0.5 * noise_psd_shape)

	r_values[np.where(pos_freq == 0)] = 0
	i_values[np.where(pos_freq == 0)] = 0

	FT_pos = r_values + i_values*1j
	FT_neg = np.conj(FT_pos[1:-1]) 
	FT_neg = FT_neg[::-1]

	FT = np.append(FT_pos, FT_neg)
	FT = inv_frac_rms2_norm(FT, dt, n_bins, noise_mean_rate)
# 	FT = inv_leahy_norm(FT, dt, n_bins, noise_mean_rate)
	
# 	noise_unnorm_power = np.absolute(FT) ** 2
# 	noise_unnorm_power[np.where(noise_unnorm_power < 0)] = 0

	return FT
## End of function 'make_noise_seg'
	
	
################################################################################
def main(n_bins, dt, noise_mean_rate, num_seg, noise_psd_scalefactor, pl_scale,\
	qpo_scale, fake_e_spec_file, psd_file, ccf_file):
	"""
			main
			
	Summary.
	
	"""
	###################
	## Initializations
	###################
	
	detchans=64
	df = 1.0 / dt / float(n_bins)
	num_seconds = dt * n_bins
	
	total_noise_lc = np.asarray([])
	sum_power = np.zeros(n_bins)
	cs_sum = np.zeros((n_bins, detchans))
	sum_power_ci = np.zeros((n_bins, detchans), dtype=np.float64)
	sum_power_ref = np.zeros(n_bins, dtype=np.float64)
	sum_rate_whole_ci = np.zeros(detchans, dtype=np.float64)
	sum_rate_whole_ref = 0
	cs_sum = np.zeros((n_bins, detchans), dtype=np.complex128)
	
# 	beta = -0.32240812  ## Slope of power law (include negative here if needed)
# 	beta = -7.48026728e-01
	beta = -2.95778095e-01
	## For QPOs, Q factor is w_0 / gamma
	w_1 = 5.46710256  ## Centroid frequency of QPO
	gamma_1 = 0.80653875  ## FWHM of QPO
	mean = 5.45060810e+00
	std_dev = 4.71222233e-01
	FWHM = 2 * np.sqrt(2 * np.log(2))*std_dev
# 	print FWHM
	
	fake_e_spec = simlc.read_fakeit_spectra(fake_e_spec_file)
	exposure = num_seg * num_seconds
	
	##########################################
	## Making an array of Fourier frequencies
	##########################################
	
	frequencies = np.arange(float(-n_bins/2)+1, float(n_bins/2)+1) * df
	pos_freq = frequencies[np.where(frequencies >= 0)]  
	## positive should have 2 more than negative, because of the 0 freq and the 
	## nyquist freq
	neg_freq = frequencies[np.where(frequencies < 0)]
	nyquist = np.max(np.abs(frequencies))

# 	noise_psd_shape = noise_psd_scalefactor * lorentzian(pos_freq, w_1, gamma_1)
# 	noise_psd_shape = noise_psd_scalefactor * (qpo_scale * lorentzian(pos_freq,\
# 		w_1, gamma_1) + pl_scale * powerlaw(pos_freq, beta))
	noise_psd_shape = noise_psd_scalefactor * (qpo_scale * gaussian(pos_freq,\
		mean, std_dev) + pl_scale * powerlaw(pos_freq, beta))
	
	##########################################################################
	## Making a new power spectrum from new random variables for each segment
	##########################################################################

	for i in range(num_seg):
	
		rate_ci = np.zeros((n_bins, detchans))
		FT = make_noise_seg(pos_freq, noise_psd_shape, dt, n_bins,\
			noise_mean_rate)
		
		noise_lc = fftpack.ifft(FT).real + noise_mean_rate
		noise_lc[np.where(noise_lc < 0)] = 0.0

		## Poisson stats are applied to counts, not count rate
		for i in range(0, detchans):
			rate_ci[:,i] = np.random.poisson(noise_lc * fake_e_spec[i] * dt / exposure) / dt  
		rate_ref = np.random.poisson(noise_lc * np.mean(fake_e_spec[3:32]) * dt / exposure) / dt  
# 		print np.shape(rate_ci)
# 		exit()
		
		
		cs_segment, mean_rate_segment_ci, mean_rate_segment_ref, power_ci, \
			power_ref = xcf.make_cs(rate_ci, rate_ref, n_bins, detchans)
		
		## Sums across segments -- arrays, so it adds by index
		sum_rate_whole_ci += mean_rate_segment_ci
		sum_rate_whole_ref += mean_rate_segment_ref
		sum_power_ci += power_ci
		sum_power_ref += power_ref
		cs_sum += cs_segment
		total_noise_lc = np.append(total_noise_lc, noise_lc)
	## End of for-loop through the segments

	mean_rate_whole_ci = sum_rate_whole_ci / float(num_seg)
	mean_rate_whole_ref = sum_rate_whole_ref / float(num_seg)
	mean_power_ci = sum_power_ci / float(num_seg)
	mean_power_ref = sum_power_ref / float(num_seg)
	cs_avg = cs_sum / float(num_seg)

	power = sum_power_ref / float(num_seg)
	power *= 2.0 * dt / float(n_bins) / (noise_mean_rate ** 2)  ## Fractional rms^2 normalization
	power -= 2.0 / noise_mean_rate
	power = power[0:len(pos_freq)]
# 	power = 2.0 * power * dt / float(n_bins) / (noise_mean_rate)  ## Leahy normalization
	
	ccf_end, ccf_error = xcf.UNFILT_cs_to_ccf_w_err(cs_avg, dt, n_bins, \
    	detchans, num_seconds, num_seg, mean_rate_whole_ci, \
    	mean_rate_whole_ref, mean_power_ci, mean_power_ref, True)
	
	t = np.arange(0, n_bins)

	##########
	## Output
	##########
	
	power_out(psd_file, pos_freq, power, dt, n_bins, nyquist, num_seg, \
		noise_mean_rate)
	
	lc_out(dt, total_noise_lc)
	
	ccf_out(ccf_file, dt, n_bins, detchans, num_seg, mean_rate_whole_ci, \
		mean_rate_whole_ref, t, ccf_end, ccf_error)	
	
	sim_table_out(ccf_end)
		
## End of function 'main'


################################################################################
if __name__=='__main__':
	
	############################################
	## Parsing input arguments and calling main
	############################################
	parser = argparse.ArgumentParser(description='This program makes a \
lightcurve from a Timmer and Koenig simulated power spectrum.')
	
	parser.add_argument('n_bins', type=int, help="Number of bins per segment.")
	parser.add_argument('dt', type=float, default=0.5, help="Time step between bins, in seconds. [0.5]")
	parser.add_argument('noise_mean_rate', type=float, default=1000.0, help="Mean count rate of the light curve. [1000]")
	parser.add_argument('num_seg', type=int, default=1.0, help="Number of segments to compute. [1.0]")
	parser.add_argument('variance', type=float, default=1.0, help="Variance of the frac rms2 power spectrum.")
	parser.add_argument('pl_scale', type=float, default=1e-4, help="Scale factor for power law component of the power spectrum. [1e-4]")
	parser.add_argument('qpo_scale', type=float, default=1e-2, help="Scale factor for qpo component of the power spectrum. [1e-2]")
	parser.add_argument('fake_e_spec', help="Name of fake energy spectrum.")
	parser.add_argument('psd_file', help="Name of power spectrum output (fits) file.")
	parser.add_argument('ccf_file', help="Name of CCF output (fits) file.")
	
	args = parser.parse_args()

	main(args.n_bins, args.dt, args.noise_mean_rate, args.num_seg, \
		args.variance, args.pl_scale, args.qpo_scale, args.fake_e_spec, \
		args.psd_file, args.ccf_file)

################################################################################
