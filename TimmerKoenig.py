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
import tools

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
	

# 	fig, ax = plt.subplots()
# 	ax.plot(time_bins, lightcurve, lw=2)
# 	ax.set_xlabel("Time bins")
# 	ax.set_ylabel("Count rate")
# 	ax.set_xlim(0,np.max(time_bins))
# 	plt.savefig('./out_sim/TK_lightcurve.png', dpi=200)
# # 	plt.show()
# 	plt.close()
	
	
	out_data = np.column_stack((time, lightcurve))
	np.savetxt('./out_sim/TK_lightcurve.dat', out_data)
	
	detchans=64
	
	out_file = "./out_sim/TK_lightcurve.lc"
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
		
	Writes the first 100 elements of the ccf to a table.
	
	"""
	out_file="./TK_simulation_table.dat"
	with open(out_file, 'a') as out:
		for element in ccf[0:100,6]:
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
	pl[0] = 0
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
def inv_frac_rms2_norm(amplitudes, dt, n_bins):
	"""
	fracrms_power = power * 2.0 * dt / float(n_bins) / (mean_rate ** 2)
	if i don't put mean rate in here, it's like the mean is 1, so then i could
	multiply by whatever mean i want at the end
	"""
	
	inv_fracrms = np.zeros(len(amplitudes))
	inv_fracrms[1:] = amplitudes[1:] * n_bins  / 2.0 / dt
	inv_fracrms[0] = 0
	return inv_fracrms


################################################################################
def inv_leahy_norm(amplitudes, dt, n_bins, mean_rate):
	"""
	leahy_power = power * 2.0 * dt / float(n_bins) / mean_rate
	"""
	inv_leahy = amplitudes * n_bins * mean_rate / 2.0 / dt
	return inv_leahy
	
	
################################################################################
def make_FT_seg(pos_freq, psd_shape, dt, n_bins, FT_seed):
	"""
			make_FT_seg
	
	Makes a Fourier segment of the defined shape.
	
	"""
	np.random.seed(seed=FT_seed)  ## to ensure that we get the same 'random' 
		## numbers every realization; only want the poisson to change
	rand_r = np.random.standard_normal(len(pos_freq))
	rand_i = np.random.standard_normal(len(pos_freq)-1)
	rand_i = np.append(rand_i, 0.0)  ## because the nyquist frequency should 
		## only have a real value

	## Creating the real and imaginary values from the lists of random numbers
	r_values = rand_r * np.sqrt(0.5 * psd_shape)
	i_values = rand_i * np.sqrt(0.5 * psd_shape)

	r_values[0] = 0
	i_values[0] = 0
	
	## Putting them together to make the Fourier transform
	FT_pos = r_values + i_values*1j
	FT_neg = np.conj(FT_pos[1:-1]) 
	FT_neg = FT_neg[::-1]
	FT = np.append(FT_pos, FT_neg)

	return FT
## End of function 'make_FT_seg'
	
	
################################################################################
def main(n_bins, dt, num_seg, num_sim, psd_variance, exposure, \
	pl_scale, qpo_scale, fake_e_spec_file, psd_file, ccf_file):
	"""
			main
			
	Summary.
	
	"""
	###################
	## Initializations
	###################
	
	detchans=64
	num_seconds = dt * n_bins
# 	ks_lc_len = int(2048.0 / dt)  ## Generating 2048 seconds of lc at a time
	ks_lc_len = n_bins
	df_ks_lc = 1.0 / dt / float(ks_lc_len)
	
	beta = -1.0  ## Slope of power law (include negative here if needed)
	## For QPOs, Q factor is w_0 / gamma
# 	w_1 = 5.46710256  ## Centroid frequency of QPO
# 	gamma_1 = 0.80653875  ## FWHM of QPO
	mean = 5.40940085493
	std_dev = 0.473032436922
	FWHM = 2 * np.sqrt(2 * np.log(2))*std_dev
# 	print FWHM

	num_seeds = np.ceil(exposure / ks_lc_len / dt)
	FT_seed = np.random.randint(200, high=5000, size=num_seeds)

	#######################################
	## Reading in the fake energy spectrum
	#######################################
	
	fake_e_spec = simlc.read_fakeit_spectra(fake_e_spec_file)
	mean_countrate = np.sum(fake_e_spec) / exposure
	print "Mean countrate in channel 6:", fake_e_spec[6]/exposure
	ci_countrate = fake_e_spec / exposure
	ref_countrate = np.sum(fake_e_spec[2:26])/exposure
	
	countrate_ratio = 1170.0 / ref_countrate
	ci_countrate *= countrate_ratio
	ref_countrate *= countrate_ratio
# 	print "CI countrate:", ci_countrate
# 	print "Ref band countrate:", ref_countrate
	
	############################################################
	## Looping through the independent realizations/simulations
	############################################################
	
	for i in range(num_sim):
	
		print " "
		total_lc = np.asarray([])
		cs_sum = np.zeros((n_bins, detchans))
		sum_power_ci = np.zeros((n_bins, detchans), dtype=np.float64)
		sum_power_ref = np.zeros(n_bins, dtype=np.float64)
		sum_rate_whole_ci = np.zeros(detchans, dtype=np.float64)
		sum_rate_whole_ref = 0
		cs_sum = np.zeros((n_bins, detchans), dtype=np.complex128)
		mean_pow = np.zeros(n_bins)
		count_numseg = 0
		count_lc = 0
		mean_rate = 0
	
		########################################################################
		## Generating the light curve in 2048-second chunks and splitting it up
		########################################################################
		
		while count_numseg < num_seg:

			##########################################
			## Making an array of Fourier frequencies
			##########################################
						
			frequencies = np.arange(float(-ks_lc_len/2)+1, \
				float(ks_lc_len/2)+1) * df_ks_lc
			pos_freq = frequencies[np.where(frequencies >= 0)]  
			## positive should have 2 more than negative, because of the 0 freq
			## and the nyquist freq
			neg_freq = frequencies[np.where(frequencies < 0)]
			nyquist = np.max(np.abs(frequencies))
			
			##################################################
			## Defining the shape of the noise power spectrum
			##################################################		
						
			psd_shape = ((qpo_scale * gaussian(pos_freq,\
				mean, std_dev)) + (pl_scale * powerlaw(pos_freq, beta)))  ## this has power density units; 
			
			psd_shape_var = np.sum(psd_shape * df_ks_lc)
			
# 			print "Input PSD var:", psd_variance,"(frac rms)"
# 			print "PSD shape var:", psd_shape_var
# 			print psd_variance / psd_shape_var
			psd_shape *= (psd_variance / psd_shape_var)
			
			psd_shape = inv_frac_rms2_norm(psd_shape, dt, n_bins)
			
			#########################################################
			## Make Fourier transform using Timmer and Koenig method
			#########################################################
			
			FT = make_FT_seg(pos_freq, psd_shape, dt, ks_lc_len, \
				FT_seed[count_lc])
				
			temp_pow = np.abs(FT) ** 2
			mean_pow += temp_pow
# 			fracrms_power = temp_pow * 2.0 * dt / float(n_bins)
# 			fracrms_power = fracrms_power[0:len(pos_freq)]  ## don't subtract noise b/c no poisson noise to subtract!		
# 			fracrms_var = np.sum(fracrms_power * df_ks_lc)
# 			fracrms_rms = np.sqrt(fracrms_var)
# 			print "Var of orig FT:", fracrms_var, "(frac rms)"
# 			print "RMS of orig FT:", fracrms_rms, "(frac rms)"
						
			## Get light curve from Fourier transform
			light_curve = fftpack.ifft(FT).real
			
			###########################################################
			## Looping through the segments per simulation/realization
			###########################################################
		
			while len(light_curve) >= n_bins and count_numseg < num_seg:
				
				rate_ci = np.zeros((n_bins, detchans))
				rate_ref = np.zeros((n_bins, detchans))
				
				## Slicing this segment out of the larger 2048ks-long light curve
				this_lc = light_curve[0:n_bins]
				np.random.seed(seed=None)
				
				## Poisson stats are applied to counts, not count rate
				## Making a lightcurve to plot for checking
				
				## Adding Poisson noise to the interest and reference light curves
				for j in range(0, detchans):
					counts = (this_lc + ci_countrate[j]) * dt
					counts[np.where(counts < 0)] = 0.0
					rate_ci[:,j] = np.random.poisson(counts) / dt
				
				rate_ref = np.random.poisson((this_lc + ref_countrate) * dt) / dt
				
				## Compute cross spectrum
				cs_seg, mean_rate_seg_ci, mean_rate_seg_ref, power_ci, \
					power_ref = xcf.make_cs(rate_ci, rate_ref, n_bins, detchans)

				count_numseg += 1
				light_curve = light_curve[n_bins:]
				
				## Sums across segments -- arrays, so it adds by index
				sum_rate_whole_ci += mean_rate_seg_ci
				sum_rate_whole_ref += mean_rate_seg_ref
				sum_power_ci += power_ci
				sum_power_ref += power_ref
				cs_sum += cs_seg
# 				total_lc = np.append(total_lc, mean_lc)
# 				mean_rate += mean_rate_seg

			## End of for-loop through the segments
			count_lc += 1
			
		## End of while loop through ks lightcurves
		
# 		mean_rate /= float(num_seg)
# 		print "Mean count rate:", mean_rate
		tot_freq = fftpack.fftfreq(n_bins, d=dt)
		nyq_index = np.argmax(tot_freq) + 2  # +1 for nyquist, +1 for list cutoff 
		tot_freq = np.abs(tot_freq[0:nyq_index])
		nyq_seg = np.max(tot_freq)
		mean_pow /= float(num_seg)
		
		mean_rate_whole_ci = sum_rate_whole_ci / float(num_seg)
		mean_rate_whole_ref = sum_rate_whole_ref / float(num_seg)
		mean_power_ci = sum_power_ci / float(num_seg)
		mean_power_ref = sum_power_ref / float(num_seg)
		cs_avg = cs_sum / float(num_seg)
		
		df_seg  = 1.0 / float(num_seconds)
		## Fractional rms^2 normalization
# 		fracrms_power = mean_pow * 2.0 * dt / float(n_bins) / (ref_countrate ** 2)  
# 		fracrms_power = fracrms_power[0:nyq_index]
		
		fracrms_power = mean_power_ref * 2.0 * dt / float(n_bins) / (mean_rate_whole_ref ** 2)
		fracrms_power = fracrms_power[0:nyq_index]
		fracrms_power -= (2.0 / mean_rate_whole_ref)
		
		total_variance = np.sum(fracrms_power * df_seg)
		print "\nVariance of plotted:", total_variance, "(frac rms)"
		print "Variance of plotted:", total_variance * (ref_countrate ** 2), "(abs rms)"
		rms_total = np.sqrt(total_variance)
		print "Rms of plotted: ", rms_total, "(frac rms)"
		
		###################################
		## Computing the cross correlation
		###################################
		
		ccf_end, ccf_error = xcf.UNFILT_cs_to_ccf_w_err(cs_avg, dt, n_bins, \
			detchans, num_seconds, num_seg, mean_rate_whole_ci, \
			mean_rate_whole_ref, mean_power_ci, mean_power_ref, True)
		
		print "Mean rate for all of ci:", np.sum(mean_rate_whole_ci)
		print "Mean rate for ref:", mean_rate_whole_ref
		
		## Printing the first 100 ccf values to a table
		sim_table_out(ccf_end)
		
		t = np.arange(0, n_bins)

		##########
		## Output
		##########
		
		if i == 0:
			power_out(psd_file, tot_freq, fracrms_power, dt, n_bins, nyq_seg, num_seg,\
				mean_rate_whole_ref)
	
# 			lc_out(dt, total_lc)
	
			ccf_out(ccf_file, dt, n_bins, detchans, num_seg, mean_rate_whole_ci, \
				mean_rate_whole_ref, t, ccf_end, ccf_error)	
		## End of 'if i = 0 then do output'
		
	## End of loop through simulations/realizations
	
## End of function 'main'


################################################################################
if __name__=='__main__':
	
	############################################
	## Parsing input arguments and calling main
	############################################
	
	parser = argparse.ArgumentParser(description='This program makes a \
lightcurve from a Timmer and Koenig simulated power spectrum.', epilog='For \
optional arguments, default values are given in square brackets at end of \
description.')
	
	parser.add_argument('n_bins', type=int, help="Number of bins per segment.")
	
	parser.add_argument('dt', type=float, default=0.5, help="Time step between \
bins, in seconds. [0.5]")

	parser.add_argument('num_seg', type=int, default=1.0, help="Number of \
segments to compute. [1.0]")

	parser.add_argument('num_sim', type=int, default=1.0, help="Number of \
simulations to run. [1.0]")

	parser.add_argument('variance', type=float, default=1.0, help="Variance of \
the frac rms2 power spectrum.")

	parser.add_argument('exposure', type=float, default=10000, help="Exposure \
time of the fake energy spectrum, in seconds. [10000]")

	parser.add_argument('pl_scale', type=float, default=1e-4, help="Scale \
factor for power law component of the power spectrum. [1e-4]")

	parser.add_argument('qpo_scale', type=float, default=1e-2, help="Scale \
factor for qpo component of the power spectrum. [1e-2]")

	parser.add_argument('fake_e_spec', help="Name of fake energy spectrum.")

	parser.add_argument('psd_file', help="Root name of power spectrum output \
(fits) file.")

	parser.add_argument('ccf_file', help="Root name of CCF output (fits) file.")
	
	args = parser.parse_args()

	main(args.n_bins, args.dt, args.num_seg, args.num_sim, \
		args.variance, args.exposure, args.pl_scale, args.qpo_scale, \
		args.fake_e_spec, args.psd_file, args.ccf_file)

################################################################################
