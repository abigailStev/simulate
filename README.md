# simulate
	
Simulates light curves with specific energy spectra for testing spectral-timing
techniques in the repos power_spectra, cross_correlation, and lag_spectra.

## Contents

### error_covariance.ipynb
An iPython notebook that tests the covariance of many runs of TimmerKoenig.py
and the resulting cross-correlation functions. Used for determining errors in
cross_correlation.

### run_simulate.sh
Bash script to run simulate_lightcurves.py.

### run_TK.sh
Bash script to run TimmerKoenig.py and plotting scripts.

### simulate_lightcurves.py
Simulates a light curve that has time-domain periodic variability an an energy 
spectrum composed of a blackbody and power law that vary out of phase.

### TimmerKoenig.py
Does a Timmer & Koenig simulation of a power spectrum, converts that to a light
curve with or without Poisson noise, then splits that up into segments and 
treats it like data for testing power_spectra and cross_correlation.



#### Disclaimer: This code come with no legal guarantees.