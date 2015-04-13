		simulate
	
Two main programs in here: one simulates a light curve that has time-domain 
periodic variability an an energy spectrum composed of a blackbody and power 
law that vary out of phase. The other does a Timmer & Koenig simulation of a 
power spectrum, converts that to a light curve w/ Poisson noise, then splits 
that up into segments and treats it like data for testing the ccf and errors.

CONTENTS:

simulate_lightcurves.py -- Simulates a light curve that has time-domain 
periodic variability an an energy spectrum composed of a blackbody and power 
law that vary out of phase

run_simulate.sh -- Bash script to run simulate_lightcurves.py.

TimmerKoenig.py -- Does a Timmer & Koenig simulation of a power spectrum, 
converts that to a light curve w/ Poisson noise, then splits that up into 
segments and treats it like data for testing the ccf and errors.

run_TK.sh -- Bash script to run TimmerKoenig.py and plotting scripts.


THIS CODE COMES WITH NO GUARANTEES.