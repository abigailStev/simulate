# simulate
	
Simulates interesting light curves with specific energy spectra for testing 
spectral-timing techniques in the repos power_spectra, cross_correlation, and 
lag_spectra. Please see [Stevens & Uttley 2016](https://ui.adsabs.harvard.edu/#abs/2016MNRAS.460.2796S/abstract)
for reference.

The functionality of this software will be folded into [Stingray](http://stingraysoftware.github.io/),
so please get involved over there!

## Contents

### error_covariance.ipynb
An iPython notebook that tests the covariance of many runs of TimmerKoenig.py
and the resulting cross-correlation functions. Used for determining errors in
cross_correlation.

### fake_qpo_spectra.py
Simulates a QPO from a set of energy spectral parameters and computes power 
spectrum, cross-correlation, and lag-energy spectrum. The input file 
'parfit_file' lists the model string in the first column, then either the 
tied energy spectral fit parameter or the three parameters of the best-fitting 
sine wave to the variation with QPO phase in the form [amplitude,phase,y-offset] 
(with no spaces).

### LICENSE.md
The software in this repository is licensed under the MIT license. Details can
be found in this document.

### run_fake_qpo_spectra.sh
Bash script to run fake_qpo_spectra.py.

### run_simulate.sh
Bash script to run simulate_lightcurves.py.

### run_TK.sh
Bash script to run TimmerKoenig.py and plotting scripts.

### simulate_lightcurves.py
Simulates a light curve that has time-domain periodic variability an an energy 
spectrum composed of a blackbody and power law that vary out of phase.

### TimmerKoenig.py
Does a [Timmer & Koenig](https://ui.adsabs.harvard.edu/#abs/1995A&A...300..707T/abstract)
simulation of a power spectrum, converts that to a light curve with or without
Poisson noise, then splits that up into segments and treats it like data for
testing power_spectra and cross_correlation.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Authors
* Abigail Stevens (UvA)

## Copyright

All content © 2014-2017 the Authors, and is distributed under the MIT
license. See LICENSE.md for details.

If you use this code, please cite [Stevens & Uttley 2016](https://ui.adsabs.harvard.edu/#abs/2016MNRAS.460.2796S/abstract).

The functionality of this software will be folded into [Stingray](http://stingraysoftware.github.io/),
so please get involved over there!