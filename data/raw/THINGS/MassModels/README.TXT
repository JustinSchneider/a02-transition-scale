
MASS MODELS


This directory contains the mass models as published in de Blok et al
2008 and Oh et al 2008.

File name explanation:

ISO = pseudo-isothermal halo model
NFW = NFW halo model

fix = M/L* fixed to the value predicted by stellar population
free= M/L* was a free parameter in the fit

Kr  = assumes a Kroupa IMF, all other models assume a diet-Salpeter IMF


File format: The first 12 lines (starting with "!#") contain information
on the fit parameters. A "Parameter Error" of 0.00000 (line 8) generally
means that the parameter was kept fixed during the fitting.

col 1.   Radius in kpc
col 2.   Rotation velocity of the gas component (km/s)
col 3.   Rotation velocity of the main disk component (km/s). This includes any M/L* scaling.
col 4.   Rotation velocity of the inner disk component (km/s). This includes any M/L* scaling.
col 5.   Observed rotation velocity (i.e. the input rotation curve) in km/s.
col 6.   Associated uncertainty in observed rotation velocity (km/s).
col 7.   Fitted rotation velocity associated with the dark matter halo (km/s)
col 8.   Total model rotation velocity (km/s)



