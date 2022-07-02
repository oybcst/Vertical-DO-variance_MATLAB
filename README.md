# Vertical-DO-variance_MATLAB
This repository collects MATLAB code to compute vertical dissolved oxygen variance and its budget terms using ROMS Fennel model outputs.

To use this code, one should have result files of ROMS and biology diagonose files. 

The correct order to run the code is

DOVar.m

O2_prod.m

air_sea_exchange_SOD.m

DOVar_bio_air_sea_SOD_NEM.m

get_bay_wide_xxx.m    Extract specfic terms for plotting. 



Required external funtions to run the code

sigma2z.m          Convert ROMS sigma coordiate to z coordiate

lanczosfilter.m    Low/high pass filter

simps.m            Numerical intergration using Simpson's rule

derivative.m       Numerical dirivative, can be used to calculate VDOV tendecy term.
