# Vertical-DO-variance_MATLAB
This repository collects MATLAB code to compute vertical dissolved oxygen variance and its budget terms using ROMS Fennel model outputs.

To use this code, one should have result files of ROMS and biology files. 

The correct order to run the code is:
1. DOVar.m
2. O2_prod.m
3. air_sea_exchange_SOD.m
4. DOVar_bio_air_sea_SOD_NEM.m
5. get_bay_wide_xxx.m    Extract specfic terms for plotting 


Required external funtions:
- sigma2z.m          Convert ROMS sigma coordinate to z coordinate
- lanczosfilter.m    Low/high pass filter
- simps.m            Numerical integration using Simpson's rule
- derivative.m       Numerical derivative, can be used to calculate VDOV tendency term

