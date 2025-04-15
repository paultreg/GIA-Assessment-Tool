# A tool to assess the accuracy of Glacial Isostatic Adjustment predictions of present-day crustal uplift rates
### Guadalupe Alvarez Rodriguez, Paul Tregoning and Rebecca McGirr
### _Geophysical Research Letters_, submitted April 2025.

### Abstract
> Ongoing glacial isostatic adjustment (GIA) is detectable in geodetic height time series and changes in the temporal gravity field. Global GIA models are often used to remove these signals from data but quantifying the errors in such models is difficult due to insufficient knowledge of Earth rheology and past ice history on Earth. Here we describe how estimates of Earth's temporal gravity field and observed height time series can be used to quantify errors in GIA models. We tested several GIA models in Fennoscandia, Laurentia and Greenland and found spatial coherence in the pattern of errors, with RMS velocity errors of ~1 mm/yr, ~2 mm/yr and ~5 mm/yr, respectively. Surprising, there are substantial similarities in the errors of the models tested. Our diagnostic tool can be used to identify regions where ice histories and/or Earth rheology parameters are deficient in the GIA models.


## How is it done?
By subtracting viscoelastic spherical harmonic coefficients of a GIA model from observed temporal gravity field spherical harmonics (i.e. from GRACE/FO data), it becomes possible to calculate the corresponding elastic and viscoelastic vertical deformation that would result. Summing these two and comparing the reconstructed height time series with observed height changes (i.e. from GPS height time series) permits a direct assessment of the accuracy of the GIA model.

If the reconstructed velocity is too high then the GIA model has over-estimated the viscoelastic velocity, meaning tha the model has either melted too much ice in that region or has a rheology profile that is too stiff, and vice versa if the reconstructed velocity is too low.
![Figure showing reconstructed vs observed height time series. From Alvarez Rodriguez et al (2025)](https://rses.anu.edu.au/geodynamics/GIA_assess/Fig1.jpeg)

## Software
We provide here three python scripts to permit users to perform the computations:

![compute.py](compute.py) : generates the reconstructed height time series. Inputs are 1) GRACE/FO mascon solution on 0,25 degree global grid, 2) files of GIA model(s) spherical harmonic coefficients, 3) list of GPS height time series.

![example_plot_h5.py](example_plot_h5.py): sample plotting script to visualise computations generated by compute.py

![example_reconstruction.py](example_reconstruction.py): sample computation and plotting script to generate a reconstruction at one location using one GIA model and one GRACE/FO mascon solution.

![SI_software.tar.gz](SI_software.tar.gz): tarfile containing scripts and sample input data for running these python scripts.

![GIA_assess_functions.tar](GIA_assess_functions.tar): tarfile of just the python functions required to run the above scripts

![repair_offsets.f90](repair_offsets.f90): fortran90 program to identify and repair offsets in GPS time series.


## Citation
This work is described in detail in Alvarez Rodriguez et al (2025). Please cite this paper if you use this code.


Alvarez Rodriguez, G., P. Tregoning and R. McGirr (2025), "A tool to assess the accuracy of Glacial Isostatic Adjustment predictions of present-day crustal uplift rates", _Geophysical Research Letters_, submitted April 2025.

