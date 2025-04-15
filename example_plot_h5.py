#!/usr/bin/env python
# coding: utf-8

# In[19]:


import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import GIA_assess as ga
import sys

# ## Example Output File and Plotting  
# 
# This script demonstrates how to read an HDF5 output file produced by `compute.py` and generate example plots.  
# 
# The example file [`output/example_output.h5`](output/example_output.h5) contains observed and reconstructed vertical deformation for four GNSS sites in Southeast Greenland: **MIK2, KULU, TREO,** and **NNVN**.  
# 
# - **Observed GNSS vertical deformation** data were obtained from [UNR Geodesy](http://geodesy.unr.edu/).  
# - **Reconstructed vertical deformation** was calculated using the **ANU RL02 solutions** (McGirr et al., 2024) from [ANU Data Commons](https://datacommons.anu.edu.au/DataCommons/rest/records/anudc:6133/data/Solution%20Files/ANU_mascons_RL02_2002-2023_grid.nc) and present-day GIA from:  
#     1. **ICE-6G_D** (Peltier et al., 2018) – [DOI: 10.1002/2016JB013844](https://doi.org/10.1002/2016JB013844)  
#     2. **ICE-6G_ANU** (Purcell et al., 2016) – [DOI: 10.1002/2015JB012742](https://doi.org/10.1002/2015JB012742)  
# 
# The script reads the example HDF5 file and generates two example plots.

if len(sys.argv[:]) < 2 :
    print("Sample script to plot reconstructed GPS time series and map of velocity errors over Greenland")
    print("Runstring: example_plot_h5.py file.h5")
    exit()
    
# Read the output h5 file
infile = sys.argv[1]
GPS_sites, GIA_models, recon = ga.read_reconstruction_h5(infile)


# plot time series for each GPS site
for site in GPS_sites.keys():
    # plot the GPS data
    plt.figure(figsize=(6,5))
    plt.title(f'{site}: {GPS_sites[site].lat:.2f}, {GPS_sites[site].lon:.2f}')
    plt.scatter(GPS_sites[site].decyear, GPS_sites[site].dh * 1000, s=0.5, color='grey',label='GPS')
    plt.xlabel('Time [years]')
    plt.ylabel('Vertical displacement [mm]')
    plt.xlim(2000,2026)
    # plot the reconstructed data
    for model in recon[site].keys():
        pp=plt.plot(recon[site][model].decyear, recon[site][model].dh * 1000, label=model)
        # get color of pp
        col = pp[0].get_color()
        # reconstructed elastic component
        plt.plot(recon[site][model].decyear,recon[site][model].dh_e * 1000,col,linewidth=0.5,alpha=0.5)
        # reonstructed viscoelastic component (ie the GIA model)
        plt.plot(recon[site][model].decyear,recon[site][model].dh_ve * 1000,col,linewidth=1,linestyle='--')
    plt.legend()
    plt.grid()

plt.show()


# plot a map of the residual velocities in Greenland
for model in GIA_models:
    fig, ax = plt.subplots(1,1,figsize=(5,4),
                           subplot_kw={'projection':ccrs.Stereographic(central_latitude=70, central_longitude=-42)})
    ax.set_title(f'{model}')
    ax.set_extent((-65, -15, 58, 84), crs=ccrs.PlateCarree())
    ax.coastlines()

    # create an array of site latitudes and longitudes
    lats = np.zeros(len(GPS_sites))
    lons = np.zeros(len(GPS_sites))
    residuals = np.zeros(len(GPS_sites))
    for i,site in enumerate(GPS_sites.keys()):
        lats[i] = GPS_sites[site].lat
        lons[i] = GPS_sites[site].lon
        residuals[i] = recon[site][model].dh_res

    # get vmin and vmax, centered on zero
    vrange = np.ceil(np.max(np.abs(residuals * 1000)))
    vmin, vmax = -vrange, vrange    
    sc = ax.scatter(lons, lats, s=30, c=residuals * 1000, cmap='seismic', edgecolors='k', 
                    linewidths=0.5, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), zorder=10)
    fig.colorbar(sc, ax=ax, orientation='vertical', label='Residual velocity (mm/yr)')

plt.show()


# In[ ]:




