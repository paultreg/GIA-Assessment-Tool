#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pyshtools as pysh
import GIA_assess as ga
import sys

# ## Example Output File and Plotting  
# 
# This notebook demonstrates the use of `compute.py` to reconstruct vertical deformation for selected GPS stations. The reconstruction is based on corrected gridded GRACE/GRACE-FO mascon solutions and present-day uplift from GIA.  
# 
# In this example, we produce reconstructed vertical deformation for the **KULU** site in Southeast Greenland:  
# 
# - **Observed GNSS vertical deformation** data were obtained from [UNR Geodesy](http://geodesy.unr.edu/).  
# - **Reconstructed vertical deformation** was calculated using the **ANU RL02 solutions** (McGirr et al., 2024) from [ANU Data Commons](https://datacommons.anu.edu.au/DataCommons/rest/records/anudc:6133/data/Solution%20Files/ANU_mascons_RL02_2002-2023_grid.nc) and present-day GIA from **ICE-6G_D** (Peltier et al., 2018), available at [DOI: 10.1002/2016JB013844](https://doi.org/10.1002/2016JB013844).  
# 
# The results of the reconstruction, compared to the GPS site data, are written to an HDF5 file.

if len(sys.argv[:]) < 2 :
    print("Sample script to compute reconstructed GPS time series for a single site, single GIA model and single GRACE/FO mascon solution")
    print("Runstring: example_reconstruction.py SITE [e.g. KULU]")
    exit()

# In[ ]:
site = sys.argv[1]

# define arguments
GRACE_file = 'sample_data/GRACE_data/ANU_mascons_RL02_2002-2023_grid.nc'
#GRACE_file = 'sample_data/GRACE_data/CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'
AOD_path = 'sample_data/AOD1B/AOD1B_YYYY-MM_atm.sph'
GIA_file = 'sample_data/GIA_models/ICE-6G_D_GEOID.hs'
GPS_file = "sample_data/GPS_data/"+site+".tenv3"
lmax = 180

# define default input files
love_numbers = 'files/Load_Love2_CM.dat'
GIA_corr_file = 'files/ICE-6G_D_GRACE.hs'


# ### Step 1: Load and process GRACE data 
# - Load the corrected GRACE mascon grid
# - Convert to dimensionless stokes coefficients
# - Restore the ICE-6G_D GIA model used to correct the mascons
# - Add atmospheric pressure using monthly averaged AOD1B ATM GAA product

# In[ ]:


# load mascon grid
GRACE = ga.GRACE_data(GRACE_file)

# remove a mean field
mean_idx = (GRACE.decyear > 2004.0) * (GRACE.decyear < 2010.0)
mean_ewh = np.mean(GRACE.ewh[mean_idx], axis=0)
GRACE.ewh -= mean_ewh


# In[ ]:


# map to Driscoll and Healy grid
GRACE.map_to_DH(lmax=lmax)

# convert GRACE to spherical harmonics
GRACE.DHgrid_to_SH(lmax=lmax)

# convert SH to dimensionless stokes
Fn_ewh = ga.get_function(love_numbers,lmax=lmax,type='EWH',frame='figure')
for iepoch in range(GRACE.clm.shape[0]):
    GRACE.clm[iepoch, :, 0, 0] = 0
    GRACE.clm[iepoch, 0, 1:, :] /= Fn_ewh[1:, None]
    GRACE.clm[iepoch, 1, 1:, :] /= Fn_ewh[1:, None]


# In[ ]:


# load in ICE6G_D and add back to GRACE
GIA_corr = ga.GIA_model(GIA_corr_file,lmax=lmax)
GIA_corr.rate_to_dGIA(GRACE.decyear - np.mean(GRACE.decyear[mean_idx]))

# Add ICE6G_D back to GRACE
GRACE.clm += GIA_corr.clm


# In[ ]:


# Add AOD1B ATM monthly average (GAA)
GRACE.read_AOD1B_SH(AOD_path, lmax=lmax)
GRACE.aod_clm[:, :, 0, 0] = 0

# Add AOD1B GAA to GRACE
GRACE.clm += GRACE.aod_clm


# ### Step 2: Load and process GIA model and GPS site data
# - Read in GIA rate and convert to monthly GIA stokes coefficients at GRACE epochs
# - Read the GPS data
# 

# In[ ]:


# read in GIA rate in dimensionless stokes coefficients
GIA_models = {}
key = GIA_file.split('/')[-1].split('.')[0]
GIA_models[key] = ga.GIA_model(GIA_file,lmax=lmax)

# turn into monthly GIA stokes coefficients at GRACE epochs
GIA_models[key].rate_to_dGIA(GRACE.decyear - np.mean(GRACE.decyear[mean_idx]))


# In[ ]:


# read the GPS data, store in a dictionary
# Bec's code for a single site
GPS_sites = {}
key = GPS_file.split('/')[-1].split('.')[0]
GPS_sites[key] = ga.GPS_data(GPS_file)


# ### Step 3: Reconstruct vertical deformation at GPS site location
# - Get functions to convert dimensionless stokes to elastic vertical and viscoelastic coefficients
# - Remove GIA model from GRACE stokes and convert to viscoeleatic coefficients
# - Convert monthly GIA from dimensionless stokes to elastic vertical coefficients
# - Get reconstructed vertical deformation (elastic + viscoelastic) at GPS site lat, lon

# In[ ]:


# get the conversion factors for elastic and viscoelastic coefficients
Fn_e = ga.get_function(love_numbers,lmax=lmax,type='elastic_vertical',frame='figure')
Fn_ve = ga.get_function(love_numbers,lmax=lmax,type='viscoelastic',frame='figure')


# In[ ]:


# get GRACE clm in elastic vertical displacement 
clm_elastic = np.zeros_like(GRACE.clm)

# set up dictionary to store the reconstructed height changes
recon = {}

# loop through the GIA models, get the elastic and viscoelastic coefficients
for model, GIA in GIA_models.items(): 

    # store reconstructed height changes for each GIA model
    recon[model] = {}

    for iepoch in range(GRACE.clm.shape[0]):
        # remove the GIA model
        clm_elastic[iepoch] = GRACE.clm[iepoch] - GIA.clm[iepoch]
        # convert to elastic coefficients
        clm_elastic[iepoch, 0, 1:, :] *= Fn_e[1:, None]
        clm_elastic[iepoch, 1, 1:, :] *= Fn_e[1:, None]

    # get GIA clm in viscoelastic vertical displacement
    clm_visco = np.zeros_like(GIA.clm)
    for iepoch in range(GIA.clm.shape[0]):
        clm_visco[iepoch] = GIA.clm[iepoch]
        # convert to viscoelastic coefficients
        clm_visco[iepoch, 0, 1:, :] *= Fn_ve[1:, None]   
        clm_visco[iepoch, 1, 1:, :] *= Fn_ve[1:, None]

    # reconstruct the height changes at each GPS station
    for site, GPS in GPS_sites.items():
        # align GPS timeseries with first GRACE epoch
        GPS.GRACE_aligned(GRACE.decyear[0])
    
        # get the reconstructed dh at the GPS site
        recon[model][site] = ga.reconstruct_dh(GPS.lat, GPS.lon)
        recon[model][site].get_reconstructed_dh(GRACE.decyear, clm_elastic, clm_visco)


# compute the value of the extrapolated/interpolated GPS time series at the first epoch of the GRACE time series
# so that we can remove that offset and, therefore, align the GPS and GRACE time series in the plot


# plot a time series of the GPS data
for site in GPS_sites.keys():
    # plot the GPS data
    plt.figure(figsize=(5,4))
    plt.title(f'{site}: {GPS_sites[site].lat:.2f}, {GPS_sites[site].lon:.2f}')
    plt.scatter(GPS_sites[site].decyear, GPS_sites[site].dh * 1000, s=0.5, color='grey',label='GPS')
    plt.xlabel('Time [years]')
    plt.ylabel('Vertical displacement [mm]')
    plt.xlim(2000,2026)
    # plot the reconstructed data
    for model in recon.keys():
        plt.plot(recon[model][site].decyear, recon[model][site].dh * 1000, label='reconstructed')
        plt.plot(recon[model][site].decyear,recon[model][site].dh_e * 1000,'orange',linewidth=1.5,alpha=0.5,label='Elastic')
        plt.plot(recon[model][site].decyear,recon[model][site].dh_ve * 1000,'red',linewidth=1,linestyle='--',label='Viscoelastic')
    plt.grid()
    plt.legend()
    plt.show()


# ### Step 4: Calculate the velocity of GPS and reconstructed vertical displacement
# - Calculate rate of GPS vertical displacement over common data period
# - Calculate rate of reconstructed vertical displacement over common data period
# - Calculate the residual velocity

# In[ ]:


# calculate the velocity of GPS dh and the reconstructed dh
for site, GPS in GPS_sites.items():
    # find the common time period
    min_yr = np.max([GPS.decyear[0], GRACE.decyear[0]])
    max_yr = np.min([GPS.decyear[-1], GRACE.decyear[-1]])
    GPS_idx = (GPS.decyear > min_yr) * (GPS.decyear < max_yr)

    # get the velocity of the GPS data
    x = GPS.decyear[GPS_idx] - np.mean(GPS.decyear[GPS_idx])
    y = GPS.dh[GPS_idx]
    ### Degree-1 polynomial model:
    # GPS.model = ga.linear_regression(x,y,degree=1,annual=True,semiannual=True,verbose=False)
    # GPS.dh_vel = GPS.model.m
    ### Degree-2 polynomial model:
    GPS.model = ga.linear_regression(x,y,degree=2,annual=True,semiannual=True,verbose=False)
    GPS.dh_vel = GPS.model.b
    gps_str = f'GPS ({site})'
    print(f'{gps_str.ljust(20)}: {GPS.dh_vel * 1000:.2f} mm/yr')

    # get the velocity of the reconstructed data for each GIA model
    for model in recon.keys():
        recon_idx = (GRACE.decyear > min_yr) * (GRACE.decyear < max_yr)
        x = recon[model][site].decyear[recon_idx] - np.mean(recon[model][site].decyear[recon_idx])
        y = recon[model][site].dh[recon_idx]
        ### Degree-1 polynomial model:
        # recon[model][site].model = ga.linear_regression(x,y,degree=1,annual=True,semiannual=True,verbose=False)
        # recon[model][site].dh_vel = recon[model][site].model.m
        ### Degree-2 polynomial model:
        recon[model][site].model = ga.linear_regression(x,y,degree=2,annual=True,semiannual=True,verbose=False)
        recon[model][site].dh_vel = recon[model][site].model.b
        recon[model][site].dh_res = GPS.dh_vel - recon[model][site].dh_vel
        print(f'{model.ljust(20)}: {recon[model][site].dh_vel * 1000:.2f} mm/yr, residual: {recon[model][site].dh_res * 1000:.2f} mm/yr')


# In[ ]:


# plot a map of the residual velocities in Greenland
for model in GIA_models.keys():
    fig, ax = plt.subplots(1,1,subplot_kw={'projection':ccrs.Stereographic(central_latitude=70, central_longitude=-42)})
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
        residuals[i] = recon[model][site].dh_res

    # get vmin and vmax, centered on zero
    vrange = np.ceil(np.max(np.abs(residuals * 1000)))
    vmin, vmax = -vrange, vrange    
    sc = ax.scatter(GPS.lon, GPS.lat, s=30, c=residuals * 1000, cmap='seismic', edgecolors='k', 
                    linewidths=0.5, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), zorder=10)
    fig.colorbar(sc, ax=ax, orientation='vertical', label='Residual velocity (mm/yr)')
    plt.show()


# ### Step 5: Write results to HDF5

# In[ ]:


# Write test to h5
outfile = f'output/{site}_{model}_test.h5' 
# if file already exists, delete it
if os.path.exists(outfile):
    os.remove(outfile)
ga.write_reconstruction_to_h5(outfile, GPS_sites, recon)


# In[ ]:




