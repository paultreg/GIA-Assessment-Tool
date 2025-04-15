'''
Program to compute reconstructed height changes from corrected GRACE mascons and GIA models and
compares these to GPS height data for a list of GIA models and GPS stations.

Author: Rebecca McGirr, 2025
Research School of Earth Sciences
The Australian National University
'''

import os
import argparse 
import numpy as np
import GIA_assess as ga

############################################################################### 
runstr = "compute.py --GRACE corrected_mascon_solution.nc --AOD_path AOD1B_YYYY-MM_atm.sph --GIA_list GIA_filenames --GPS_list GPS_filenames.txt --lmax 180 --out outfile"
parser = argparse.ArgumentParser(prog='compute.py',formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''This program computes reconstructed height changes from corrected GRACE mascons and GIA models and
                                                compares these to atmospheric pressure corrected GPS data for a list of GIA models and GPS stations.
Runstring: %s'''%runstr)

# Define arguments
parser.add_argument("--GRACE",default='sample_data/GRACE_data/ANU_mascons_RL02_2002-2023_grid.nc',type=str,help="netcdf file containing corrected GRACE mascon solution")
parser.add_argument("--AOD_path",default='sample_data/AOD1B/AOD1B_YYYY-MM_atm.sph',type=str,help="filepath to monthly average ATM AOD1B files (GAA)")
parser.add_argument("--GIA_list",default='GIA_models.txt',type=str,help="text file containing list of GIA files in dimensionless SH coefficients")
parser.add_argument("--GPS_list",default='GPS_data.txt',type=str,help="text file containing list of GPS tenv3 files corrected for offsets")
parser.add_argument("--lmax",default=180,type=int,help="maximum degree of SH coefficients [default=180]")
parser.add_argument("--out",default='output/out.h5',type=str,help="output h5 file path [default=output/out.h5]")

# read in the arguments
args = parser.parse_args()
GRACE_file = args.GRACE
AOD_path = args.AOD_path
GPS_list = args.GPS_list
GIA_list = args.GIA_list
lmax = args.lmax
outfile = args.out

###############################################################################
## DEFINE DEFAULT INPUT FILES
###############################################################################

# Love number file
love_numbers = 'files/Load_Love2_CM.dat'
# GIA correction stokes coefficients (default ICE6G_D_GRACE.hs)
GIA_corr_file = 'files/ICE-6G_D_GRACE.hs'

###############################################################################
## LOAD AND PROCESS GRACE DATA 
###############################################################################

# Load the GRACE data
GRACE = ga.GRACE_data(GRACE_file)

# remove a mean field
mean_idx = (GRACE.decyear > 2004.0) * (GRACE.decyear < 2010.0)
mean_ewh = np.mean(GRACE.ewh[mean_idx], axis=0)
GRACE.ewh -= mean_ewh

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


# load in GIA model used to correct GRACE mascons
if not os.path.exists(GIA_corr_file):
    print("File not found")
    exit()
else:
    GIA_corr = ga.GIA_model(GIA_corr_file,lmax=lmax)
    GIA_corr.rate_to_dGIA(GRACE.decyear - np.mean(GRACE.decyear[mean_idx]))

# Add ICE6G_D back to GRACE
GRACE.clm += GIA_corr.clm

# Add AOD1B ATM monthly average (GAA)
GRACE.read_AOD1B_SH(AOD_path, lmax=lmax)
GRACE.aod_clm[:, :, 0, 0] = 0
GRACE.clm += GRACE.aod_clm

###############################################################################
## LOAD AND PROCESS GIA MODELS 
###############################################################################

# Read the list of GIA model file locations
with open(GIA_list) as f:
    GIA_files = f.readlines()
GIA_files = [x.strip() for x in GIA_files]

# read the GIA models, store in a dictionary
GIA_models = {}
for fname in GIA_files:
    key = fname.split('/')[-1].split('.')[0]
    # read in GIA rate in dimensionless stokes coefficients
    GIA_models[key] = ga.GIA_model(fname,lmax=lmax)
    # turn into monthly GIA stokes coefficients at GRACE epochs
    GIA_models[key].rate_to_dGIA(GRACE.decyear - np.mean(GRACE.decyear[mean_idx]))

###############################################################################
## LOAD AND PROCESS GPS DATA 
###############################################################################

# Read the list of GPS data file locations
with open(GPS_list) as f:
    GPS_files = f.readlines()
GPS_files = [x.strip() for x in GPS_files]

# read the GPS data, store in a dictionary
GPS_sites = {}
for fname in GPS_files:
    key = fname.split('/')[-1].split('.')[0]
    GPS_sites[key] = ga.GPS_data(fname)

###############################################################################
## RECONSTRUCT VERTICAL DEFORMATION
###############################################################################

# get the conversion factors for elastic and viscoelastic coefficients
Fn_e = ga.get_function(love_numbers,lmax=lmax,type='elastic_vertical',frame='figure')
Fn_ve = ga.get_function(love_numbers,lmax=lmax,type='viscoelastic',frame='figure')

# get GRACE clm in elastic vertical displacement 
clm_elastic = np.zeros_like(GRACE.clm)

# set up dictionary to store the reconstructed height changes
recon = {}

# loop through the GIA models, get the elastic and viscoelastic coefficients
print("GRACE Solution: ",GRACE_file)
for model, GIA in GIA_models.items(): 
    print("Computing reconstructed heights using GIA model: ",model)
    
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
        # remove the modelled dh at the first GRACE epoch
        GPS.GRACE_aligned(GRACE.decyear[0])
    
        # get the reconstructed dh at the GPS site
        recon[model][site] = ga.reconstruct_dh(GPS.lat, GPS.lon)
        recon[model][site].get_reconstructed_dh(GRACE.decyear, clm_elastic, clm_visco)
        
###############################################################################
## CALCULATE VELOCITY OF GPS AND RECONSTRUCTED DH
###############################################################################

# calculate the velocity of GPS dh and the reconstructed dh
#RMS_values = np.zeros((len(GPS_sites.items),len(recon.keys)))

nsites = 0
for site, GPS in GPS_sites.items():
    
    # only model common data period
    min_yr = np.max([GPS.decyear[0], GRACE.decyear[0]])
    max_yr = np.min([GPS.decyear[-1], GRACE.decyear[-1]])
    GPS_idx = (GPS.decyear > min_yr) * (GPS.decyear < max_yr)

    # get the velocity of the GPS data
    x = GPS.decyear[GPS_idx] - np.mean(GPS.decyear[GPS_idx])
    y = GPS.dh[GPS_idx]
    GPS.model = ga.linear_regression(x,y,degree=2,annual=True,semiannual=True,verbose=False)
    GPS.dh_vel = GPS.model.b
    gps_str = f'GPS ({site})'
    print(f'{gps_str.ljust(20)}     : {GPS.dh_vel * 1000:6.2f} mm/yr')
    

    # get the velocity of the reconstructed data for each GIA model
    for model in recon.keys():
        recon_idx = (GRACE.decyear > min_yr) * (GRACE.decyear < max_yr)
        x = recon[model][site].decyear[recon_idx] - np.mean(recon[model][site].decyear[recon_idx])
        y = recon[model][site].dh[recon_idx]
        # 
        recon[model][site].model = ga.linear_regression(x,y,degree=2,annual=True,semiannual=True,verbose=False)
        recon[model][site].dh_vel = recon[model][site].model.b
        recon[model][site].dh_res = GPS.dh_vel - recon[model][site].dh_vel
        print(f'{model.ljust(20)} {site}: {recon[model][site].dh_vel * 1000:6.2f} mm/yr, residual: {recon[model][site].dh_res * 1000:6.2f} mm/yr')

        # add the square of the residual velocity
        GIA_models[model].model_rms += (recon[model][site].dh_res * 1000)**2
        #imodel += 1
        
    nsites += 1

# print RMS values
for model in recon.keys():
    GIA_models[model].model_rms = np.sqrt(GIA_models[model].model_rms/nsites)
    print(GRACE_file," Velocity RMS for GIA model ",model," : ",GIA_models[model].model_rms)
        
###############################################################################
## WRITE OUT RESULTS
###############################################################################

# if file already exists, delete it
if os.path.exists(outfile):
    os.remove(outfile)
# write to file
ga.write_reconstruction_to_h5(outfile, GPS_sites, recon)
