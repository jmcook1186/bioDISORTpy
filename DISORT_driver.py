#!/usr/bin/env python
"""
author: J Cook, Feb 2020

Driver script for Python wrapper to the DISORT library that allows SNICAR style input parameter
definitions, calculation of optical thicknesses for mixed ice/impurity layers and operation
over a wavelength range rather than monochrome calculations.

SNICAR inputs re used to define an optical thickness, ssa and gg for each vertical layer - these
are then fed into DISORT.

Module '_disort' is auto-generated with f2py (version:2).

RUN THIS SCRIPT FROM THE TERMINAL!

i.e
cd /home/joe/Code/pyDISORT/test/
python DISORT_driver.py


"""

import disort
import numpy as np
import matplotlib.pyplot as plt
import setup_SNICAR
import pandas as pd
import collections as c


inputs = c.namedtuple('inputs',['dir_base',\
    'rf_ice', 'incoming_i', 'DIRECT', 'layer_type',\
    'APRX_TYP', 'DELTA', 'solzen', 'R_sfc', 'dz', 'rho_layers', 'grain_rds',\
    'side_length', 'depth', 'rwater', 'nbr_lyr', 'nbr_aer', 'grain_shp', 'shp_fctr', 'grain_ar', 'GA_units',\
    'Cfactor','mss_cnc_soot1', 'mss_cnc_soot2', 'mss_cnc_brwnC1', 'mss_cnc_brwnC2', 'mss_cnc_dust1',\
    'mss_cnc_dust2', 'mss_cnc_dust3', 'mss_cnc_dust4', 'mss_cnc_dust5', 'mss_cnc_ash1', 'mss_cnc_ash2',\
    'mss_cnc_ash3', 'mss_cnc_ash4', 'mss_cnc_ash5', 'mss_cnc_ash_st_helens', 'mss_cnc_Skiles_dust1', 'mss_cnc_Skiles_dust2',\
    'mss_cnc_Skiles_dust3', 'mss_cnc_Skiles_dust4', 'mss_cnc_Skiles_dust5', 'mss_cnc_GreenlandCentral1',\
    'mss_cnc_GreenlandCentral2', 'mss_cnc_GreenlandCentral3', 'mss_cnc_GreenlandCentral4',\
    'mss_cnc_GreenlandCentral5', 'mss_cnc_Cook_Greenland_dust_L', 'mss_cnc_Cook_Greenland_dust_C',\
    'mss_cnc_Cook_Greenland_dust_H', 'mss_cnc_snw_alg', 'mss_cnc_glacier_algae', 'FILE_soot1',\
    'FILE_soot2', 'FILE_brwnC1', 'FILE_brwnC2', 'FILE_dust1', 'FILE_dust2', 'FILE_dust3', 'FILE_dust4', 'FILE_dust5',\
    'FILE_ash1', 'FILE_ash2', 'FILE_ash3', 'FILE_ash4', 'FILE_ash5', 'FILE_ash_st_helens', 'FILE_Skiles_dust1', 'FILE_Skiles_dust2',\
    'FILE_Skiles_dust3', 'FILE_Skiles_dust4', 'FILE_Skiles_dust5', 'FILE_GreenlandCentral1',\
    'FILE_GreenlandCentral2', 'FILE_GreenlandCentral3', 'FILE_GreenlandCentral4', 'FILE_GreenlandCentral5',\
    'FILE_Cook_Greenland_dust_L', 'FILE_Cook_Greenland_dust_C', 'FILE_Cook_Greenland_dust_H', 'FILE_snw_alg', 'FILE_glacier_algae'])


##############################
## 2) Set working directory 
##############################

# set dir_base to the location of the BioSNICAR_GO_PY folder
inputs.dir_base = '/home/tothepoles/Desktop/bioDISORTpy/'
savepath = inputs.dir_base # base path for saving figures

################################
## 3) Choose plot/print options
################################

show_figs = False # toggle to display spectral albedo figure
save_figs = True # toggle to save spectral albedo figure to file
print_BBA = True # toggle to print broadband albedo to terminal
print_band_ratios = False # toggle to print various band ratios to terminal
smooth = False # apply optional smoothing function (Savitzky-Golay filter)
window_size = 9 # if applying smoothing filter, define window size
poly_order = 3 # if applying smoothing filter, define order of polynomial

#######################################
## 4) RADIATIVE TRANSFER CONFIGURATION
#######################################

inputs.DIRECT   = 1       # 1= Direct-beam incident flux, 0= Diffuse incident flux
inputs.APRX_TYP = 1        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
inputs.DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
inputs.solzen   = 40      # if DIRECT give solar zenith angle between 0 and 89 degrees (from 0 = nadir, 90 = horizon)


# CHOOSE ATMOSPHERIC PROFILE for surface-incident flux:
#    0 = mid-latitude winter
#    1 = mid-latitude summer
#    2 = sub-Arctic winter
#    3 = sub-Arctic summer
#    4 = Summit,Greenland (sub-Arctic summer, surface pressure of 796hPa)
#    5 = High Mountain (summer, surface pressure of 556 hPa)
#    6 = Top-of-atmosphere
# NOTE that clear-sky spectral fluxes are loaded when direct_beam=1,
# and cloudy-sky spectral fluxes are loaded when direct_beam=0
inputs.incoming_i = 2


inputs.dz = [10, 10] # thickness of each vertical layer (unit = m)
inputs.nbr_lyr = len(inputs.dz)  # number of snow layers
inputs.layer_type = [1,1] # Fresnel layers for the ADD_DOUBLE option, set all to 0 for the TOON option
inputs.rho_layers = [700, 700] # density of each layer (unit = kg m-3) 
inputs.nbr_wvl=480 
#inputs.R_sfc = np.array([0.1 for i in range(inputs.nbr_wvl)]) # reflectance of underlying surface - set across all wavelengths
inputs.R_sfc = np.genfromtxt('./Data/rain_polished_ice_spectrum.csv', delimiter = 'csv') # import underlying ice from file

###############################################################################
## 5) SET UP OPTICAL & PHYSICAL PROPERTIES OF SNOW/ICE GRAINS
# For hexagonal plates or columns of any size choose GeometricOptics
# For sphere, spheroids, koch snowflake with optional water coating choose Mie
###############################################################################

inputs.rf_ice = 2 # define source of ice refractive index data. 0 = Warren 1984, 1 = Warren 2008, 2 = Picard 2016
# Ice grain shape can be 0 = sphere, 1 = spheroid, 2 = hexagonal plate, 3 = koch snowflake, 4 = hexagonal prisms
# For 0,1,2,3:
inputs.grain_shp =[0,0] # grain shape(He et al. 2016, 2017)
inputs.grain_rds = [130,130] # effective grain radius of snow/bubbly ice (becomes bubble rds when layer_type==1)
inputs.rwater = [0, 0] # radius of optional liquid water coating

# For 4:
inputs.side_length = [10000,10000] 
inputs.depth = [10000,10000]


# Shape factor = ratio of nonspherical grain effective radii to that of equal-volume sphere
### only activated when sno_shp > 1 (i.e. nonspherical)
### 0=use recommended default value (He et al. 2017)
### use user-specified value (between 0 and 1)
inputs.shp_fctr = [0,0] 

# Aspect ratio (ratio of width to length)
inputs.grain_ar = [0,0] 

#######################################
## 5) SET LAP CHARACTERISTICS
#######################################

# Define total number of different LAPs/aerosols in model
inputs.nbr_aer = 30

# define units for glacier algae MAC input file
# 0 = m2/kg
# 1 = m2/cell
inputs.GA_units = 0

# determine C_factor (can be None or a number)
# this is the concentrating factor that accounts for
# resolution difference in field samples and model layers
inputs.Cfactor = 10

# Set names of files containing the optical properties of these LAPs:
inputs.FILE_soot1  = 'mie_sot_ChC90_dns_1317.nc'
inputs.FILE_soot2  = 'miecot_slfsot_ChC90_dns_1317.nc'
inputs.FILE_brwnC1 = 'brC_Kirch_BCsd.nc'
inputs.FILE_brwnC2 = 'brC_Kirch_BCsd_slfcot.nc'
inputs.FILE_dust1  = 'dust_balkanski_central_size1.nc'
inputs.FILE_dust2  = 'dust_balkanski_central_size2.nc'
inputs.FILE_dust3  = 'dust_balkanski_central_size3.nc'
inputs.FILE_dust4  = 'dust_balkanski_central_size4.nc'
inputs.FILE_dust5 = 'dust_balkanski_central_size5.nc'
inputs.FILE_ash1  = 'volc_ash_eyja_central_size1.nc'
inputs.FILE_ash2 = 'volc_ash_eyja_central_size2.nc'
inputs.FILE_ash3 = 'volc_ash_eyja_central_size3.nc'
inputs.FILE_ash4 = 'volc_ash_eyja_central_size4.nc'
inputs.FILE_ash5 = 'volc_ash_eyja_central_size5.nc'
inputs.FILE_ash_st_helens = 'volc_ash_mtsthelens_20081011.nc'
inputs.FILE_Skiles_dust1 = 'dust_skiles_size1.nc'
inputs.FILE_Skiles_dust2 = 'dust_skiles_size2.nc'
inputs.FILE_Skiles_dust3 = 'dust_skiles_size3.nc'
inputs.FILE_Skiles_dust4 = 'dust_skiles_size4.nc'
inputs.FILE_Skiles_dust5 = 'dust_skiles_size5.nc'
inputs.FILE_GreenlandCentral1 = 'dust_greenland_central_size1.nc'
inputs.FILE_GreenlandCentral2 = 'dust_greenland_central_size2.nc'
inputs.FILE_GreenlandCentral3 = 'dust_greenland_central_size3.nc'
inputs.FILE_GreenlandCentral4 = 'dust_greenland_central_size4.nc'
inputs.FILE_GreenlandCentral5  = 'dust_greenland_central_size5.nc'
inputs.FILE_Cook_Greenland_dust_L = 'dust_greenland_Cook_LOW_20190911.nc'
inputs.FILE_Cook_Greenland_dust_C = 'dust_greenland_Cook_CENTRAL_20190911.nc'
inputs.FILE_Cook_Greenland_dust_H = 'dust_greenland_Cook_HIGH_20190911.nc'
inputs.FILE_snw_alg  = 'snw_alg_r025um_chla020_chlb025_cara150_carb140.nc'
inputs.FILE_glacier_algae = 'Cook2020_glacier_algae_4_40.nc'


# Indicate mass mixing ratios scenarios for each impurity (units: ng(species)/g(ice), or ppb)
# glacier algae in cells/mL if GA_units ==1, ppb if GA_units == 0.
# The script will loop over the different mixing scenarios


inputs.mss_cnc_soot1 = [0]*len(inputs.dz)    # uncoated black carbon (Bohren and Huffman, 1983)
inputs.mss_cnc_soot2 = [0]*len(inputs.dz)    # coated black carbon (Bohren and Huffman, 1983)
inputs.mss_cnc_brwnC1 = [0]*len(inputs.dz)   # uncoated brown carbon (Kirchstetter et al. (2004).)
inputs.mss_cnc_brwnC2 = [0]*len(inputs.dz)   # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
inputs.mss_cnc_dust1 = [0]*len(inputs.dz)    # dust size 1 (r=0.05-0.5um) (Balkanski et al 2007)
inputs.mss_cnc_dust2 = [0]*len(inputs.dz)    # dust size 2 (r=0.5-1.25um) (Balkanski et al 2007)
inputs.mss_cnc_dust3 = [0]*len(inputs.dz)    # dust size 3 (r=1.25-2.5um) (Balkanski et al 2007)
inputs.mss_cnc_dust4 = [0]*len(inputs.dz)    # dust size 4 (r=2.5-5.0um)  (Balkanski et al 2007)
inputs.mss_cnc_dust5 = [0]*len(inputs.dz)    # dust size 5 (r=5.0-50um)  (Balkanski et al 2007)
inputs.mss_cnc_ash1 = [0]*len(inputs.dz)    # volcanic ash size 1 (r=0.05-0.5um) (Flanner et al 2014)
inputs.mss_cnc_ash2 = [0]*len(inputs.dz)    # volcanic ash size 2 (r=0.5-1.25um) (Flanner et al 2014)
inputs.mss_cnc_ash3 = [0]*len(inputs.dz)    # volcanic ash size 3 (r=1.25-2.5um) (Flanner et al 2014)
inputs.mss_cnc_ash4 = [0]*len(inputs.dz)    # volcanic ash size 4 (r=2.5-5.0um) (Flanner et al 2014)
inputs.mss_cnc_ash5 = [0]*len(inputs.dz)    # volcanic ash size 5 (r=5.0-50um) (Flanner et al 2014)
inputs.mss_cnc_ash_st_helens = [0]*len(inputs.dz)   # ashes from Mount Saint Helen's
inputs.mss_cnc_Skiles_dust1 = [0]*len(inputs.dz)   # Colorado dust size 1 (Skiles et al 2017)
inputs.mss_cnc_Skiles_dust2 = [0]*len(inputs.dz)    # Colorado dust size 2 (Skiles et al 2017)
inputs.mss_cnc_Skiles_dust3 = [0]*len(inputs.dz)    # Colorado dust size 3 (Skiles et al 2017)
inputs.mss_cnc_Skiles_dust4 = [0]*len(inputs.dz)  # Colorado dust size 4 (Skiles et al 2017)
inputs.mss_cnc_Skiles_dust5 = [0]*len(inputs.dz)  # Colorado dust size 5 (Skiles et al 2017)
inputs.mss_cnc_GreenlandCentral1 = [0]*len(inputs.dz) # Greenland Central dust size 1 (Polashenski et al 2015)
inputs.mss_cnc_GreenlandCentral2 = [0]*len(inputs.dz) # Greenland Central dust size 2 (Polashenski et al 2015)
inputs.mss_cnc_GreenlandCentral3 = [0]*len(inputs.dz) # Greenland Central dust size 3 (Polashenski et al 2015)
inputs.mss_cnc_GreenlandCentral4 = [0]*len(inputs.dz) # Greenland Central dust size 4 (Polashenski et al 2015)
inputs.mss_cnc_GreenlandCentral5 = [0]*len(inputs.dz) # Greenland Central dust size 5 (Polashenski et al 2015)
inputs.mss_cnc_Cook_Greenland_dust_L = [0]*len(inputs.dz) # GRIS dust (Cook et al. 2019 "LOW")
inputs.mss_cnc_Cook_Greenland_dust_C = [0]*len(inputs.dz) # GRIS dust 1 (Cook et al. 2019 "mean")
inputs.mss_cnc_Cook_Greenland_dust_H = [0]*len(inputs.dz) # GRIS dust 1 (Cook et al. 2019 "HIGH")
inputs.mss_cnc_snw_alg = [0]*len(inputs.dz)    # Snow Algae (spherical, C nivalis) (Cook et al. 2017)
inputs.mss_cnc_glacier_algae = [0.000001,0]    # glacier algae in cells/ml or ppb depending on GA_units (Cook et al. 2020)



# DISORT CONFIG
prnt = np.array([True, True, True, True, True]) # determines what info to print to console
umu = [0.2,0.4,0.6,0.8,1.0]  # cosine of viewing zenith angle
coszen = [0.57] #[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] # cosine of solar zenith angle
phi0 = 0  # solar azimuth angle
phi = 0 # viewing azimuth angle
uTau = 0.0  # optical thickness where fluxes are calculated
Nstr = 16 # number of streams to include in model

# RETRIEVE OPTICAL THICKNESS, SSA, ASYMMETRY PARAMETER FROM SNICAR

[tau, g, SSA, mu_not, nbr_wvl, wvl, L_snw, flx_slr] = setup_SNICAR.setup_SNICAR(inputs)

fbeam = flx_slr # incoming irradiance (output by SNICAR)

# G_STAR, SSA_STAR, G_STAR, fbeam are all arrays generated by the snicar function
# containing values at each wavelength. Each call to DISORT takes
# a single value representing one wavelength. The monochrome versions
# of each variable are aliased gg, w0 and dTau for feeding into disort.

def run_disort(gg, w0, dTau, umu0, umu, phi0, phi, albedo, uTau, fbeam, prnt):

    [rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed] =\
                                    disort.run(dTau = dTau, w0=w0, iphas=iphas, gg=gg, temp=273.15,
                                            umu0=umu0, phi0=phi0, albedo=albedo, fbeam=fbeam,
                                                utau=uTau, umu=umu, phi=phi, prnt=prnt, UsrAng=True,
                                                    onlyFl=True, Nstr=Nstr)
                                                    
    alb = flup/(rfldir+rfldn) # calculate albedo from upwards and downards fluxes

    return alb

### SET UP OUTPUT ARRAYS
alb_out = np.zeros(shape=(len(coszen),len(wvl)))
BBA_out = np.zeros(shape=(len(coszen),1))
counter = 0

#################################
# CALL DISORT
#################################

for ang in coszen:   # iterate over cosine of solar zenith angle
    
    for i in range(len(wvl)): # iterate over wavelength

        dTau = tau[:,i]
        w0 = SSA[:,i]
        iphas = np.ones(len(inputs.dz),dtype='int')*3
        #iphas[:] = 3 # 3 = use Henyey-Greenstein asymmetry param
        gg = g[:,i] # define asymmetry param
        albedo = inputs.R_sfc[i]
        umu0 = ang
        alb = run_disort(gg, w0, dTau, umu0, umu, phi0, phi, albedo, uTau, fbeam, prnt)
        alb_out[counter,i] = alb
        

    # calculcate broadband albedo
    alb_mean = np.mean(alb_out,axis=0)
        
    albDF = pd.DataFrame(columns=['Alb'])
    albDF['Alb'] = alb_out[counter,:]
    albDF['Alb'].loc[albDF['Alb']<0]=np.nan
    albDF['Alb'] = albDF['Alb'].interpolate(axis=0)

    alb_out[counter,:] = albDF['Alb']
    BBA_out[counter,:] = np.sum(albDF['Alb'].loc[albDF['Alb']>=0]*flx_slr[albDF['Alb']>=0])/np.sum(flx_slr[albDF['Alb']>=0])

    counter+=1

################################
# PLOT FIGURES
################################

# plot albedo at each solar zenith
plt.figure()
alb_mean = np.mean(alb_out,axis=0)

smooth = False

if smooth:
    from scipy.signal import savgol_filter
    yhat = savgol_filter(alb_mean, 15, 3)
    alb_mean = yhat


for i in range(len(coszen)):
    plt.plot(wvl[0:270],alb_out[i,0:270],label='coszen: {}'.format(str(coszen[i])))

plt.plot(wvl[0:270],alb_mean[0:270],linestyle='dashed',color='k',label='zenMean')
plt.ylabel('Albedo')
plt.xlabel('Wavelength')
plt.xlim(0.2,2.5)
plt.legend()
plt.savefig(str(inputs.dir_base+'albedo_DISORT.png'))

if print_band_ratios:

    I2DBA = alb_mean[40]/alb_mean[36]
    I3DBA = (1/alb_mean[36] - 1/alb_mean[40]) / alb_mean[45]
    NDCI = ((alb_mean[40]-alb_mean[38])-(alb_mean[45]-alb_mean[38]))*((alb_mean[40]-alb_mean[38])/(alb_mean[45]-alb_mean[38]))
    MCI = (alb_mean[40]-alb_mean[36])/(alb_mean[40]+alb_mean[36])
    II = np.log(alb_mean[26])/np.log(alb_mean[56])

    print('wvls = {},{},{},{}'.format(wvl[40],wvl[36],wvl[38],wvl[45]))

    print("\nINDEX VALUES")
    print("2DBA Index: ",I2DBA)
    print("3DBA index: ", I3DBA)
    print("NDCI index: ", NDCI)
    print("MCI index: ", MCI)
    print("Impurity Index: ", II)

print(np.mean(BBA_out))