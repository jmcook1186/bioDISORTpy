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


######################################
## 3. RADIATIVE TRANSFER CONFIGURATION
#######################################

DIRECT = 1        # 1= Direct-beam incident flux, 0= Diffuse incident flux
DELTA = 1        # 1= Apply Delta approximation, 0= No delta
MIE = True          # use single scattering optical properties generating using Mie scattering
GO = False          # use single scattering optical properties generating using geometric optics
coszen   = 0.57      # if DIRECT give cosine of solar zenith angle 
save_path = '/home/tothepoles/Desktop/bioDISORTpy/' # path to save figures to
print_band_ratios = False


#############################################
## 4. SET PHYSICAL PROPERTIES OF THE ICE/SNOW
#############################################

# grain shapes are only active in Mie scattering mode - GO assumes hexagonal columns.
# snw_shp can be 0 = sphere, 1 = spheroid, 2 = hexagonal plate, 3= koch snowflake
# shp_fctr = ratio of nonspherical grain effective radii to that of equal-volume sphere
    # 0=use recommended default value (He et al. 2017);
    # use user-specified value (between 0 and 1)
    # only activated when sno_shp > 1 (i.e. nonspherical)

dz = [0.001, 0.05, 0.2, 0.2, 0.5] # thickness of each vertical layer (unit = m)
nbr_lyr = len(dz)  # number of snow layers
layer_type = [0,0,1,1,1]
R_sfc = 0.15 # reflectance of undrlying surface - set across all wavelengths
rho_snw = [600, 800, 900, 900, 900] # density of each layer (unit = kg m-3)
rds_snw = [2000,2000,2000,2000,2000] # effective grain radius of snow/bubbly ice
rwater = [0, 0, 0, 0, 0] # if  using Mie calculations, add radius of optional liquid water coating
snw_shp =[2,2,2,2,2] # grain shape(He et al. 2016, 2017)
shp_fctr = [0,0,0,0,0] # shape factor (ratio of aspherical grain radii to that of equal-volume sphere)
snw_ar = [0.9,0.9,0.9,0.9,0.9] # aspect ratio (ratio of width to length)

# if using GeometricOptics, set side_length and depth
side_length = [30000,30000,30000,30000,30000] 
depth = [30000,30000,30000,30000,30000]

# SET IMPURITY MASS CONCENTRATIONS IN EACH LAYER
# units are ppb or ng/g i.e. 1e3 = 1 ppm or 1 ug/g, 1e6 = 1 ppt or 1 mg/g

mss_cnc_soot1 = [0,0,0,0,0]    # uncoated black carbon
mss_cnc_soot2 = [0,0,0,0,0]    # coated black carbon
mss_cnc_dust1 = [0,0,0,0,0]    # global average dust 1
mss_cnc_dust2 = [0,0,0,0,0]    # global average dust 2
mss_cnc_dust3 = [0,0,0,0,0]    # global average dust 3
mss_cnc_dust4 = [0,0,0,0,0]    # global average dust 4
mss_cnc_ash1 = [0,0,0,0,1]    # volcanic ash species 1
mss_cnc_GRISdust1 = [0,0,0,0,0]    # GRIS dust 1 (Cook et al. 2019 "mean")
mss_cnc_GRISdust2 = [0,0,0,0,0]    # GRIS dust 2 (Cook et al. 2019 HIGH)
mss_cnc_GRISdust3 = [0,0,0,0,0]    # GRIS dust 3 (Cook et al. 2019 LOW)
mss_cnc_GRISdustP1 = [0,0,0,0,0]  # GRIS dust 1 (Polashenki2015: low hematite)
mss_cnc_GRISdustP2 = [0,0,0,0,0]  # GRIS dust 1 (Polashenki2015: median hematite)
mss_cnc_GRISdustP3 = [0,0,0,0,0]  # GRIS dust 1 (Polashenki2015: median hematite)
mss_cnc_snw_alg = [0,0,0,0,0]# Snow Algae (spherical, C nivalis)
mss_cnc_glacier_algae1 = [0,0,0,0,0]    # glacier algae type1
mss_cnc_glacier_algae2 = [0,0,0,0,0]    # glacier algae type2

# DISORT CONFIG
prnt = np.array([True, True, True, True, True]) # determines what info to print to console
umu = [0.2,0.4,0.6,0.8,1.0]  # cosine of viewing zenith angle
coszen = [0.57] #[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] # cosine of solar zenith angle
phi0 = 0  # solar azimuth angle
phi = 0 # viewing azimuth angle
albedo = R_sfc  # albedo of underlying surface
uTau = 0.0  # optical thickness where fluxes are calculated
Nstr = 16 # number of streams to include in model

# RETRIEVE OPTICAL THICKNESS, SSA, ASYMMETRY PARAMETER FROM SNICAR

[flx_slr, g_star, SSA_star, tau_star, wvl] = setup_SNICAR.setup_SNICAR(
    MIE, GO, DIRECT, DELTA, layer_type, snw_shp, shp_fctr, snw_ar, R_sfc, dz, rho_snw, rds_snw, side_length, depth, nbr_lyr,
    mss_cnc_soot1, mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2,
    mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1,
    mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1,
    mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, mss_cnc_snw_alg, mss_cnc_glacier_algae1,
    mss_cnc_glacier_algae2)

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

        dTau = tau_star[:,i]
        w0 = SSA_star[:,i]
        iphas = np.ones(len(dz),dtype='int')*3
        iphas[:] = 3
        gg = g_star[:,i]
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
    plt.plot(wvl[0:250],alb_out[i,0:250],label='coszen: {}'.format(str(coszen[i])))

plt.plot(wvl[0:250],alb_mean[0:250],linestyle='dashed',color='k',label='zenMean')
plt.ylabel('Albedo')
plt.xlabel('Wavelength')
plt.xlim(0.35,2.5)
plt.legend()
plt.savefig(str(save_path+'albedo_DISORT.png'))

plt.figure()
plt.plot(coszen,BBA_out)
plt.xlabel('Cosine of solar zenith')
plt.ylabel('BBA')
plt.ylim(0,1)
plt.savefig(str(save_path+'BBA_DISORT.png'))

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