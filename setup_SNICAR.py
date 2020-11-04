"""
This function takes user-defined values from the driver script, feeds them into
a reduced form of SNICAR that outputs the optical thickness, single scattering albedo
and asymmetry parameter that then feed into disort

"""
import xarray as xr
import numpy as np

def setup_SNICAR(MIE, GO, DIRECT, DELTA, layer_type, snw_shp, shp_fctr, snw_ar, R_sfc, dz, rho_snw, rds_snw, side_length, depth,
                 nbr_lyr, mss_cnc_soot1, mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2,
                    mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1,
                     mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1,
                      mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, mss_cnc_snw_alg, mss_cnc_glacier_algae1,
                       mss_cnc_glacier_algae2):

    if MIE and GO:
        print("*** ERROR: BOTH MIE AND GO SELECTED: PLEASE SET ONE TO FALSE *** ")
   
    ##############################################
    ## 4) SET LAP CHARACTERISTICS
    else:
        nbr_aer = 16  # Define total number of different LAPs/aerosols in model

        # set filename stubs
        stb1 = 'RealPhenol_algae_geom_'  # %name stub 1
        stb2 = '.nc'  # file extension
        wrkdir2 = '/data/home/tothepoles/Desktop/bioDISORTpy/Data/Algal_Optical_Props/'  # working directory
        snw_stb1 = 'snw_alg_'  # name stub for snow algae

        # CHOOSE DIMENSIONS OF GLACIER ALGAE 1
        algae_r = 6  # algae radius
        algae_l = 120  # algae length
        glacier_algae1 = str(wrkdir2 + stb1 + str(algae_r) + '_' + str(algae_l) + stb2)  # create filename string

        # CHOOSE DIMENSIONS OF GLACIER ALGAE 2
        algae2_r = 6  # algae radius
        algae2_l = 20  # algae length
        glacier_algae2 = str(wrkdir2 + stb1 + str(algae2_r) + '_' + str(algae2_l) + stb2)  # create filename string

        # CHOOSE SNOW ALGAE DIAMETER
        snw_algae_r = 1  # snow algae diameter
        snw_alg = str(wrkdir2 + snw_stb1 + str(snw_algae_r) + stb2)  # create filename string

        # SET FILE NAMES CONTAINING OPTICAL PARAMETERS FOR ALL IMPURITIES:

        FILE_soot1  = 'mie_sot_ChC90_dns_1317.nc'
        FILE_soot2  = 'miecot_slfsot_ChC90_dns_1317.nc'
        FILE_dust1  = 'aer_dst_bln_20060904_01.nc'
        FILE_dust2  = 'aer_dst_bln_20060904_02.nc'
        FILE_dust3  = 'aer_dst_bln_20060904_03.nc'
        FILE_dust4  = 'aer_dst_bln_20060904_04.nc'
        FILE_ash1  = 'volc_ash_mtsthelens_20081011.nc'
        FILE_GRISdust1 = 'dust_greenland_Cook_CENTRAL_20190911.nc'
        FILE_GRISdust2 = 'dust_greenland_Cook_HIGH_20190911.nc'
        FILE_GRISdust3 = 'dust_greenland_Cook_LOW_20190911.nc'
        FILE_GRISdustP1 = 'dust_greenland_L_20150308.nc'
        FILE_GRISdustP2 = 'dust_greenland_C_20150308.nc'
        FILE_GRISdustP3 = 'dust_greenland_H_20150308.nc'
        FILE_snw_alg  = snw_alg # snow algae (c nivalis)
        FILE_glacier_algae1 = glacier_algae1 # Glacier algae
        FILE_glacier_algae2 = glacier_algae2 # Glacier algae

        # set working directory (location of netcdf library)
        dir_base = "/data/home/tothepoles/Desktop/bioDISORTpy/"

        if MIE:

            dir_files = "Data/Mie_files/"
            # retrieve wavelength from arbitrary choice of netcdf file
            temp = xr.open_dataset(str(dir_base + dir_files + "ice_wrn_0500.nc"))
            wvl = np.array(temp['wvl'].values)
            wvl = wvl * 1e6
            nbr_wvl = len(wvl)

        if GO:
            dir_files = "Data/GO_files/"
            # retrieve wavelength from arbitrary choice of netcdf file
            temp = xr.open_dataset(str(dir_base+dir_files+"ice_geom_5000_5000.nc"))
            wvl = np.array(temp['wvl'].values)
            wvl = wvl*1e6
            nbr_wvl = len(wvl)

        # set reflectance of underlying surface
        R_sfc = [R_sfc for _ in range(nbr_wvl)]
        R_sfc = np.array(R_sfc)

        # Incoming Irradiance
        flx_slr = []

        if DIRECT:

            with open(str(dir_base + dir_files + "mlw_sfc_flx_frc_clr.txt")) as file:
                for line in file:
                    line = float(line.rstrip("\n"))
                    flx_slr.append(line)
            flx_slr = np.array(flx_slr)
            flx_slr[flx_slr==0]=1e-30


        else:

            with open(str(dir_base + dir_files + "mlw_sfc_flx_frc_cld.txt")) as file:
                for line in file:
                    line = float(line.rstrip("\n"))
                    flx_slr.append(line)

            flx_slr = np.array(flx_slr)
            flx_slr[flx_slr==0]=1e-30

        # Read in ice optical properties
        # set string stubs for reading in ice optical files

        #set up empty arrays
        SSA_snw = np.empty([nbr_lyr, nbr_wvl])
        MAC_snw = np.empty([nbr_lyr, nbr_wvl])
        g_snw = np.empty([nbr_lyr, nbr_wvl])

        if MIE:

            fl_stb1 = "ice_wrn_"
            fl_stb2 = ".nc"

            for i in np.arange(0,nbr_lyr,1):

                if rds_snw[i] ==0:
                    print("ERROR: ICE GRAIN RADIUS SET TO ZERO")

                else:
                    s1 = "{}".format(str(rds_snw[i]))
                    s1 = s1.rjust(4,'0')
                    FILE_ice = str(dir_base+dir_files+fl_stb1+s1+fl_stb2)
                
                    # read in single scattering albedo, MAC and g for ice crystals in each layer
                    file =  xr.open_dataset(FILE_ice)


                    SSA = file['ss_alb'].values
                    SSA_snw[i,:] = SSA

                    ext_cff_mss = file['ext_cff_mss'].values
                    MAC_snw[i,:] = ext_cff_mss

                    asm_prm = file['asm_prm'].values
                    g_snw[i,:] = asm_prm

            
                    

        if GO:
            fl_stb1 = "ice_geom_"
            fl_stb2 = ".nc"

            for i in np.arange(0, nbr_lyr, 1):

                if (side_length[i] == 0) | (depth[i] == 0):

                    print("ERROR: ICE GRAIN LENGTH AND/OR DEPTH SET TO ZERO")

                else:
                    s1 = str(side_length[i])
                    s2 = str(depth[i])
                    FILE_ice = str(dir_base + dir_files + fl_stb1 + s1 + "_" + s2 + fl_stb2)


                # read in single scattering albedo, MAC and g for ice crystals in each layer
                file =  xr.open_dataset(FILE_ice)


                SSA = file['ss_alb'].values
                SSA_snw[i,:] = SSA

                ext_cff_mss = file['ext_cff_mss'].values
                MAC_snw[i,:] = ext_cff_mss

                asm_prm = file['asm_prm'].values
                g_snw[i,:] = asm_prm

       
        ###############################################################
        ## NEW GRAIN SHAPE ROUTINE
        ################################################################

        # Constants for aspherical ice particles
        # g_snw asymmetry factor parameterization coefficients (6 bands) from
        # Table 3 & Eqs. 6-7 in He et al. (2017)
        # assume same values for 4-5 um band, which leads to very small biases (<3%)
        g_wvl = np.array([0.25,0.70,1.41,1.90,2.50,3.50,4.00,5.00]) # wavelength (um) division point
        g_wvl_center = np.array(g_wvl[1:8])/2 + np.array(g_wvl[0:7])/2  # center point for wavelength band
        g_b0 = np.array([9.76029E-01,9.67798E-01,1.00111E+00,1.00224E+00,9.64295E-01,9.97475E-01,9.97475E-01])
        g_b1 = np.array([5.21042E-01,4.96181E-01,1.83711E-01,1.37082E-01,5.50598E-02,8.48743E-02,8.48743E-02])
        g_b2 = np.array([-2.66792E-04,1.14088E-03,2.37011E-04,-2.35905E-04,8.40449E-04,-4.71484E-04,-4.71484E-04])
        
        # Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007
        g_F07_c2 = np.array([1.349959e-1,1.115697e-1,9.853958e-2,5.557793e-2,-1.233493e-1,0.0,0.0])
        g_F07_c1 = np.array([-3.987320e-1,-3.723287e-1,-3.924784e-1,-3.259404e-1,4.429054e-2,-1.726586e-1,-1.726586e-1])
        g_F07_c0 = np.array([7.938904e-1,8.030084e-1,8.513932e-1,8.692241e-1,7.085850e-1,6.412701e-1,6.412701e-1])
        g_F07_p2 = np.array([3.165543e-3,2.014810e-3,1.780838e-3,6.987734e-4,-1.882932e-2,-2.277872e-2,-2.277872e-2])
        g_F07_p1 = np.array([1.140557e-1,1.143152e-1,1.143814e-1,1.071238e-1,1.353873e-1,1.914431e-1,1.914431e-1])
        g_F07_p0 = np.array([5.292852e-1,5.425909e-1,5.601598e-1,6.023407e-1,6.473899e-1,4.634944e-1,4.634944e-1])
        fs_hex = 0.788 # shape factor for hexagonal plate (reference)
        

        for i in np.arange(0,nbr_lyr,1):
            
            if layer_type[i] == 0:
                
                if snw_shp[i] == 0: # snow

                    pass # if layer type is for spheres, no changes required

                
                elif snw_shp[i] == 1: # 1 = spheroid, He et al. (2017) parameterization
                    
                    diam_ice = 2.0 * rds_snw[i] # effective snow grain diameter
                    
                    if shp_fctr[i] == 0:

                        fs_sphd = 0.929 # default shape factor for spheroid; He et al. (2017), Table 1
                    
                    else:
                        
                        fs_sphd = shp_fctr[i] # if shp_factor not 0, then use user-defined value
                    
                    
                    if snw_ar[i] == 0:

                        AR_tmp = 0.5 # default aspect ratio for spheroid; He et al. (2017), Table 1
                    
                    else:

                        AR_tmp = snw_ar[i]
                    
                    g_snw_Cg_tmp = g_b0 * (fs_sphd/fs_hex)**g_b1 * diam_ice**g_b2 # Eq.7, He et al. (2017)
                    gg_snw_F07_tmp = g_F07_c0 + g_F07_c1 * AR_tmp + g_F07_c2 * AR_tmp**2 # Eqn. 3.1 in Fu (2007)


                elif snw_shp[i] == 2: # 3=hexagonal plate, He et al. 2017 parameterization

                    diam_ice = 2.0 * rds_snw[i] # effective snow grain diameter
                    
                    if shp_fctr[i] == 0:
                        
                        fs_hex0 = 0.788 # default shape factor for hexagonal plates; He et al. (2017), Table 1
                    
                    else:

                        fs_hex0 = shp_fctr[i]
                    

                    if snw_ar[i] == 0:

                        AR_tmp = 2.5 # default aspect ratio for hexagonal plate; He et al. (2017), Table 1
                    
                    else:
                    
                        AR_tmp = snw_ar[i]
                            
                    g_snw_Cg_tmp = g_b0 * (fs_hex0/fs_hex)**g_b1 * diam_ice**g_b2 # Eq.7, He et al. (2017)
                    gg_snw_F07_tmp = g_F07_p0 + g_F07_p1 * np.log(AR_tmp) + g_F07_p2 * (np.log(AR_tmp))**2   # Eqn. 3.3 in Fu (2007)
                    

                elif snw_shp[i] == 3: # 4=koch snowflake, He et al. (2017) parameterization

                    diam_ice = 2.0 * rds_snw[i] / 0.544 # effective snow grain diameter
                    
                    if shp_fctr[i] == 0:
                    
                        fs_koch = 0.712 # default shape factor for koch snowflake; He et al. (2017), Table 1
                    
                    else:
                        
                        fs_koch = shp_fctr[i]
                    

                    if snw_ar[i] == 0:

                        AR_tmp = 2.5 # default aspect ratio for koch snowflake; He et al. (2017), Table 1
                    
                    else:

                        AR_tmp = snw_ar[i]
                    
                    g_snw_Cg_tmp = g_b0 * (fs_koch/fs_hex)**g_b1 * diam_ice**g_b2 # Eq.7, He et al. (2017)
                    gg_snw_F07_tmp = g_F07_p0 + g_F07_p1 * np.log(AR_tmp) + g_F07_p2 * (np.log(AR_tmp))**2  # Eqn. 3.3 in Fu (2007)
                    
                
                if snw_shp[i] > 0:

                    from scipy.interpolate import pchip
                    # 6 wavelength bands for g_snw to be interpolated into 480-bands of SNICAR
                    # shape-preserving piecewise interpolation into 480-bands
                    g_Cg_intp = pchip(g_wvl_center,g_snw_Cg_tmp)(wvl)
                    gg_F07_intp = pchip(g_wvl_center,gg_snw_F07_tmp)(wvl)
                    g_snw_F07 = gg_F07_intp + (1.0 - gg_F07_intp) / SSA_snw[i,:] / 2 # Eq.2.2 in Fu (2007)
                    g_snw[i,:] = g_snw_F07 * g_Cg_intp # Eq.6, He et al. (2017)
                    g_snw[i,381:480] = g_snw[i,380] # assume same values for 4-5 um band, with very small biases (<3%)
                
                
                g_snw[g_snw > 0.99] = 0.99 # avoid unreasonable values (so far only occur in large-size spheroid cases)


        # open netcdf files
        FILE_soot1 = xr.open_dataset(str(dir_base + dir_files+ FILE_soot1))
        FILE_soot2 = xr.open_dataset(str(dir_base + dir_files+ FILE_soot2))
        FILE_dust1 = xr.open_dataset(str(dir_base + dir_files+ FILE_dust1))
        FILE_dust2 = xr.open_dataset(str(dir_base + dir_files+ FILE_dust2))
        FILE_dust3 = xr.open_dataset(str(dir_base + dir_files+ FILE_dust3))
        FILE_dust4 = xr.open_dataset(str(dir_base + dir_files+ FILE_dust4))
        FILE_ash1 = xr.open_dataset(str(dir_base + dir_files+ FILE_ash1))
        FILE_GRISdust1 = xr.open_dataset(str(dir_base + dir_files+ FILE_GRISdust1))
        FILE_GRISdust2 = xr.open_dataset(str(dir_base + dir_files+ FILE_GRISdust2))
        FILE_GRISdust3 = xr.open_dataset(str(dir_base + dir_files+ FILE_GRISdust3))
        FILE_GRISdustP1 = xr.open_dataset(str(dir_base + dir_files+ FILE_GRISdustP1))
        FILE_GRISdustP2 = xr.open_dataset(str(dir_base + dir_files+ FILE_GRISdustP2))
        FILE_GRISdustP3 = xr.open_dataset(str(dir_base + dir_files+ FILE_GRISdustP3))
        FILE_snw_alg = xr.open_dataset(FILE_snw_alg)
        FILE_glacier_algae1 = xr.open_dataset(FILE_glacier_algae1)
        FILE_glacier_algae2 = xr.open_dataset(FILE_glacier_algae2)

        # read in aerosol optical properties
        SSAaer = np.zeros([nbr_aer,nbr_wvl])

        SSAaer[0,:] = FILE_soot1['ss_alb'].values
        SSAaer[1,:] = FILE_soot2['ss_alb'].values
        SSAaer[2,:] = FILE_dust1['ss_alb'].values
        SSAaer[3,:] = FILE_dust2['ss_alb'].values
        SSAaer[4,:] = FILE_dust3['ss_alb'].values
        SSAaer[5,:] = FILE_dust4['ss_alb'].values
        SSAaer[6,:] = FILE_ash1['ss_alb'].values
        SSAaer[7,:] = FILE_GRISdust1['ss_alb'].values
        SSAaer[8,:] = FILE_GRISdust2['ss_alb'].values
        SSAaer[9,:] = FILE_GRISdust3['ss_alb'].values
        SSAaer[10,:] = FILE_GRISdustP1['ss_alb'].values
        SSAaer[11,:] = FILE_GRISdustP2['ss_alb'].values
        SSAaer[12,:] = FILE_GRISdustP3['ss_alb'].values
        SSAaer[13,:] = FILE_snw_alg['ss_alb'].values
        SSAaer[14,:] = FILE_glacier_algae1['ss_alb'].values
        SSAaer[15,:] = FILE_glacier_algae2['ss_alb'].values

        MACaer = np.zeros([nbr_aer, nbr_wvl])

        MACaer[0,:] = FILE_soot1['ext_cff_mss'].values
        MACaer[1,:] = FILE_soot2['ext_cff_mss'].values
        MACaer[2,:] = FILE_dust1['ext_cff_mss'].values
        MACaer[3,:] = FILE_dust2['ext_cff_mss'].values
        MACaer[4,:] = FILE_dust3['ext_cff_mss'].values
        MACaer[5,:] = FILE_dust4['ext_cff_mss'].values
        MACaer[6,:] = FILE_ash1['ext_cff_mss'].values
        MACaer[7,:] = FILE_GRISdust1['ext_cff_mss'].values
        MACaer[8,:] = FILE_GRISdust2['ext_cff_mss'].values
        MACaer[9,:] = FILE_GRISdust3['ext_cff_mss'].values
        MACaer[10,:] = FILE_GRISdustP1['ext_cff_mss'].values
        MACaer[11,:] = FILE_GRISdustP2['ext_cff_mss'].values
        MACaer[12,:] = FILE_GRISdustP3['ext_cff_mss'].values
        MACaer[13,:] = FILE_snw_alg['ext_cff_mss'].values
        MACaer[14,:] = FILE_glacier_algae1['ext_cff_mss'].values
        MACaer[15,:] = FILE_glacier_algae2['ext_cff_mss'].values

        Gaer = np.zeros([nbr_aer,nbr_wvl])

        Gaer[0,:] = FILE_soot1['asm_prm'].values
        Gaer[1,:] = FILE_soot2['asm_prm'].values
        Gaer[2,:] = FILE_dust1['asm_prm'].values
        Gaer[3,:] = FILE_dust2['asm_prm'].values
        Gaer[4,:] = FILE_dust3['asm_prm'].values
        Gaer[5,:] = FILE_dust4['asm_prm'].values
        Gaer[6,:] = FILE_ash1['asm_prm'].values
        Gaer[7,:] = FILE_GRISdust1['asm_prm'].values
        Gaer[8,:] = FILE_GRISdust2['asm_prm'].values
        Gaer[9,:] = FILE_GRISdust3['asm_prm'].values
        Gaer[10,:] = FILE_GRISdustP1['asm_prm'].values
        Gaer[11,:] = FILE_GRISdustP2['asm_prm'].values
        Gaer[12,:] = FILE_GRISdustP3['asm_prm'].values
        Gaer[13,:] = FILE_snw_alg['asm_prm'].values
        Gaer[14,:] = FILE_glacier_algae1['asm_prm'].values
        Gaer[15,:] = FILE_glacier_algae2['asm_prm'].values


        # load mass concentrations per layer into numpy array (one row per layer, one column per umpurity)
        # and convert to kg/kg unit

        MSSaer = np.zeros([nbr_lyr, nbr_aer])
        MSSaer[0:nbr_lyr,0] = mss_cnc_soot1
        MSSaer[0:nbr_lyr,1] = mss_cnc_soot2
        MSSaer[0:nbr_lyr,2] = mss_cnc_dust1
        MSSaer[0:nbr_lyr,3] = mss_cnc_dust2
        MSSaer[0:nbr_lyr,4] = mss_cnc_dust3
        MSSaer[0:nbr_lyr,5] = mss_cnc_dust4
        MSSaer[0:nbr_lyr,6] = mss_cnc_ash1
        MSSaer[0:nbr_lyr,7] = mss_cnc_GRISdust1
        MSSaer[0:nbr_lyr,8] = mss_cnc_GRISdust2
        MSSaer[0:nbr_lyr,9] = mss_cnc_GRISdust3
        MSSaer[0:nbr_lyr,10] = mss_cnc_GRISdustP1
        MSSaer[0:nbr_lyr,11] = mss_cnc_GRISdustP2
        MSSaer[0:nbr_lyr,12] = mss_cnc_GRISdustP3
        MSSaer[0:nbr_lyr,13] = mss_cnc_snw_alg
        MSSaer[0:nbr_lyr,14] = mss_cnc_glacier_algae1
        MSSaer[0:nbr_lyr,15] = mss_cnc_glacier_algae2

        MSSaer = MSSaer*1e-9


        #####################################
        # Begin solving Radiative Transfer
        #####################################

        """
        #1. Calculate effective tau (optical depth), SSA (single scattering albedo) and 
        # g (assymetry parameter) for the ice + impurities mixture.
    
        # SSA and g for the individual components has been calculated using Mie theory and
        # stored in a netcdf file. Here, these values are combined to give an overall
        # SSA and g for the ice + impurity mixture
    
        """

        # initialize arrays
        g_sum = np.zeros([nbr_lyr, nbr_wvl])
        tau = np.zeros([nbr_lyr, nbr_wvl])
        SSA = np.zeros([nbr_lyr, nbr_wvl])
        SSA_star = np.zeros([nbr_lyr, nbr_wvl])
        g = np.zeros([nbr_lyr, nbr_wvl])
        g_star = np.zeros([nbr_lyr, nbr_wvl])
        L_aer = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
        tau_aer = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
        tau_star = np.zeros([nbr_lyr, nbr_wvl])
        tau_sum = np.zeros([nbr_lyr, nbr_wvl])
        SSA_sum = np.zeros([nbr_lyr, nbr_wvl])
        L_snw = np.zeros(nbr_lyr)
        tau_snw = np.zeros([nbr_lyr,nbr_wvl])

        # for each layer, the layer mass (L) is density * layer thickness
        # for each layer the optical depth is the layer mass * the mass extinction coefficient
        # first for the ice in each layer
        
        for i in range(nbr_lyr):

            L_snw[i] = rho_snw[i] * dz[i]
            
            tau_snw[i, :] = L_snw[i] * MAC_snw[i, :]

        # then for the LAPs in each layer
        for i in range(nbr_lyr):

            for j in range(nbr_aer):
                
                L_aer[i, j, :] = L_snw[i]*MSSaer[i, j]

                tau_aer[i, j, :] = L_aer[i, j, :]*MACaer[j, :]

                tau_sum = tau_sum + tau_aer[i, j, :]
                SSA_sum = SSA_sum + (tau_aer[i, j, :] * SSAaer[j, :])
                g_sum = g_sum + (tau_aer[i, j, :] * SSAaer[j, :] * Gaer[j, :])


        # finally, for each layer calculate the effective SSA, tau and g for the snow+LAP
        for i in range(nbr_lyr):

            tau[i,:] = tau_sum[i,:] + tau_snw[i,:]
            SSA[i,:] = (1/tau[i,:]) * (SSA_sum[i,:] + SSA_snw[i,:] * tau_snw[i,:])
            g[i, :] = (1 / (tau[i, :] * (SSA[i, :]))) * (g_sum[i,:] + (g_snw[i, :] * SSA_snw[i, :] * tau_snw[i, :]))

        ############################################
        # PERFORM DELTA TRANSFORMATION IF REQUIRED
        ############################################
        # The star represents the delta transformed quantity
        # if no delta transformation is applied, the starred quantity
        # is equal to the unstarred quantity

        if DELTA:
            g_star = g/(1+g)
            SSA_star = ((1-(g**2))*SSA)/(1-(SSA*(g**2)))
            tau_star = (1-(SSA*(g**2)))*tau

        else:
            g_star = g
            SSA_star = SSA
            tau_star = tau



    return flx_slr, g_star, SSA_star, tau_star, wvl