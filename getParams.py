import numpy as np
import xarray as xr

"""
takes input values from GenerateTrainingData.py and returns spectral asymmetry parameter, single scattering albedo and
optical thickness in each layer to feed into DISORT.

"""
def getParams_mie(DIRECT, R_sfc, dz, rho_snw, rds_snw, nbr_lyr, nbr_aer,
                 mss_cnc_soot1, mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2,
                 mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1,
                 mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1,
                 mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, mss_cnc_snw_alg, mss_cnc_glacier_algae1,
                 mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2, FILE_dust3, FILE_dust4,
                 FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2,
                 FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):


    # set working directory (location of netcdf library)
    dir_base = "/home/joe/Code/Albedo_NeuralNet/"
    dir_alg = "Data/Algal_Optical_Props/"
    dir_mie_files = "Data/Mie_files/"

    # retrieve wavelength from arbitrary choice of netcdf file
    temp = xr.open_dataset(str(dir_base + dir_mie_files + "ice_wrn_0500.nc"))
    wvl = np.array(temp['wvl'].values)
    wvl = wvl * 1e6
    nbr_wvl = len(wvl)

    # set reflectance of underlying surface
    R_sfc = [R_sfc for _ in range(nbr_wvl)]
    R_sfc = np.array(R_sfc)


    flx_slr = []

    if DIRECT:

        with open(str(dir_base + dir_mie_files + "mlw_sfc_flx_frc_clr.txt")) as file:
            for line in file:
                line = float(line.rstrip("\n"))
                flx_slr.append(line)
        flx_slr = np.array(flx_slr)
        flx_slr[flx_slr == 0] = 1e-30


    else:

        with open(str(dir_base + dir_mie_files + "mlw_sfc_flx_frc_cld.txt")) as file:
            for line in file:
                line = float(line.rstrip("\n"))
                flx_slr.append(line)

        flx_slr = np.array(flx_slr)
        flx_slr[flx_slr == 0] = 1e-30


    # Read in ice optical properties
    # set string stubs for reading in ice optical files

    fl_stb1 = "ice_wrn_"
    fl_stb2 = ".nc"

    # set up empty arrays
    SSA_snw = np.empty([nbr_lyr, nbr_wvl])
    MAC_snw = np.empty([nbr_lyr, nbr_wvl])
    g_snw = np.empty([nbr_lyr, nbr_wvl])

    for i in np.arange(0, nbr_lyr, 1):

        if rds_snw[i] == 0:
            print("ERROR: ICE GRAIN RADIUS SET TO ZERO")

        else:
            s1 = "{}".format(str(rds_snw[i]))
            s1 = s1.rjust(4, '0')
            FILE_ice = str(dir_base + dir_mie_files + fl_stb1 + s1 + fl_stb2)

        # read in single scattering albedo, MAC and g for ice crystals in each layer
        with xr.open_dataset(FILE_ice) as temp:
            SSA = temp['ss_alb'].values
            SSA_snw[i, :] = SSA

            ext_cff_mss = temp['ext_cff_mss'].values
            MAC_snw[i, :] = ext_cff_mss

            asm_prm = temp['asm_prm'].values
            g_snw[i, :] = asm_prm

    # open netcdf files
    FILE_soot1 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_soot1))
    FILE_soot2 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_soot2))
    FILE_dust1 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_dust1))
    FILE_dust2 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_dust2))
    FILE_dust3 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_dust3))
    FILE_dust4 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_dust4))
    FILE_ash1 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_ash1))
    FILE_GRISdust1 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_GRISdust1))
    FILE_GRISdust2 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_GRISdust2))
    FILE_GRISdust3 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_GRISdust3))
    FILE_GRISdustP1 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_GRISdustP1))
    FILE_GRISdustP2 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_GRISdustP2))
    FILE_GRISdustP3 = xr.open_dataset(str(dir_base + dir_mie_files + FILE_GRISdustP3))
    FILE_snw_alg = xr.open_dataset(FILE_snw_alg)
    FILE_glacier_algae1 = xr.open_dataset(FILE_glacier_algae1)
    FILE_glacier_algae2 = xr.open_dataset(FILE_glacier_algae2)

    # read in aerosol optical properties
    SSAaer = np.zeros([nbr_aer, nbr_wvl])

    SSAaer[0, :] = FILE_soot1['ss_alb'].values
    SSAaer[1, :] = FILE_soot2['ss_alb'].values
    SSAaer[2, :] = FILE_dust1['ss_alb'].values
    SSAaer[3, :] = FILE_dust2['ss_alb'].values
    SSAaer[4, :] = FILE_dust3['ss_alb'].values
    SSAaer[5, :] = FILE_dust4['ss_alb'].values
    SSAaer[6, :] = FILE_ash1['ss_alb'].values
    SSAaer[7, :] = FILE_GRISdust1['ss_alb'].values
    SSAaer[8, :] = FILE_GRISdust2['ss_alb'].values
    SSAaer[9, :] = FILE_GRISdust3['ss_alb'].values
    SSAaer[10, :] = FILE_GRISdustP1['ss_alb'].values
    SSAaer[11, :] = FILE_GRISdustP2['ss_alb'].values
    SSAaer[12, :] = FILE_GRISdustP3['ss_alb'].values
    SSAaer[13, :] = FILE_snw_alg['ss_alb'].values
    SSAaer[14, :] = FILE_glacier_algae1['ss_alb'].values
    SSAaer[15, :] = FILE_glacier_algae2['ss_alb'].values

    MACaer = np.zeros([nbr_aer, nbr_wvl])

    MACaer[0, :] = FILE_soot1['ext_cff_mss'].values
    MACaer[1, :] = FILE_soot2['ext_cff_mss'].values
    MACaer[2, :] = FILE_dust1['ext_cff_mss'].values
    MACaer[3, :] = FILE_dust2['ext_cff_mss'].values
    MACaer[4, :] = FILE_dust3['ext_cff_mss'].values
    MACaer[5, :] = FILE_dust4['ext_cff_mss'].values
    MACaer[6, :] = FILE_ash1['ext_cff_mss'].values
    MACaer[7, :] = FILE_GRISdust1['ext_cff_mss'].values
    MACaer[8, :] = FILE_GRISdust2['ext_cff_mss'].values
    MACaer[9, :] = FILE_GRISdust3['ext_cff_mss'].values
    MACaer[10, :] = FILE_GRISdustP1['ext_cff_mss'].values
    MACaer[11, :] = FILE_GRISdustP2['ext_cff_mss'].values
    MACaer[12, :] = FILE_GRISdustP3['ext_cff_mss'].values
    MACaer[13, :] = FILE_snw_alg['ext_cff_mss'].values
    MACaer[14, :] = FILE_glacier_algae1['ext_cff_mss'].values
    MACaer[15, :] = FILE_glacier_algae2['ext_cff_mss'].values

    Gaer = np.zeros([nbr_aer, nbr_wvl])

    Gaer[0, :] = FILE_soot1['asm_prm'].values
    Gaer[1, :] = FILE_soot2['asm_prm'].values
    Gaer[2, :] = FILE_dust1['asm_prm'].values
    Gaer[3, :] = FILE_dust2['asm_prm'].values
    Gaer[4, :] = FILE_dust3['asm_prm'].values
    Gaer[5, :] = FILE_dust4['asm_prm'].values
    Gaer[6, :] = FILE_ash1['asm_prm'].values
    Gaer[7, :] = FILE_GRISdust1['asm_prm'].values
    Gaer[8, :] = FILE_GRISdust2['asm_prm'].values
    Gaer[9, :] = FILE_GRISdust3['asm_prm'].values
    Gaer[10, :] = FILE_GRISdustP1['asm_prm'].values
    Gaer[11, :] = FILE_GRISdustP2['asm_prm'].values
    Gaer[12, :] = FILE_GRISdustP3['asm_prm'].values
    Gaer[13, :] = FILE_snw_alg['asm_prm'].values
    Gaer[14, :] = FILE_glacier_algae1['asm_prm'].values
    Gaer[15, :] = FILE_glacier_algae2['asm_prm'].values

    # load mass concentrations per layer into numpy array (one row per layer, one column per umpurity)
    # and convert to kg/kg unit

    MSSaer = np.zeros([nbr_lyr, nbr_aer])
    MSSaer[0:nbr_lyr, 0] = mss_cnc_soot1
    MSSaer[0:nbr_lyr, 1] = mss_cnc_soot2
    MSSaer[0:nbr_lyr, 2] = mss_cnc_dust1
    MSSaer[0:nbr_lyr, 3] = mss_cnc_dust2
    MSSaer[0:nbr_lyr, 4] = mss_cnc_dust3
    MSSaer[0:nbr_lyr, 5] = mss_cnc_dust4
    MSSaer[0:nbr_lyr, 6] = mss_cnc_ash1
    MSSaer[0:nbr_lyr, 7] = mss_cnc_GRISdust1
    MSSaer[0:nbr_lyr, 8] = mss_cnc_GRISdust2
    MSSaer[0:nbr_lyr, 9] = mss_cnc_GRISdust3
    MSSaer[0:nbr_lyr, 10] = mss_cnc_GRISdustP1
    MSSaer[0:nbr_lyr, 11] = mss_cnc_GRISdustP2
    MSSaer[0:nbr_lyr, 12] = mss_cnc_GRISdustP3
    MSSaer[0:nbr_lyr, 13] = mss_cnc_snw_alg
    MSSaer[0:nbr_lyr, 14] = mss_cnc_glacier_algae1
    MSSaer[0:nbr_lyr, 15] = mss_cnc_glacier_algae2

    MSSaer = MSSaer * 1e-9

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
    SSA_sum = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau = np.zeros([nbr_lyr, nbr_wvl])
    SSA = np.zeros([nbr_lyr, nbr_wvl])
    g = np.zeros([nbr_lyr, nbr_wvl])
    L_aer = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau_aer = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau_sum = np.zeros([nbr_lyr, nbr_wvl])
    SSA_sum = np.zeros([nbr_lyr, nbr_wvl])
    L_snw = np.zeros(nbr_lyr)
    tau_snw = np.zeros([nbr_lyr, nbr_wvl])
    direct = np.zeros([nbr_lyr, nbr_wvl])
    F_net = np.zeros([nbr_lyr, nbr_wvl])
    F_btm_net = np.zeros([1, nbr_wvl])
    F_top_net = np.zeros([1, nbr_wvl])
    intensity = np.zeros([nbr_lyr, nbr_wvl])
    F_top_pls = np.zeros([1, nbr_wvl])
    F_up = np.zeros([nbr_lyr, nbr_wvl])
    F_down = np.zeros([nbr_lyr, nbr_wvl])
    F_net2 = np.zeros([nbr_lyr, nbr_wvl])
    intensity2 = np.zeros([nbr_lyr, nbr_wvl])
    intensity2_top = np.zeros(nbr_wvl)
    F_abs = np.zeros([nbr_lyr, nbr_wvl])
    abs_vis = np.zeros(nbr_lyr)
    abs_nir = np.zeros(nbr_lyr)

    # for each layer, the layer mass (L) is density * layer thickness
    # for each layer the optical depth is the layer mass * the mass extinction coefficient
    # first for the ice in each layer
    for i in range(nbr_lyr):
        L_snw[i] = rho_snw[i] * dz[i]
        tau_snw[i, :] = L_snw[i] * MAC_snw[i, :]

    # then for the LAPs in each layer
    for i in range(nbr_lyr):
        for j in range(nbr_aer):
            L_aer[i, j, :] = np.multiply(L_snw[i], MSSaer[i, j])
            tau_aer[i, j, :] = np.multiply(L_aer[i, j, :], MACaer[j, :])

            tau_sum = tau_sum + tau_aer[i, j, :]
            SSA_sum = SSA_sum + (tau_aer[i, j, :] * SSAaer[j, :])
            g_sum = g_sum + (tau_aer[i, j, :] * SSAaer[j, :] * Gaer[j, :])

    # finally, for each layer calculate the effective SSA, tau and g for the snow+LAP
    for i in range(nbr_lyr):
        tau[i, :] = tau_sum[i, :] + tau_snw[i, :]
        SSA[i, :] = (1 / tau[i, :]) * (SSA_sum[i, :] + SSA_snw[i, :] * tau_snw[i, :])
        g[i, :] = (1 / (tau[i, :] * (SSA[i, :]))) * (g_sum[i, :] + (g_snw[i, :] * SSA_snw[i, :] * tau_snw[i, :]))

    ############################################
    # PERFORM DELTA TRANSFORMATION IF REQUIRED
    ############################################
    # The star represents the delta transformed quantity
    # if no delat transformation is applied, the starred quantity
    # is equal to the unstarred quantity

    g_star = g
    SSA_star = SSA
    tau_star = tau


    return flx_slr, g_star, SSA_star, tau_star, wvl


def getParams_GO(DIRECT, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer,
                mss_cnc_soot1, mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2,
                mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1,
                mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1,
                mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, mss_cnc_snw_alg, mss_cnc_glacier_algae1,
                mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2, FILE_dust3, FILE_dust4,
                FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2,
                FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):


    # set working directory (location of netcdf library)
    dir_base = "/home/joe/Code/Albedo_NeuralNet/"
    dir_alg = "Data/Algal_Optical_Props/"
    dir_GO_files = "Data/GO_files/"

    # retrieve wavelength from arbitrary choice of netcdf file
    temp = xr.open_dataset(str(dir_base + dir_GO_files + "ice_geom_5000_5000.nc"))
    wvl = np.array(temp['wvl'].values)
    wvl = wvl * 1e6
    nbr_wvl = len(wvl)

    # set reflectance of underlying surface
    R_sfc = [R_sfc for _ in range(nbr_wvl)]
    R_sfc = np.array(R_sfc)

    # Incoming Irradiance
    flx_slr = []

    if DIRECT:

        with open(str(dir_base + dir_GO_files + "mlw_sfc_flx_frc_clr.txt")) as file:
            for line in file:
                line = float(line.rstrip("\n"))
                flx_slr.append(line)
        flx_slr = np.array(flx_slr)
        flx_slr[flx_slr == 0] = 1e-30


    else:

        with open(str(dir_base + dir_GO_files + "mlw_sfc_flx_frc_cld.txt")) as file:
            for line in file:
                line = float(line.rstrip("\n"))
                flx_slr.append(line)

        flx_slr = np.array(flx_slr)
        flx_slr[flx_slr == 0] = 1e-30


    # Read in ice optical properties
    # set string stubs for reading in ice optical files

    fl_stb1 = "ice_geom_"
    fl_stb2 = ".nc"

    # set up empty arrays
    SSA_snw = np.empty([nbr_lyr, nbr_wvl])
    MAC_snw = np.empty([nbr_lyr, nbr_wvl])
    g_snw = np.empty([nbr_lyr, nbr_wvl])

    for i in np.arange(0, nbr_lyr, 1):

        if (side_length[i] == 0) | (depth[i] == 0):

            print("ERROR: ICE GRAIN LENGTH AND/OR DEPTH SET TO ZERO")

        else:

            s1 = str(side_length[i])
            s2 = str(depth[i])
            FILE_ice = str(dir_base + dir_GO_files + fl_stb1 + s1 + "_" + s2 + fl_stb2)

        # read in single scattering albedo, MAC and g for ice crystals in each layer
        with xr.open_dataset(FILE_ice) as temp:
            SSA = temp['ss_alb'].values
            SSA_snw[i, :] = SSA

            ext_cff_mss = temp['ext_cff_mss'].values
            MAC_snw[i, :] = ext_cff_mss

            asm_prm = temp['asm_prm'].values
            g_snw[i, :] = asm_prm

    # open netcdf files
    FILE_soot1 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_soot1))
    FILE_soot2 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_soot2))
    FILE_dust1 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_dust1))
    FILE_dust2 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_dust2))
    FILE_dust3 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_dust3))
    FILE_dust4 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_dust4))
    FILE_ash1 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_ash1))
    FILE_GRISdust1 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_GRISdust1))
    FILE_GRISdust2 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_GRISdust2))
    FILE_GRISdust3 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_GRISdust3))
    FILE_GRISdustP1 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_GRISdustP1))
    FILE_GRISdustP2 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_GRISdustP2))
    FILE_GRISdustP3 = xr.open_dataset(str(dir_base + dir_GO_files + FILE_GRISdustP3))
    FILE_snw_alg = xr.open_dataset(FILE_snw_alg)
    FILE_glacier_algae1 = xr.open_dataset(FILE_glacier_algae1)
    FILE_glacier_algae2 = xr.open_dataset(FILE_glacier_algae2)

    # read in aerosol optical properties
    SSAaer = np.zeros([nbr_aer, nbr_wvl])

    SSAaer[0, :] = FILE_soot1['ss_alb'].values
    SSAaer[1, :] = FILE_soot2['ss_alb'].values
    SSAaer[2, :] = FILE_dust1['ss_alb'].values
    SSAaer[3, :] = FILE_dust2['ss_alb'].values
    SSAaer[4, :] = FILE_dust3['ss_alb'].values
    SSAaer[5, :] = FILE_dust4['ss_alb'].values
    SSAaer[6, :] = FILE_ash1['ss_alb'].values
    SSAaer[7, :] = FILE_GRISdust1['ss_alb'].values
    SSAaer[8, :] = FILE_GRISdust2['ss_alb'].values
    SSAaer[9, :] = FILE_GRISdust3['ss_alb'].values
    SSAaer[10, :] = FILE_GRISdustP1['ss_alb'].values
    SSAaer[11, :] = FILE_GRISdustP2['ss_alb'].values
    SSAaer[12, :] = FILE_GRISdustP3['ss_alb'].values
    SSAaer[13, :] = FILE_snw_alg['ss_alb'].values
    SSAaer[14, :] = FILE_glacier_algae1['ss_alb'].values
    SSAaer[15, :] = FILE_glacier_algae2['ss_alb'].values

    MACaer = np.zeros([nbr_aer, nbr_wvl])

    MACaer[0, :] = FILE_soot1['ext_cff_mss'].values
    MACaer[1, :] = FILE_soot2['ext_cff_mss'].values
    MACaer[2, :] = FILE_dust1['ext_cff_mss'].values
    MACaer[3, :] = FILE_dust2['ext_cff_mss'].values
    MACaer[4, :] = FILE_dust3['ext_cff_mss'].values
    MACaer[5, :] = FILE_dust4['ext_cff_mss'].values
    MACaer[6, :] = FILE_ash1['ext_cff_mss'].values
    MACaer[7, :] = FILE_GRISdust1['ext_cff_mss'].values
    MACaer[8, :] = FILE_GRISdust2['ext_cff_mss'].values
    MACaer[9, :] = FILE_GRISdust3['ext_cff_mss'].values
    MACaer[10, :] = FILE_GRISdustP1['ext_cff_mss'].values
    MACaer[11, :] = FILE_GRISdustP2['ext_cff_mss'].values
    MACaer[12, :] = FILE_GRISdustP3['ext_cff_mss'].values
    MACaer[13, :] = FILE_snw_alg['ext_cff_mss'].values
    MACaer[14, :] = FILE_glacier_algae1['ext_cff_mss'].values
    MACaer[15, :] = FILE_glacier_algae2['ext_cff_mss'].values

    Gaer = np.zeros([nbr_aer, nbr_wvl])

    Gaer[0, :] = FILE_soot1['asm_prm'].values
    Gaer[1, :] = FILE_soot2['asm_prm'].values
    Gaer[2, :] = FILE_dust1['asm_prm'].values
    Gaer[3, :] = FILE_dust2['asm_prm'].values
    Gaer[4, :] = FILE_dust3['asm_prm'].values
    Gaer[5, :] = FILE_dust4['asm_prm'].values
    Gaer[6, :] = FILE_ash1['asm_prm'].values
    Gaer[7, :] = FILE_GRISdust1['asm_prm'].values
    Gaer[8, :] = FILE_GRISdust2['asm_prm'].values
    Gaer[9, :] = FILE_GRISdust3['asm_prm'].values
    Gaer[10, :] = FILE_GRISdustP1['asm_prm'].values
    Gaer[11, :] = FILE_GRISdustP2['asm_prm'].values
    Gaer[12, :] = FILE_GRISdustP3['asm_prm'].values
    Gaer[13, :] = FILE_snw_alg['asm_prm'].values
    Gaer[14, :] = FILE_glacier_algae1['asm_prm'].values
    Gaer[15, :] = FILE_glacier_algae2['asm_prm'].values

    # load mass concentrations per layer into numpy array (one row per layer, one column per umpurity)
    # and convert to kg/kg unit

    MSSaer = np.zeros([nbr_lyr, nbr_aer])
    MSSaer[0:nbr_lyr, 0] = mss_cnc_soot1
    MSSaer[0:nbr_lyr, 1] = mss_cnc_soot2
    MSSaer[0:nbr_lyr, 2] = mss_cnc_dust1
    MSSaer[0:nbr_lyr, 3] = mss_cnc_dust2
    MSSaer[0:nbr_lyr, 4] = mss_cnc_dust3
    MSSaer[0:nbr_lyr, 5] = mss_cnc_dust4
    MSSaer[0:nbr_lyr, 6] = mss_cnc_ash1
    MSSaer[0:nbr_lyr, 7] = mss_cnc_GRISdust1
    MSSaer[0:nbr_lyr, 8] = mss_cnc_GRISdust2
    MSSaer[0:nbr_lyr, 9] = mss_cnc_GRISdust3
    MSSaer[0:nbr_lyr, 10] = mss_cnc_GRISdustP1
    MSSaer[0:nbr_lyr, 11] = mss_cnc_GRISdustP2
    MSSaer[0:nbr_lyr, 12] = mss_cnc_GRISdustP3
    MSSaer[0:nbr_lyr, 13] = mss_cnc_snw_alg
    MSSaer[0:nbr_lyr, 14] = mss_cnc_glacier_algae1
    MSSaer[0:nbr_lyr, 15] = mss_cnc_glacier_algae2

    MSSaer = MSSaer * 1e-9

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
    SSA_sum = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau = np.zeros([nbr_lyr, nbr_wvl])
    SSA = np.zeros([nbr_lyr, nbr_wvl])
    g = np.zeros([nbr_lyr, nbr_wvl])
    L_aer = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau_aer = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau_sum = np.zeros([nbr_lyr, nbr_wvl])
    SSA_sum = np.zeros([nbr_lyr, nbr_wvl])
    L_snw = np.zeros(nbr_lyr)
    tau_snw = np.zeros([nbr_lyr, nbr_wvl])
    direct = np.zeros([nbr_lyr, nbr_wvl])
    F_net = np.zeros([nbr_lyr, nbr_wvl])
    F_btm_net = np.zeros([1, nbr_wvl])
    F_top_net = np.zeros([1, nbr_wvl])
    intensity = np.zeros([nbr_lyr, nbr_wvl])
    F_top_pls = np.zeros([1, nbr_wvl])
    F_up = np.zeros([nbr_lyr, nbr_wvl])
    F_down = np.zeros([nbr_lyr, nbr_wvl])
    F_net2 = np.zeros([nbr_lyr, nbr_wvl])
    intensity2 = np.zeros([nbr_lyr, nbr_wvl])
    intensity2_top = np.zeros(nbr_wvl)
    F_abs = np.zeros([nbr_lyr, nbr_wvl])
    abs_vis = np.zeros(nbr_lyr)
    abs_nir = np.zeros(nbr_lyr)

    # for each layer, the layer mass (L) is density * layer thickness
    # for each layer the optical depth is the layer mass * the mass extinction coefficient
    # first for the ice in each layer
    for i in range(nbr_lyr):
        L_snw[i] = rho_snw[i] * dz[i]
        tau_snw[i, :] = L_snw[i] * MAC_snw[i, :]

    # then for the LAPs in each layer
    for i in range(nbr_lyr):
        for j in range(nbr_aer):
            L_aer[i, j, :] = np.multiply(L_snw[i], MSSaer[i, j])
            tau_aer[i, j, :] = np.multiply(L_aer[i, j, :], MACaer[j, :])

            tau_sum = tau_sum + tau_aer[i, j, :]
            SSA_sum = SSA_sum + (tau_aer[i, j, :] * SSAaer[j, :])
            g_sum = g_sum + (tau_aer[i, j, :] * SSAaer[j, :] * Gaer[j, :])

    # finally, for each layer calculate the effective SSA, tau and g for the snow+LAP
    for i in range(nbr_lyr):
        tau[i, :] = tau_sum[i, :] + tau_snw[i, :]
        SSA[i, :] = (1 / tau[i, :]) * (SSA_sum[i, :] + SSA_snw[i, :] * tau_snw[i, :])
        g[i, :] = (1 / (tau[i, :] * (SSA[i, :]))) * (g_sum[i, :] + (g_snw[i, :] * SSA_snw[i, :] * tau_snw[i, :]))


    g_star = g
    SSA_star = SSA
    tau_star = tau

    return flx_slr, g_star, SSA_star, tau_star, wvl