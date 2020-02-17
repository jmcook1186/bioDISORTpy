def snicar8d_GO(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer,
mss_cnc_soot1, mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, 
mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, 
mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, mss_cnc_snw_alg, mss_cnc_glacier_algae1, 
mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2, FILE_dust3, FILE_dust4, 
FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):


    import numpy as np
    import xarray as xr
    import matplotlib.pyplot as plt

    # set working directory (location of netcdf library)
    dir_base = "/home/joe/Code/SNICAR_NeuralNet/"
    dir_alg = "Data/Algal_Optical_Props/"
    dir_GO_files = "Data/GO_files/"

    # retrieve wavelength from arbitrary choice of netcdf file
    temp = xr.open_dataset(str(dir_base+dir_GO_files+"ice_geom_5000_5000.nc"))
    wvl = np.array(temp['wvl'].values)
    wvl = wvl*1e6
    nbr_wvl = len(wvl)

    # set reflectance of underlying surface
    R_sfc = [R_sfc for _ in range(nbr_wvl)]
    R_sfc = np.array(R_sfc)

    # Incoming Irradiance
    # calculate mu_not
    # mu_not = np.cos((slr_znt / 360) * 2 * np.pi() # convert radians if required
    mu_not = coszen
    flx_slr = []

    if DIRECT:

        with open(str(dir_base + dir_GO_files + "mlw_sfc_flx_frc_clr.txt")) as file:
            for line in file:
                line = float(line.rstrip("\n"))
                flx_slr.append(line)
        flx_slr = np.array(flx_slr)
        flx_slr[flx_slr==0]=1e-30
        Fs = flx_slr / (mu_not * np.pi)

        Fd = np.zeros(nbr_wvl)

    else:

        with open(str(dir_base + dir_GO_files + "mlw_sfc_flx_frc_cld.txt")) as file:
            for line in file:
                line = float(line.rstrip("\n"))
                flx_slr.append(line)
        
        flx_slr = np.array(flx_slr)
        flx_slr[flx_slr==0]=1e-30

        Fd = [flx_slr[i]/mu_not*np.pi for i in range(nbr_wvl)]
        Fs = np.zeros(nbr_wvl)

    # Read in ice optical properties
    # set string stubs for reading in ice optical files

    fl_stb1 = "ice_geom_"
    fl_stb2 = ".nc"

    #set up empty arrays
    SSA_snw = np.empty([nbr_lyr, nbr_wvl])
    MAC_snw = np.empty([nbr_lyr, nbr_wvl])
    g_snw = np.empty([nbr_lyr, nbr_wvl])

    for i in np.arange(0,nbr_lyr,1):

        if (side_length[i] == 0) | (depth[i] == 0):

            print("ERROR: ICE GRAIN LENGTH AND/OR DEPTH SET TO ZERO")

        else:

            s1 = str(side_length[i])
            s2 = str(depth[i])
            FILE_ice = str(dir_base+dir_GO_files+fl_stb1+s1+"_"+s2+fl_stb2)
   

    # read in single scattering albedo, MAC and g for ice crystals in each layer
        with xr.open_dataset(FILE_ice) as temp:
            SSA = temp['ss_alb'].values
            SSA_snw[i,:] = SSA

            ext_cff_mss = temp['ext_cff_mss'].values
            MAC_snw[i,:] = ext_cff_mss

            asm_prm = temp['asm_prm'].values
            g_snw[i,:] = asm_prm

    # open netcdf files
    FILE_soot1 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_soot1))
    FILE_soot2 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_soot2))
    FILE_dust1 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_dust1))
    FILE_dust2 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_dust2))
    FILE_dust3 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_dust3))
    FILE_dust4 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_dust4))
    FILE_ash1 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_ash1))
    FILE_GRISdust1 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_GRISdust1))
    FILE_GRISdust2 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_GRISdust2))
    FILE_GRISdust3 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_GRISdust3))
    FILE_GRISdustP1 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_GRISdustP1))
    FILE_GRISdustP2 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_GRISdustP2))
    FILE_GRISdustP3 = xr.open_dataset(str(dir_base + dir_GO_files+ FILE_GRISdustP3))
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
    SSA_sum = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau = np.zeros([nbr_lyr, nbr_wvl])
    SSA = np.zeros([nbr_lyr, nbr_wvl])
    g = np.zeros([nbr_lyr, nbr_wvl])
    L_aer = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau_aer = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau_sum = np.zeros([nbr_lyr, nbr_wvl])
    SSA_sum = np.zeros([nbr_lyr, nbr_wvl])
    L_snw = np.zeros(nbr_lyr)
    tau_snw = np.zeros([nbr_lyr,nbr_wvl])
    direct = np.zeros([nbr_lyr,nbr_wvl])
    F_net = np.zeros([nbr_lyr, nbr_wvl])
    F_btm_net = np.zeros([1,nbr_wvl])
    F_top_net = np.zeros([1,nbr_wvl])
    intensity = np.zeros([nbr_lyr, nbr_wvl])
    F_top_pls = np.zeros([1,nbr_wvl])
    F_up = np.zeros([nbr_lyr,nbr_wvl])
    F_down = np.zeros([nbr_lyr,nbr_wvl])
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
        tau[i,:] = tau_sum[i,:] + tau_snw[i,:]
        SSA[i,:] = (1/tau[i,:]) * (SSA_sum[i,:] + SSA_snw[i,:] * tau_snw[i,:])
        g[i, :] = (1 / (tau[i, :] * (SSA[i, :]))) * (g_sum[i,:] + (g_snw[i, :] * SSA_snw[i, :] * tau_snw[i, :]))


    ############################################
    # PERFORM DELTA TRANSFORMATION IF REQUIRED
    ############################################
    # The star represents the delta transformed quantity
    # if no delat transformation is applied, the starred quantity
    # is equal to the unstarred quantity

    if DELTA:
        g_star = g/(1+g)
        SSA_star = ((1-(g**2))*SSA)/(1-(SSA*(g**2)))
        tau_star = (1-(SSA*(g**2)))*tau

    else:
        g_star = g
        SSA_star = SSA
        tau_star = tau


    # CALCULATE TOTAL OPTICAL DEPTH OF ENTIRE COLUMN
    # i.e. tau_clm = total optical depth from upper boundary 
    # to upper boundary of layer n. This is therefore a cumulative
    # quantity - subsequently lower layers contain the sum of the
    # # optical depth of all overlying layers

    tau_clm = np.zeros([nbr_lyr,nbr_wvl])
    for i in np.arange(1,nbr_lyr,1):
        #start loop from 2nd layer, i.e. index = 1
        tau_clm[i,:] = tau_clm[i-1,:]+tau_star[i-1,:]

    # SET BOUNDARY CONDITION: BOTTOM BOUNDARY
    # calculate radiation reflected skywards by underlying surface (i.e. lower model boundary)
    # remainder is lost

    S_sfc = R_sfc * mu_not * np.exp(-(tau_clm[nbr_lyr-1,:] + tau_star[nbr_lyr-1,:])/mu_not)*np.pi * Fs

    ######################################################
    # Apply Two-Stream Approximation (Toon et al, table 1)
    ######################################################
    """
    Three 2-stream approximations are available: Eddington,
    Quadrature and hemispheric mean. The equations for each
    approximation are provided in Toon et al. (1989) Table 1.

    The hemispheric mean scheme is derived by assuming that the
    phase function is equal to 1  + g  in the forward scattering 
    hemisphere and to 1  - g  in the backward scattering hemisphere. 
    The asymmetry parameter is g. The hemispheric mean is only
    useful for infrared wavelengths

    """

    if APRX_TYP == 1:
        #apply Eddington approximation
        gamma1 = (7-(SSA_star * (4+(3*g_star))))/4
        gamma2 = -(1-(SSA_star*(4-(3*g_star))))/4
        gamma3 = (2-(3*g_star*mu_not))/4
        gamma4 = 1-gamma3
        mu_one = 0.5

    elif APRX_TYP==2:
        #apply quadrature approximation
        gamma1 = np.sqrt(3)*(2-(SSA_star*(1+g_star)))/2
        gamma2 = SSA_star * np.sqrt(3)*(1-g_star)/2
        gamma3 = (1-(np.sqrt(3)*g_star*mu_not))/2
        gamma4 = 1-gamma3
        mu_one = 1/np.sqrt(3)

    elif APRX_TYP==3:
        #apply hemispheric mean approximation
        gamma1 = 2 - (SSA_star*(1+g_star))
        gamma2 = SSA_star*(1-g_star)
        gamma3 = (1-(np.sqrt(3) * g_star*mu_not))/2
        gamma4 = 1-gamma3
        mu_one = 0.5


    # Toon et al equation 21 and 22
    # Note that the values of lam and GAMMA depend upon gamma1 and gamma2, which
    # vary depending upon the two-stream approximation used
    # variable "lambda" renamed "lam" to avoid confusion with lambda function
    lam = np.sqrt(abs((gamma1**2)-(gamma2**2)))
    GAMMA = gamma2/(gamma1+lam)

    # calculate coefficients required for tridiagonal matrix calculation
    # (Toon et al Equation 44)
    e1 = 1+(GAMMA*np.exp(-lam*tau_star))
    e2 = 1-(GAMMA*np.exp(-lam*tau_star))
    e3 = GAMMA+np.exp(-lam*tau_star)
    e4 = GAMMA-np.exp(-lam*tau_star)


    ######################################
    # Calculate C-functions
    ######################################

    # C is the direct beam flux calculated at the top and bottom of each layer, i,
    # see Toon equations 23 and 24

    """ N.B. consider adding in stability check here as per Flanner's Matlab code """

    C_pls_btm = np.zeros([nbr_lyr,nbr_wvl])
    C_mns_btm = np.zeros([nbr_lyr,nbr_wvl])
    C_pls_top = np.zeros([nbr_lyr, nbr_wvl])
    C_mns_top = np.zeros([nbr_lyr, nbr_wvl])

    for i in np.arange(0,nbr_lyr,1):

        if np.sum(Fs) > 0.0:

            C_pls_btm[i,:] = (SSA_star[i,:]*np.pi*Fs*np.exp(-(tau_clm[i,:]+tau_star[i,:])/mu_not)*
                                (((gamma1[i,:]-(1/mu_not))*gamma3[i,:])+(gamma4[i,:]*gamma2[i,:])))\
                                /((lam[i,:]**2)-(1/(mu_not**2)))

            C_mns_btm[i,:] = (SSA_star[i,:]*np.pi*Fs*
                                np.exp(-(tau_clm[i,:]+tau_star[i,:])/mu_not) * (((gamma1[i,:]+(1/mu_not))*gamma4[i,:])+
                                (gamma2[i,:]*gamma3[i,:])))/((lam[i,:]**2)-(1/mu_not**2))

            C_pls_top[i,:] = (SSA_star[i,:] * np.pi * Fs * np.exp(-tau_clm[i,:]/mu_not)* ((gamma1[i,:] - (1/mu_not))
                                * gamma3[i,:] + (gamma4[i,:]*gamma2[i,:])))/((lam[1,:]**2)-(1/mu_not**2))

            C_mns_top[i,:] = (SSA_star[i,:] * np.pi * Fs * np.exp(-tau_clm[i,:]/mu_not) * ((gamma1[i,:]+(1/mu_not))
                                * gamma4[i,:] + (gamma2[i,:]*gamma3[i,:])))/((lam[i,:]**2)-(1/mu_not**2))

        else:
            # no direct-beam flux:
            C_pls_btm[i,:] = 0
            C_mns_btm[i,:] = 0
            C_pls_top[i,:] = 0
            C_mns_top[i,:] = 0


    # Toon equations 41-43.
    # Boundary values for i=1 and i=2nbr_lyr, specifics for i=odd and i=even
    # Set up lists
    A = np.zeros([2*nbr_lyr,nbr_wvl])
    B = np.zeros([2*nbr_lyr,nbr_wvl])
    D = np.zeros([2*nbr_lyr,nbr_wvl])
    E = np.zeros([2*nbr_lyr,nbr_wvl])

    ###########################################
    # Initialize tridiagonal matrix solution
    ###########################################

    # expanding the number of layers to 2*nbr_lyr so that fluxes at upper and lower
    # layer boundaries can be resolved. This section was confusing to code - for each layer
    # index (n) a second pair of indices (2 x i) are required. Different solutions are
    # applied depending upon whether i is even or odd. To translate the indexing for this 
    # from FORTRAN/MATLAB into Python, it was necessary to assert n = (i/2)-1 for even layers
    # and n = floor(i/2) for odd layers, with specific rules for the boundaries i = 0 and 
    # i = nbr_lyrs-1 (i.e. top surface and bottom surface).

    for i in np.arange(0,2*nbr_lyr,1):
 
        #TOP LAYER    
        if i==0:
            A[0,:] = 0.0
            B[0,:] = e1[0,:]
            D[0,:] = -e2[0,:]
            E[0,:] = Fd-C_mns_top[0,:]

        # BOTTOM LAYER
        elif i== 2*nbr_lyr-1:
            A[i,:] = e1[nbr_lyr-1,:]-(R_sfc * e3[nbr_lyr-1,:])
            B[i,:] = e2[nbr_lyr-1,:]-(R_sfc * e4[nbr_lyr-1,:])
            D[i,:] = 0.0
            E[i,:] = S_sfc[:] - C_pls_btm[nbr_lyr-1,:] + (R_sfc * C_mns_btm[nbr_lyr-1,:])


        # EVEN NUMBERED LAYERS
        elif i%2==0:
            n = int(i/2)-1
            A[i,:] = (e2[n,:] * e3[n,:])-(e4[n,:] * e1[n,:])
            B[i,:] = (e1[n,:] * e1[n+1,:])-(e3[n,:] * e3[n+1,:])
            D[i,:] = (e3[n,:] * e4[n+1,:])-(e1[n,:] * e2[n+1,:])
            E[i,:] = (e3[n,:] * (C_pls_top[n+1,:] - C_pls_btm[n,:])) +  (e1[n,:] * (C_mns_btm[n,:] - C_mns_top[n+1,:]))

        # ODD NUMBERED LAYERS
        elif (i%2 ==1) and (i < 2*nbr_lyr-1):

            n = int(np.floor(i/2))
            A[i,:] = (e2[n+1,:] * e1[n,:])-(e3[n,:] * e4[n+1,:])
            B[i,:] = (e2[n,:] * e2[n+1,:])-(e4[n,:] * e4[n+1,:])
            D[i,:] = (e1[n+1,:] * e4[n+1,:])-(e2[n+1,:] * e3[n+1,:])
            E[i,:] = (e2[n+1,:] * (C_pls_top[n+1,:] - C_pls_btm[n,:])) + (e4[n+1,:] * (C_mns_top[n+1,:] - C_mns_btm[n,:]))

    # Now the actual tridiagonal matrix solving. Simply dividing A/B and E/B 
    # throws an exception due to division by zero. Here we use numpy's nan_to_num
    # function to achieve the division where possible and replace nans with zeros.
    # We also set numpy to ignore the division error.

    # for bottom layer only
    # Toon et al Eq 45
    AS = np.zeros([2*nbr_lyr,nbr_wvl])
    DS = np.zeros([2*nbr_lyr,nbr_wvl])

    np.seterr(divide='ignore',invalid='ignore')
    AS[2*nbr_lyr-1,:] = np.nan_to_num(A[2*nbr_lyr-1,:]/B[2*nbr_lyr-1,:])
    DS[2*nbr_lyr-1,:] = np.nan_to_num(E[2*nbr_lyr-1,:]/B[2*nbr_lyr-1,:])

    # for all layers above bottom layer, starting at second-to-bottom and progressing towards
    # surface:
    # Toon et al Eq 46
    X = np.zeros([nbr_lyr*2,nbr_wvl])
    for i in np.arange(2*nbr_lyr-2,-1, -1):
        X[i,:] = 1/(B[i,:]-(D[i,:] * AS[i+1,:]))
        AS[i,:] = np.nan_to_num(A[i,:]*X[i,:])
        DS[i,:] = np.nan_to_num((E[i,:]-(D[i,:]*DS[i+1,:]))*X[i,:])

    # then for all layers, progressing from surface to bottom
    # Toon et al Eq 47
    Y = np.zeros([nbr_lyr*2,nbr_wvl])

    for i in np.arange(0,2*nbr_lyr,1):
        if i ==0:
            Y[0,:] = DS[0,:]
        else:
            Y[i,:] = DS[i,:] - (AS[i,:]*Y[i-1,:])

    
    #############################################################
    # CALCULATE DIRECT BEAM FLUX AT BOTTOM OF EACH LAYER

    # loop through layers
    for i in np.arange(0,nbr_lyr,1):

        # (Toon et al. eq 50)
        direct[i,:] = mu_not * np.pi * Fs * np.exp(-(tau_clm[i,:] + tau_star[i,:]) / mu_not)

        # net flux (positive upward = F_up - F_down) at the base of each layer (Toon et al. Eq 48)
        F_net[i,:] = (Y[2*i,:] * (e1[i,:]-e3[i,:])) + (Y[2*i+1,:] * (e2[i,:] - e4[i,:])) + C_pls_btm[i,:] - C_mns_btm[i,:] - direct[i,:]

        # mean intensity at the base of each layer (Toon et al. Eq 49)
        intensity[i,:] = (1/mu_one) * (Y[2*i,:] * (e1[i,:] + e3[i,:]) + Y[2*i+1,:] * (e2[i,:] + e4[i,:]) + C_pls_btm[i,:] + C_mns_btm[i,:]) + (direct[i,:]/mu_not)
        intensity[i, :] = intensity[i, :] / (4 * np.pi)


    # Upward flux at upper model boundary (Toon et al Eq 31)
    F_top_pls = (Y[0,:] * (np.exp(-lam[0,:] * tau_star[0,:]) + GAMMA[0,:])) + (Y[1,:] * (np.exp(-lam[0,:] * tau_star[0,:])-GAMMA[0,:])) + C_pls_top[0,:]


    for i in np.arange(0,nbr_lyr,1):
        # Upward flux at the bottom of each layer interface (Toon et al. Eq31)
        F_up[i,:] = Y[2*i,:] * (np.exp(0) + GAMMA[i,:] * np.exp(-lam[i,:] * tau_star[i,:])) + Y[2*i+1,:] * (np.exp(0) - GAMMA[i,:] * np.exp(-lam[i,:] * tau_star[i,:])) + C_pls_btm[i,:]

        # Downward flux at the bottom of each layer interface (Toon et al. Eq32) plus direct beam component
        F_down[i,:] = Y[2*i,:] * (GAMMA[i,:] * np.exp(0) + np.exp(-lam[i,:] * tau_star[i,:])) + Y[2*i+1,:] * (GAMMA[i,:] * np.exp(0) - np.exp(-lam[i,:] * tau_star[i,:])) + C_mns_btm[i,:] + direct[i,:]

        # Derived net (upward-downward) flux (should equal F_net)
        F_net2[i,:] = F_up[i,:] - F_down[i,:]

        intensity2[i,:] = F_up[i,:] + F_down[i,:]
    
    # surface planar intensity
    intensity2_top[:] = F_top_pls + ((mu_not * np.pi * Fs) + Fd)

    # Net flux at lower model boundary = bulk transmission through entire media
    # = energy absorbed by underlying surface
    F_btm_net[0,:] = -F_net[nbr_lyr-1,:]

    # Hemispheric wavelength-dependent albedo
    albedo = F_top_pls/ ((mu_not * np.pi * Fs)+ Fd)

    # Net flux at upper model boundary
    F_top_net[0,:] = F_top_pls - ((mu_not * np.pi * Fs) + Fd)

    # absorbed flux in each layer (negative if there is net emission (bnd_typ = 4))
    for i in np.arange(0,nbr_lyr,1):
        if i ==0:
            F_abs[0,:] = F_net[0,:]-F_top_net
        else:
            F_abs[i,:] = F_net[i,:] - F_net[i-1,:]

    # set indices for constraining calculations to VIS and NIR bands
    vis_max_idx = 39
    nir_max_idx = len(wvl)

    # Spectrally-integrated absorption in each layer:
    abs_slr = np.sum(F_abs,axis=1)

    for i in np.arange(0,nbr_lyr,1):
        abs_vis[i] = np.sum(F_abs[i,0:vis_max_idx])
        abs_nir[i] = np.sum(F_abs[i,vis_max_idx:nir_max_idx])

    # Spectrally - integrated absorption by underlying surface:
    abs_slr_btm = sum(np.squeeze(F_btm_net))
    abs_vis_btm = sum(np.squeeze(F_btm_net[0:vis_max_idx]))
    abs_nir_btm = sum(np.squeeze(F_btm_net[0,vis_max_idx:nir_max_idx]))

    # Calculate radiative heating rate in kelvin per second.
    # Multiply by 3600 to convert to K per hour
    # specfic heta capacity of ice = 2117 J kg-1 K-1
    heat_rt = abs_slr / (L_snw * 2117) # [K / s]
    heat_rt = heat_rt * 3600 # [K / hr]

    # Energy conservation check:
    # % Incident direct + diffuse radiation equals(absorbed + transmitted + bulk_reflected)
    energy_sum = (mu_not * np.pi * Fs) + Fd - (sum(F_abs) + F_btm_net + F_top_pls)

    # spectrally-integrated terms:
    # energy conservation total error
    energy_error = abs(np.sum(energy_sum))

    if energy_error > 1e-10:
        energy_conservation_error = np.sum(abs(energy_sum))
        print(f"CONSERVATION OF ENERGY ERROR OF {energy_conservation_error}")

    ######################################
    # Re-alias results for outputting
    ######################################

    # total incident insolation(Wm - 2)
    total_insolation = np.sum((mu_not * np.pi * Fs) + Fd)

    # energy absorbed by all snow layers
    abs_slr_tot = np.sum(np.sum(F_abs))

    # energy absorbed by underlying substrate
    energy_abs_under_sfc = np.sum(F_btm_net)

    # Spectrally - integrated solar, visible, and NIR albedos:
    BBA = np.sum(flx_slr * albedo) / np.sum(flx_slr)

    BBAVIS = sum(flx_slr[0:vis_max_idx]*albedo[0:vis_max_idx])/ sum(flx_slr[0:vis_max_idx])

    BBANIR = sum(flx_slr[vis_max_idx:nir_max_idx]*albedo[vis_max_idx: nir_max_idx]) / sum(flx_slr[vis_max_idx:nir_max_idx])

    # % Spectrally - integrated VIS and NIR total snowpack absorption:
    abs_vis_tot = sum(flx_slr[0:vis_max_idx]*(1 - albedo[0:vis_max_idx]))
    abs_nir_tot = sum(flx_slr[vis_max_idx:nir_max_idx]*(1 - albedo[vis_max_idx:nir_max_idx]))

    return wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation
