import numpy as np
import glob
import os
import os.path
import math
import constant
import csv
import errno


######################################
#   read_model_input()
#   - read "model_info.txt" file containing the initial setup for key parameters
#   - inputdata_file_name : input file containing the observational data for TDEs, either in .txt or .csv format
#                           read "2. Input data format" in READ_ME file
#   - output_file_name     : the name of the output file name consisting the inferred mass. If "-" is given, the output file name is given with the candidate names. If there is one candidate listed in the input file, output_file_name = "candidate-name.txt". If there are more than two, output_file_name = "first-candidate-name_last-candidate-name.txt"
#   - c1                  : value of c_1 (the characteristic distance scale of the emission region in units of the apocenter distance of the most tightly bound debris). Defualt = 1
#   - Del_omega           : Solid angle (in units of pi) of radiation from the emission region. Default = 2.0
#   - N_sampling          : The number of bin in mbh (black hole mass) and mstar (stellar mass) to find the range of mbh and mstar for given Lobs and Tobs. And it determines the grid resolution of the contour plots produced as output. If N_sampling is smaller, it will produce output more quickely.
######################################
def read_model_input():
    inputfile = "input_file/model_info.txt"
    if os.path.isfile(inputfile):
        print ("model_info.txt exists!")
        file=open(inputfile,"r")
        file_line = file.readlines()
        c1 = 1.0
        Del_omega = 2.0 * math.pi
        N_sampling = 300
        for line in file_line:
            if not line.isspace():
                p=line.split()
                if(p[0]=="inputdata_file_name"):
                    inputdata_file_name = p[1]
                if(p[0]=="output_file_name"):
                    if(p[1] =="#" or p[1]=="-"):
                        output_file_name = ""
                    else:
                        output_file_name = p[1]

                if(p[0]=="c1"):
                    c1 = float(p[1])
                if(p[0]=="Del_omega"):
                    Del_omega = float(p[1])*math.pi
                if(p[0]=="N_sampling"):
                    N_sampling = int(p[1])
        return inputdata_file_name,output_file_name, c1, Del_omega,N_sampling


######################################
#   read_input_data()
#   - routine that reads Lobs and Tobs from the input file designated in "model_info.txt"
######################################

def read_input_data(inputdata_file_name):
    if(".csv" in inputdata_file_name):
        format = "csv"
    elif(".txt" in inputdata_file_name):
        format = "txt"
    else:
        print ("This format is not supported in this version")
        exit()
    
    
    
    index_array = []
    Lpeak_array = [[],[]]
    Tpeak_array = [[],[]]
   # MBH_array = [[],[]]
   # Mhost_array=[[],[]]
   # time_array=[]
   # which_bulge_relation =[]
    tot_input_variable_count = 4
    min_error = 0.0
    mstar_range=[]
    mbh_range=[]

    linecount = 0
    nan_check=[]
    check_input = [0,0,0,0,0,0,0]
    
    input_array=[]
    if(format =="csv"):
        with open("input_file/"+inputdata_file_name, newline='') as csvfile:
            data = csv.reader(csvfile, delimiter=',')
            for line in data:
                p = line#.split()
                if("candidate_name" not in p[0]):
                    index_array.append(p[0])
                    mstar_range.append([float(p[1]),float(p[2])])
                    mbh_range.append([float(p[3]),float(p[4])])
    
                    Lpeak_array = add_array(5, p, Lpeak_array,min_error)
                    Tpeak_array = add_array(8, p, Tpeak_array,min_error)
                   # time_array = add_array_time(11, p, time_array,min_error)
                   # MBH_array = add_array(13, p, MBH_array, min_error)
                   # Mhost_array = add_array(16,p,Mhost_array,min_error)
                   # which_bulge_relation.append(int(p[17]))
                    linecount = linecount +1
    else:
        if os.path.isfile("input_file/"+inputdata_file_name):
            data = open("input_file/"+inputdata_file_name,"r").readlines()
            for line in data:
                if not line.isspace():
                    p = line.split()
                    if("candidate_name" not in p[0]):

                        index_array.append(p[0])
                        mstar_range.append([float(p[1]),float(p[2])])
                        mbh_range.append([float(p[3]),float(p[4])])
                        Lpeak_array = add_array(5, p, Lpeak_array,min_error)
                        Tpeak_array = add_array(8, p, Tpeak_array,min_error)
                  #      time_array = add_array_time(11, p, time_array,min_error)
                  #      MBH_array = add_array(13, p, MBH_array, min_error)
                  #      Mhost_array = add_array(16,p,Mhost_array,min_error)
                  #      which_bulge_relation.append(int(p[17]))
                        

                        linecount = linecount +1


    return index_array,mstar_range,mbh_range, Lpeak_array,Tpeak_array,linecount



######################################
#   reference_mbh()
#   - the reference for M_BH - bulge correlation
######################################

def reference_mbh(which):
    bulge_correlation_reference=["McConnell&Ma2013","Magorrian+1998","Haring&Rix(2004)", "Kormendy&Ho(2004)","Reines&Marta (2015,AGN)","Reines&Marta(2015)","avorgnan(2016)","Ferrarese&Ford(2005)", "Merritt&Ferrarese(2001)","Gultekin+2009","Rothberg&Fischer(2010)"]
    
    return bulge_correlation_reference[which]
    

######################################
#   make_sure_path_exists()
#   - make directory "output" in which output files are saved.
######################################

def make_sure_path_exists(path):
  try:
    os.makedirs(path)
  except OSError as exception:
    if exception.errno != errno.EEXIST:
      raise


######################################
#   add_array()
#   - a routine to save input data into an array
######################################

def add_array(which, p, array,min_error):
  if(p[which]=="-"):
    array[0].append(0.0)
    array[1].append([0.0,0.0])
  else:
    
    value = float(p[which])
    array[0].append(value)
    if(len(p)>= which+1):
        if(p[which+1] =="-"):
            error1 = 0.0
        else:
            error1 = max(min_error*value,float(p[which+1]))
    else:
        error1 = 0.0
    if(len(p)>= which+2):
        if(p[which+2] =="-"):
            error2 = 0.0
        else:
            error2 = max(min_error*value,float(p[which+2]))
    else:
        error2 = 0.0
    error = [error1,error2]
    array[1].append(error)


  return array


######################################
#   add_array_time()
#   - a routine to save input data for peak time into an array
######################################

def add_array_time(which, p, array,min_error):
  if(p[which]=="-" or p[which+1]=="-"):
    array.append(0.0)
  else:
    value = float(p[which+1])-float(p[which])
    array.append(value)

  return array


######################################
#   get_Xi()
#   - a routine to calculate Xi
#   - input :
#             mbh = black hole mass in units of 10^6msun
#           mstar = stellar mass in units of 1msun
######################################

def get_Xi(mbh,mstar):
    
    Xi = (1.27 - 0.3 *(mbh)**0.242 )*((0.620 +
    np.exp((mstar - 0.674)/0.212))/(1.0 +
                                 0.553 *np.exp((mstar - 0.674)/0.212)))
    return Xi
  
    
######################################
#   rtidal()
#   - a routine to calculate conventional tidal radius (mbh/mstar)^(1/3)rstar in cgs units
######################################

def rtidal(mbh,mstar,rstar):
    rt = (mbh * 1e6/mstar)**(1.0/3.0) * rstar * constant.unit_rsun#in cm
    return rt
    
    
######################################
#   get_rstar()
#   - a routine to calculate stellar radius (in units of solar radius) from stellar mass using r_star = 0.93 * m_star^(8/9) (Ryu+2020, arXiv:2001.03502)
######################################

def get_rstar(mstar):
    rs = 0.93 * (mstar)**(8.0/9.0)
    return rs #in rsun
    
    
######################################
#   epsilon()
#   - a routine to calculate the conventional estimate of debris energy, G*M_BH*rstar/rtidal^2 in cgs unit
######################################

def epsilon(mbh,mstar,rstar):
    rt = rtidal(mbh,mstar,rstar)
    ee = constant.unit_G * (mbh * 1e6 * constant.unit_msun) * (rstar * constant.unit_rsun) / rt**2
    return ee #in erg/g
    
    
######################################
#   get_a0()
#   - a routine to calculate the apocenter distance of orbit for the most tightly bound debris with energy Xi * epsilon in cgs unit
######################################

def get_a0(mbh,mstar,rstar):
    energy = get_Xi(mbh,mstar) * epsilon(mbh,mstar,rstar)
    return constant.unit_G * (mbh * 1e6 * constant.unit_msun)/ energy #in cm


######################################
#   get_t0()
#   - a routine to calculate the characteristic time scale in cgs unit
######################################

def get_t0(mbh,mstar,rstar):
    a0 = get_a0(mbh,mstar,rstar)
    t0 = math.pi / np.sqrt(2.0) * a0**1.5 / np.sqrt(constant.unit_G * (mbh * 1e6 * constant.unit_msun))
    return t0


######################################
#   get_tpeak()
#   - a routine to calculate the peak time (1.5t_0) in cgs unit
######################################

def get_tpeak(mbh,mstar,rstar):
    t0 = get_t0(mbh,mstar,rstar)
    tpeak = (3.0/2.0)*t0
    return tpeak



######################################
#   get_mdot_max()
#   - a routine to calculate the peak mass fallback rate in cgs unit
######################################

def get_mdot_max(mbh,mstar,rstar):
    t0 = get_t0(mbh,mstar,rstar)
    mdotmax = mstar * constant.unit_msun / t0 / 3.0
    return mdotmax
    
    
######################################
#   get_Lmax()
#   - a routine to calculate the peak luminosity in cgs unit
######################################

def get_Lmax(mbh,mstar,rstar,c1):
    energy = get_Xi(mbh,mstar) * epsilon(mbh,mstar,rstar)
    mdotmax = get_mdot_max(mbh,mstar,rstar)
    Lmax = mdotmax * energy / c1
    
    return Lmax

    
######################################
#   get_v()
#   - a routine to calculate the linewidth in cgs unit
######################################

def get_v(mbh,mstar,rstar,c1):
    energy = get_Xi(mbh,mstar) * epsilon(mbh,mstar,rstar)
    v = np.sqrt(2.0 * energy / c1 )
    return v


######################################
#   get_v()
#   - a routine to calculate the blackbody temperature at peak luminosity in cgs unit
######################################

def get_Tmax(mbh,mstar,rstar,c1,del_omega):
    Lmax = get_Lmax(mbh,mstar,rstar,c1)
    a0 = get_a0(mbh,mstar,rstar)
    factor_denom = del_omega * constant.BS_constant * c1**2 * a0**2
    Tmax = (Lmax / factor_denom )**(1.0/4.0)
    return Tmax



######################################
#   get_Lobs()
#   - a routine to calculate Lobs in cgs unit
#   - in the current version, Lobs = Lmax
#   - the time-dependence will be included
######################################

def get_Lobs(mbh,mstar,rstar,c1,t_obs):
    Lmax = get_Lmax(mbh,mstar,rstar,c1)
    tpeak = t_obs#get_tpeak(mbh,mstar,rstar)
    Lobs = Lmax  * (t_obs/tpeak)**(-5.0/3.0)
    
    return Lobs

######################################
#   get_Tobs()
#   - a routine to calculate Tobs in cgs unit
#   - in the current version, Tobs = Tmax
#   - the time-dependence will be included
######################################

def get_Tobs(mbh,mstar,rstar,c1,del_omega,t_obs):
  Lobs = get_Lobs(mbh,mstar,rstar,c1,t_obs)
  a0 = get_a0(mbh,mstar,rstar)
  factor_denom = del_omega * constant.BS_constant * c1**2 * a0**2
  Tmax = (Lobs / factor_denom )**(1.0/4.0)
  return Tmax



######################################
#   get_mbh_mstar_from_L_T()
#   - a routine to calculate M_BH and _Mstar for the next iteration in solver1_LT
######################################

def get_mbh_mstar_from_L_T(Lobs,Tobs,mbh,mstar,rstar,c1,t_obs,delta_t,del_omega):
    Xi = get_Xi(mbh,mstar)
    tpeak = get_tpeak(mbh,mstar,rstar)
    x = (t_obs - delta_t)/tpeak
    mbh_next = 1.097091710461842 * (Tobs/10**4.5)**(-8.0/3.0)* (c1/0.5)**-2 * Xi **3* (del_omega/4.0/math.pi)**(-2.0/3.0)
    mstar_next = 1.46455136624986 * (Lobs/1e44)**(9.0/4.0) * (Tobs/10**4.5)**(-1) * Xi **(-4.5)*(c1/0.5)**1.5 * (del_omega/4.0/math.pi)**(-1.0/4.0)
    return mbh_next,mstar_next



######################################
#   solver1_LT()
#   - a routine to solve the non-linear equations for M_BH and _Mstar for given Lobs and Tobs.
#   - With the first guess, it computes M_BH and M_star using the routine get_mbh_mstar_from_L_T(). It uses the mean value of the first guess and the first computed values to get the second values. It repeats this procedure until they converge simultaneously.
#    If this method fails to find the correct solutions, the routine "solver2_LT()", which is a 2-dimensional newton-Raphson solver, is
######################################

def solver1_LT (Lobs, Tobs, mbh,mstar,c1,del_omega):
    mstar_mean = mstar
    mstar_new = mstar_mean
    mstar_old = mstar_mean
    mbh_mean = mbh
    mbh_new = mbh_mean
    mbh_old = mbh_mean
    keep_iterating =1
    error = 1.0
    TOL = 1e-12
    N_ITER = 1000
    n_iter = 0
    count_converge=0
    while(keep_iterating == 1):
        
        rstar_mean = get_rstar(mstar_mean)
        t_obs = get_tpeak(mbh_mean,mstar_mean,rstar_mean)
        delta_t = 0.0
        mbh_next, mstar_next = get_mbh_mstar_from_L_T(Lobs, Tobs, mbh_mean,mstar_mean,rstar_mean,c1,t_obs,delta_t,del_omega)
    
        mbh_mean = (mbh_next + mbh_old) * 0.5
        mstar_mean = (mstar_next + mstar_old) * 0.5
        
    
        error_bh = np.abs((mbh_next - mbh_old)/mbh_old)
        
        error_star = np.abs((mstar_next - mstar_old)/mstar_old)
        
        if(error_bh < TOL and error_star < TOL):
            
            
            count_converge = count_converge + 1
            n_iter = n_iter - 1
            if (count_converge==2):
                keep_iterating = 0
                retv = 0
        else:
            count_converge = 0
        mbh_old = mbh_next
        mstar_old = mstar_next
        n_iter = n_iter + 1
        
        if(n_iter==N_ITER):
            keep_iterating=0
            retv = 1
        
    return retv, mbh_mean, mstar_mean
    
    
######################################
#   solver2_LT()
#   - A solver using a 2-dimensional Newtonian-Raphson method
#   - retv : 0 -> found converging solutions
#          : 1 -> fail to find the solutions
######################################
def solver2_LT (Lobs, Tobs, mbh,mstar,c1,del_omega):
    x1 = mbh
    x2 = mstar
    error = 1.0
    TOL = 1e-11
    N_ITER = 100
    n_iter = 0
    count_converge=0
    keep_iterating = 1
    while(keep_iterating == 1):

        rstar = get_rstar(x2)
        t_obs = get_tpeak(x1,x2,rstar)
        t_peak = get_tpeak(x1,x2,rstar)
        delta_t = 0.0
        x = (t_obs - delta_t)/t_peak
        
        f1 = (Lobs - get_Lobs(x1,x2,rstar,c1,t_obs))/10**44
        f2 = (Tobs - get_Tobs(x1,x2,rstar,c1,del_omega,t_obs))/10**4.5
        jac = get_jacobian_LT(x1,x2,rstar,c1,x)
        det = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1]
        dx1 = 1.0/det * ( jac[1][1] * f1 - jac[0][1] * f2) * 0.5
        dx2 = 1.0/det * (-jac[1][0] * f1 + jac[0][0] * f2) * 0.5
        if( x < 1):
            print ("WARNING tobs < tpeak")
            break
    
        x1 = x1 + dx1
        x2 = x2 + dx2
        error = np.abs(dx1/x1)
        if(error < TOL):
            retv = 0
            keep_iterating=0

        if(n_iter ==N_ITER):
            retv = 1
            keep_iterating =0
        n_iter = n_iter + 1


    return retv, x1, x2


######################################
#   Mhost_Mbh_relation()
#   - different M_BH - bulge correlation
#   - which
#        0 = McConnel&Ma2013 : M_bh - M_bulge relation
#        1 =   : 1 -> fail to find the solutions
######################################
def Mhost_Mbh_relation(which, Mhost):
  
  if(which ==0):
    #McConnel&Ma2013 : M_bh - M_bulge relation https://ui.adsabs.harvard.edu/abs/2013ApJ...764..184M/abstract
    alpha = 8.46
    alpha_scatter = 0.05
    beta = 1.05
    beta_scatter = 0.11
    median, min_val, max_val = find_median_min_max(alpha,alpha_scatter,beta,beta_scatter,Mhost/1e11)
    return median, min_val, max_val
  elif(which==1):
  #Magorrian+1998 : M_bh - M_bulge relation https://iopscience.iop.org/article/10.1086/300353/pdf
    alpha =-1.79
    alpha_scatter = 1.35
    beta = 0.96
    beta_scatter = 0.12
    median, min_val, max_val = find_median_min_max(alpha,alpha_scatter,beta,beta_scatter,Mhost)
    return median, min_val, max_val


  elif(which == 2):
    #Haring&Rix2004 : M_bh - M_bulge relation https://iopscience.iop.org/article/10.1086/383567/pdf
    alpha =8.20
    alpha_scatter = 0.10
    beta = +1.12
    beta_scatter = 0.06
    median, min_val, max_val = find_median_min_max(alpha,alpha_scatter,beta,beta_scatter,Mhost/1e11)
    return median, min_val, max_val
  elif(which == 3):
    #Kormendy&Ho2004 : M_bh - M_bulge relation https://arxiv.org/pdf/1304.7762.pdf
    alpha =8.6902
    alpha_scatter = 0.05
    beta = +1.16
    beta_scatter = 0.08
    median, min_val, max_val = find_median_min_max(alpha,alpha_scatter,beta,beta_scatter,Mhost/1e11)
    return median, min_val, max_val

  elif(which == 4):
#Reines&Marta 2015 : M_bh - M_galaxy relation (AGN) https://arxiv.org/pdf/1508.06274.pdf
    alpha =7.45
    alpha_scatter = 0.08
    beta = 1.05
    beta_scatter = 0.11
    median, min_val, max_val = find_median_min_max(alpha,alpha_scatter,beta,beta_scatter,Mhost/1e11)
    return median, min_val, max_val

  elif(which == 5):
#Reines&Marta 2015 : M_bh - M_galaxy relation (AGN) https://arxiv.org/pdf/1508.06274.pdf
    alpha =8.95
    alpha_scatter = 0.09
    beta = 1.40
    beta_scatter = 0.21
    median, min_val, max_val = find_median_min_max(alpha,alpha_scatter,beta,beta_scatter,Mhost/1e11)
    return median, min_val, max_val

  elif(which == 6):
#Savorgnan 2016 : nsph - M_BH relation (AGN) https://iopscience.iop.org/article/10.3847/0004-637X/821/2/88/pdf
    alpha =8.15
    alpha_scatter = 0.06
    beta = 3.37
    beta_scatter = 0.15
    median, min_val, max_val = find_median_min_max2(alpha,alpha_scatter,beta,beta_scatter,Mhost)
    return median, min_val, max_val
  elif(which ==7):
#Ferrarese&Ford (2005): nsph - sigma relation (AGN) https://arxiv.org/pdf/astro-ph/0411247.pdf
    #Mbh/10^8 = (1.66+-0.24)*(sigma/200km/s)^(4.86+-0.43)
    alpha =0.0
    alpha_scatter = 0.0
    beta = 0.0
    beta_scatter = 0.0
    median, min_val, max_val = find_median_min_max2(alpha,alpha_scatter,beta,beta_scatter,Mhost)
    return median, min_val, max_val
 


######################################
#   dL_dm()
#   - a routine to calculate d (Lobs) / d (mbh) and d (Lobs) / d (mstar)
#   - this is used in solver2_LT()
######################################

def dL_dm(mbh,mstar):
    aa = 0.620
    bb = 0.674
    cc = 0.212
    dd = 0.553
    a = 1.27
    b = 0.300
    term_exp = np.exp((mstar - bb)/cc)
    xi_bh = a - b *(mbh)**0.242
    xi_star = (aa +term_exp)/(1 + dd * term_exp)
    
    f_bh = xi_bh**2.5 * mbh**(-1.0/6.0)
    f_star = xi_star**2.5 * mstar**(4.0/9.0)
    
    dxi_bh = -0.605 * b * xi_bh**1.5 /mbh**0.758
    dxi_star = 2.5*(1.0 - aa * dd) * term_exp * xi_star**1.5 / cc / (1.0 + dd * term_exp)**2.0
    df_dbh = dxi_bh * mbh**(-1.0/6.0) - 1.0/6.0 * xi_bh**2.5 * (mbh)**(-7.0/6.0)
    df_dstar = dxi_star * mstar**(4.0/9.0) + xi_star**2.5 * (4.0/9.0)* mstar**(-5.0/9.0)
    
    
    return f_bh, df_dbh, f_star, df_dstar
    
######################################
#   dT_dm()
#   - a routine to calculate d (Tobs) / d (mbh) and d (Tobs) / d (mstar)
#   - this is used in solver2_LT()
######################################

def dT_dm(mbh,mstar):
    
    aa = 0.620
    bb = 0.674
    cc = 0.212
    dd = 0.553
    a = 1.27
    b = 0.300
    term_exp = np.exp((mstar - bb)/cc)
    xi_bh = a - b* (mbh)**0.242
    xi_star = (aa +term_exp)/(1 + dd * term_exp)
    
    f_bh = xi_bh**(9.0/8.0) * mbh**(-3.0/8.0)
    f_star = xi_star**(9.0/8.0)
    
    dxi_bh = -0.27225 * b * xi_bh**0.125 /mbh**0.758  #D[xi_bh^(9.0/8.0),mbh]
    dxi_star = 1.125*(1.0 - aa * dd) * term_exp * xi_star**0.125 / cc / (1.0 + dd * term_exp)**2.0 #D[xi_star^(9.0/8.0),mstar]
    
    df_dbh = dxi_bh * mbh**(-3.0/8.0) - 3.0/8.0 * xi_bh**(9.0/8.0) * (mbh)**(-11.0/8.0)
    df_dstar = dxi_star
    return f_bh, df_dbh, f_star, df_dstar
    
######################################
#   get_jacobian_LT()
#   - a routine to calculate jacobian matrix used in solver2_LT()
######################################

def get_jacobian_LT(mbh,mstar,rstar,c1,x):
    jac = np.zeros((2,2))
    tpeak_factor = 3.0/2.0
    
    L_bh, dL_dbh, L_mstar, dL_dmstar = dL_dm(mbh,mstar)
    T_bh, dT_dbh, T_mstar, dT_dmstar = dT_dm(mbh,mstar)
    factor1 = 0.8571575661411786 * (0.5/c1) #* x **(-5.0/3.0)
    jac[0][0] = factor1 * L_mstar * dL_dbh
    jac[0][1] = factor1 * L_bh * dL_dmstar
    factor2 = 3.274093669346937 * (c1/0.5)**(-3.0/4.0)#* x**(-5.0/12.0)
    jac[1][0] = factor2 * T_mstar * dT_dbh
    jac[1][1] = factor2 * T_bh * dT_dmstar
    
    return jac
    
 
######################################
#   find_mbh_mstar_from_input()
#   - a routine to find ranges of M_BH and M_star which reproduces Lobs and Tobs
#   - the central value of the ranges is used for a first guess in solver1_LT()
#   - If solver1_LT and solver2_LT fail, the central value and the uncertainties found in this routine are provided
######################################


def find_mbh_mstar_from_input(input_array,sample,mbh_mstar_array,N_sampling,mbh,mstar,c1,del_omega,input_count,mass_range_array,input_kind,double_intersection):
    LPEAK = 0
    TPEAK = 1
    RBB = 2
    LINEWIDTH = 3
    MBH = 4
    found = 0
    min_error = 0.00 #in percentile
    min_range = 0.05

   
    A_obs = input_array[0][sample]
    error1 = input_array[1][sample][0]
    error2 = input_array[1][sample][1]
    A_obs_range =[A_obs - error1, A_obs + error2]

    nan = "OK"
    if(A_obs!=0.0):
        found = 0
        peak_value = np.zeros((N_sampling,N_sampling))
        for i_bh in range(N_sampling):
            for i_star in range(N_sampling):
                bh = mbh[i_bh]
                star = mstar[i_star]
                rstar = get_rstar(star)
                t_obs = get_tpeak(bh,star,rstar)
                if(input_kind==LPEAK):
                    peak_value[i_star,i_bh] = get_Lobs(bh,star,rstar,c1,t_obs)
                elif(input_kind ==TPEAK):
                    peak_value[i_star,i_bh] = get_Tobs(bh,star,rstar,c1,del_omega,t_obs)
#                elif(input_kind==LINEWIDTH):
#                    peak_value[i_star,i_bh] = get_v(bh,star,rstar,c1)/1e5
#                elif(input_kind==MBH):
#                    peak_value[i_star,i_bh] = bh*1e6
                if(np.isfinite(peak_value[i_star,i_bh])==False):
                    if(np.isnan(peak_value[i_star,i_bh])==True):
                        nan = "NaN"
                    else:
                        nan = "Inf"
                if((A_obs_range[1] - A_obs_range[0]>0.0 and peak_value[i_star,i_bh] > A_obs_range[0] and peak_value[i_star,i_bh] < A_obs_range[1]) or  np.abs(peak_value[i_star,i_bh]-A_obs_range[0])/A_obs_range[0] < min_range):
                        
                        mbh_mstar_array[input_count][0].append(bh)
                        mbh_mstar_array[input_count][1].append(star)
                        mass_range_array[input_count,i_star,i_bh] = 1.0
                        found = 1
#                        overall_intersection[i_star,i_bh] = overall_intersection[i_star,i_bh] + 1
                        if(input_kind==LPEAK or input_kind==TPEAK):
                            double_intersection[i_star,i_bh] = double_intersection[i_star,i_bh] + 1
        mbh_mstar_array[input_count][0] = np.array(mbh_mstar_array[input_count][0])
        mbh_mstar_array[input_count][1] = np.array(mbh_mstar_array[input_count][1])
        if(found==0):
            print ("No inferred mass within the given search range (0: Lpeak, 1: Tpeak)", input_kind)
    else:
        print ("INPUT VALUE IS ZERO")


    return double_intersection,mbh_mstar_array,mass_range_array,nan



######################################
#   find_centroid_range()
#   - find a central value of the ranges found using find_mbh_mstar_from_input()
######################################

    
def find_centroid_range(N_sampling,mbh,mstar,double_intersection):
    accu_bh = 0.0
    accu_star = 0.0
    area_total = 0.0
    bh_array=[]
    star_array=[]
    int_count = 2
    for i_bh in range(N_sampling):
        for i_star in range(N_sampling):
            bh = mbh[i_bh]
            star = mstar[i_star]
            rs = get_rstar(star)

            value = double_intersection[i_star,i_bh]
            if(value ==int_count):
                if(i_star==0):
                    area_star = mstar[i_star+1]-mstar[i_star]
                elif(i_star == N_sampling-1):
                    area_star = mstar[i_star]-mstar[i_star-1]
                else:
                    area_star = (mstar[i_star+1]-mstar[i_star-1])/2.0
                if(i_bh==0):
                    area_bh = mbh[i_bh+1]-mbh[i_bh]
                elif(i_bh == N_sampling-1):
                    area_bh = mbh[i_bh]-mbh[i_bh-1]
                else:
                    area_bh = (mbh[i_bh+1]-mbh[i_bh-1])/2.0
                
                area = area_star * area_bh
                accu_bh = accu_bh + mbh[i_bh]*area
                accu_star = accu_star + mstar[i_star]*area
                area_total = area_total + area
                bh_array.append(mbh[i_bh])
                star_array.append(mstar[i_star])
    centroid_bh = accu_bh/ area_total
    centroid_star = accu_star/ area_total
    max_bh = np.amax(bh_array)
    min_bh = np.amin(bh_array)
    max_star = np.amax(star_array)
    min_star = np.amin(star_array)
    return centroid_bh,min_bh,max_bh,centroid_star,min_star,max_star



def relative_error_calc(param1, param2,val_ref1,val_ref2, mbh_sol,mstar_sol,c1,del_omega):
    param = [param1,param2]
    val_ref = [val_ref1,val_ref2]
    LPEAK = 0
    TPEAK = 1

    error = [[],[]]
    for i in range(len(param)):
        if(param[i]==LPEAK):
            val = get_Lobs(mbh_sol,mstar_sol,get_rstar(mstar_sol),c1,2.0)
            error[LPEAK] = np.abs((val_ref[LPEAK]-val)/val_ref[LPEAK])
        elif(param[i]==TPEAK):
            val = get_Tobs(mbh_sol,mstar_sol,get_rstar(mstar_sol),c1,del_omega,2.0)
            error[TPEAK] = np.abs((val_ref[TPEAK]-val)/val_ref[TPEAK])

    return error[LPEAK],error[TPEAK]
