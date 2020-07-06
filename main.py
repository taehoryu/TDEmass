import numpy as np
import scipy as sp
import glob
import os
import os.path
import errno
import colorsys
import re
import sys
import math
sys.path.insert(0, 'module/')
import module
import constant
import plot_module as pm

fileDir = os.path.dirname(os.path.realpath('__file__'))

####################################################
#   Define a few parameters
####################################################
figure_storage = "output"
module.make_sure_path_exists(figure_storage)

LPEAK = 0
TPEAK = 1
RBB = 2
LINEWIDTH = 3
MBH = 4
plot_format = ".png"
plot_quality = 200
TOL = 1e-7
tot_input_variable_count = 4
check_input = [0,0,0,0,0,0,0]



#Read key input parameter and input file name
inputdata_file_name, output_file_name, c1, del_omega, N_sampling = module.read_model_input()


#Read input data
index_array,mstar_range,mbh_range, Lpeak_array,Tpeak_array,samplesize = module.read_input_data(inputdata_file_name)
print ('{:^15} {:^15}'.format("inputdata_file_name : ",inputdata_file_name))
print ('{:^15} {:^15}'.format("outputfile name     : ",output_file_name))
print ('{:^15} {:.2f}'.format("                 c1 = ", c1))
print ('{:^15} {:.2f}'.format("    Delta omega[pi] = ", del_omega/math.pi))
print ('{:^15} {:.2f}'.format("  Total sample size = ", samplesize))
double_intersection_sample = np.zeros(samplesize)

if(output_file_name==""):
    if (samplesize == 1):
        output_file_name = index_array[0] +".txt"
    else:
        output_file_name = index_array[0]+ "_"+ output_file_name[-1]+".txt"
    
output_file=open(figure_storage+"/"+output_file_name,"w")
output_file.write("    candidate_name           Lobs[erg/s]   Lobs_1[erg/s]   Lobs_2[erg/s]   Tobs[K]   Tobs_1[K]   Tobs_2[K]   mbh[10^6msol]   mbh_1[10^6msol]   mbh_2[10^6msol]   mstar[msol]   mstar_1[msol]      mstar_2[msol]  \n")


double_intersection_array=[]
mbh_sol_array=[]
mstar_sol_array=[]
for sample in range(samplesize):
    double_intersection = np.zeros((N_sampling,N_sampling),dtype = int)
    
    mbh_min = mbh_range[sample][0]#/1e6
    mbh_max = mbh_range[sample][1]#/1e6
    mstar_min = mstar_range[sample][0]
    mstar_max = mstar_range[sample][1]
    mbh = 10**np.linspace(np.log10(mbh_min),np.log10(mbh_max),N_sampling)
    mstar = 10**np.linspace(np.log10(mstar_min),np.log10(mstar_max),N_sampling)
    
    check_input = [0,0,0,0,0,0,0]
    
    nan_check=[]
    mbh_mstar_array =[]
    mass_range_array= np.empty((tot_input_variable_count,N_sampling,N_sampling))
    input_variable_count = -1
    for i in range(tot_input_variable_count):
        mass_range_array[i,:,:] = np.NaN
        mbh_mstar_array.append([])
        mbh_mstar_array[i].append([])
        mbh_mstar_array[i].append([])



    #First find ranges of mbh and mstar which reproduces Lobs and Tobs within their uncertainties
    input_variable_count += 1
    double_intersection,mbh_mstar_array,mass_range_array,nan=module.find_mbh_mstar_from_input(Lpeak_array,sample,mbh_mstar_array,N_sampling,mbh,mstar,c1,del_omega,input_variable_count,mass_range_array,LPEAK,double_intersection)
    check_input[input_variable_count] = 1
    nan_check.append(nan)


    input_variable_count += 1
    double_intersection,mbh_mstar_array,mass_range_array,nan=module.find_mbh_mstar_from_input(Tpeak_array,sample,mbh_mstar_array,N_sampling,mbh,mstar,    c1,del_omega,input_variable_count,mass_range_array,TPEAK,    double_intersection)
    check_input[input_variable_count] = 1
    nan_check.append(nan)


    #First find centroid of the ranges found and use it for the first guess in solver1_LT()
    centroid_bh,min_bh,max_bh,centroid_star,min_star,max_star = module.find_centroid_range(N_sampling,mbh,mstar,double_intersection)
    if(centroid_bh > 0.0):
        #First try to find the solution
        
        retv, mbh_sol,mstar_sol = module.solver1_LT(Lpeak_array[0][sample],Tpeak_array[0][sample],centroid_bh,centroid_star,c1,del_omega)
        L_check =module.get_Lobs(mbh_sol,mstar_sol,module.get_rstar(mstar_sol),c1,2.0)
        T_check =module.get_Tobs(mbh_sol,mstar_sol,module.get_rstar(mstar_sol),c1,del_omega,2.0)
        error_L =np.abs(L_check-Lpeak_array[0][sample])/Lpeak_array[0][sample]
        error_T =np.abs(T_check-Tpeak_array[0][sample])/Tpeak_array[0][sample]
        #Check if the solutions are correct. If not try the second solver
        if(retv!=0 or error_L > TOL or error_T > TOL):
            print ( index_array[sample]," Solver [1] fails -> Try Solver[2]")
            
            retv, mbh_sol, mstar_sol = module.solver2_LT (Lpeak_array[0][sample], Tpeak_array[0][sample], centroid_bh,centroid_star,c1,del_omega)
            L_check =module.get_Lobs(mbh_sol,mstar_sol,module.get_rstar(mstar_sol),c1,2.0)
            T_check =module.get_Tobs(mbh_sol,mstar_sol,module.get_rstar(mstar_sol),c1,del_omega,2.0)
            error_L =np.abs(L_check-Lpeak_array[0][sample])/Lpeak_array[0][sample]
            error_T =np.abs(T_check-Tpeak_array[0][sample])/Tpeak_array[0][sample]
            
        
        if(error_L<TOL and error_T<TOL):
            error_bh_l = mbh_sol - min_bh
            error_bh_h = max_bh - mbh_sol
            error_star_l = mstar_sol - min_star
            error_star_h = max_star - mstar_sol
            print (index_array[sample]," [SUCCESS] m_bh[10^6msol]= {0:.3g}".format(mbh_sol), ", m_star[msol] = {0:.3g}".format(mstar_sol))

            mbh_sol_array.append([mbh_sol,error_bh_l,error_bh_h])
            mstar_sol_array.append([mstar_sol,error_star_l,error_star_h])
        else:
            error_bh_l = centroid_bh - min_bh
            error_bh_h = max_bh - centroid_bh
            error_star_l = centroid_star - min_star
            error_star_h = max_star - centroid_star
            print (index_array[sample]," [Solver [2] fails] using centroid values : m_bh[10^6msol] = {0:.3g}".format(mbh_sol), ", m_star[msol] = {0:.3g}".format(mstar_sol))
            mbh_sol_array.append([centroid_bh,error_bh_l,error_bh_h])
            mstar_sol_array.append([centroid_star,error_star_l,error_star_h])
    else:
        print (index_array[sample]," [No Solution] within the given mass range")
        mbh_sol_array.append([-100,-100,-100])
        mstar_sol_array.append([-100,-100,-100])
           

    double_intersection_array.append(double_intersection)
    
    
    
    mean_bh = mbh_sol_array[sample][0]
    error_bh_l = mbh_sol_array[sample][1]
    error_bh_h = mbh_sol_array[sample][2]
    mean_star = mstar_sol_array[sample][0]
    error_star_l = mstar_sol_array[sample][1]
    error_star_h = mstar_sol_array[sample][2]
    
    lpeak = Lpeak_array[0][sample]
    l_error1 = Lpeak_array[1][sample][0]
    l_error2 = Lpeak_array[1][sample][1]
    TTpeak = Tpeak_array[0][sample]
    T_error1 = Tpeak_array[1][sample][0]
    T_error2 = Tpeak_array[1][sample][1]
    if(error_bh_l >= 1.0 and error_bh_h>=1.0):
        output_file.write("{:^25}     {:.3g}       {:.3g}         {:.3g}     {:.3g}     {:.3g}       {:.3g}          {:.3g}            {:.3g}              {:.3g}            {:.3g}           {:.3g}              {:.3g} \n".format(index_array[sample],lpeak,l_error1,l_error2,TTpeak,T_error1,T_error2,mean_bh,error_bh_l,error_bh_h,mean_star,error_star_l,error_star_h))
    else:
        output_file.write("{:^25}     {:.3g}       {:.3g}         {:.3g}     {:.3g}     {:.3g}       {:.3g}          {:.3g}            {:.3g}             {:.3g}           {:.3g}           {:.3g}              {:.3g} \n".format(index_array[sample],lpeak,l_error1,l_error2,TTpeak,T_error1,T_error2,mean_bh,error_bh_l,error_bh_h,mean_star,error_star_l,error_star_h))

    #Plot the solutions on a (M_BH - M_star) grid
    pm.plotting(index_array[sample],figure_storage, double_intersection, mbh, mstar, mass_range_array,mbh_mstar_array,mean_bh, mean_star,plot_format,plot_quality)

output_file.close()
#Plot the solutions and their uncertainties on a (M_BH - M_star) grid for entire sample
pm.plot_double_intersection(figure_storage,index_array,double_intersection_array,mass_range_array,plot_format,plot_quality,c1,del_omega,samplesize,mbh_sol_array,mstar_sol_array)
