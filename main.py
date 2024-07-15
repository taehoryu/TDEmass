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
figure_storage = "output"
module.make_sure_path_exists(figure_storage)
LPEAK = 0
TPEAK = 1
RBB = 2
LINEWIDTH = 3
MBH = 4
plot_format = ".png"
plot_quality = 200
TOL = 1e-4
N_sampling = 500
tot_input_variable_count = 4
check_input = [0,0,0,0,0,0,0]
####################################################



#Read key input parameter, input file name and output file name
inputdata_file_name, output_file_name, c1, del_omega,mstar_search_range,mbh_search_range, include_tcool, h_r  = module.read_model_input()


#Read input data
index_array, mstar_range, mbh_range, Lpeak_array, Tpeak_array, samplesize = module.read_input_data(inputdata_file_name, mstar_search_range, mbh_search_range)
print ('{:^20} {:^15}'.format("    inputdata_file_name = ",inputdata_file_name))
print ('{:^20} {:^15}'.format("    outputfile name     = ",output_file_name))
print ('{:^20} {:.2f}'.format("                     c1 = ", c1))
print ('{:^20} {:.2f}'.format("        Delta omega[pi] = ", del_omega/math.pi))
print ('{:^20} {:d}'.format("      Total sample size = ", samplesize))
if( include_tcool==1):
    print ('{:^20}'.format("Include cooling effect? = yes"))
    print ("           Aspect ratio = %5.3e"%(h_r))

else:
    print ('{:^20}'.format("Include cooling effect? = no"))

double_intersection_sample = np.zeros(samplesize)
if(output_file_name==""):
    if (samplesize == 1):
        output_file_name = index_array[0] +".txt"
    else:
        output_file_name = index_array[0]+ "_"+ output_file_name[-1]+".txt"
    
#open the output file and write a header line.
output_file=open(figure_storage+"/"+output_file_name,"w")
output_file.write("{:^25}".format("candidate_name")+"{:^11}".format("Lobs")+"{:^9}".format("dLobs-")+"{:^10}".format("dLobs+")+
"{:^12}".format("Tobs")+"{:^7}".format("dTobs-")+"{:^11}".format("dTobs+")+"{:^10}".format("mbh")+
"{:^11}".format("dmbh-")+"{:^12}".format("dmbh+")+"{:^8}".format("mstar")+"{:^8}".format("dmstar-")+
"{:^9}".format("dmstar+")+"{:^8}".format("t0")+"{:^8}".format("dt0-")+"{:^8}".format("dt0+")+
"{:^9}".format("a0")+"{:^11}".format("da0-")+"{:^9}".format("da0+")+"\n")

output_file.write("{:^25}".format("          ")+"{:^11}".format("[erg/s]")+"{:^10}".format("[erg/s]")+"{:^9}".format("[erg/s]")+
"{:^12}".format("[K]")+"{:^8}".format("[K]")+"{:^9}".format("[K]")+"{:^12}".format("[10^6msol]")+
"{:^10}".format("[10^6msol]")+"{:^13}".format("[10^6msol]")+"{:^6}".format("[msol]")+"{:^10}".format("[msol]")+
"{:^9}".format("[msol]")+"{:^7}".format("[days]")+"{:^8}".format("[days]")+"{:^8}".format("[days]")+
"{:^10}".format("[10^14cm]")+"{:^10}".format("[10^14cm]")+"{:^10}".format("[10^14cm]")+"\n")


double_intersection_array=[]
mbh_sol_array=[]
mstar_sol_array=[]
t0_sol_array =[]
solution_exist=[]
for sample in range(samplesize):
    double_intersection = np.zeros((N_sampling, N_sampling), dtype = int)
    
    mbh_min = mbh_range[sample][0]#/1e6
    mbh_max = mbh_range[sample][1]#/1e6
    mstar_min = mstar_range[sample][0]
    mstar_max = mstar_range[sample][1]
    mbh = 10**np.linspace(np.log10(mbh_min), np.log10(mbh_max), N_sampling)
    mstar = 10**np.linspace(np.log10(mstar_min), np.log10(mstar_max), N_sampling)
    check_input = [0, 0, 0, 0, 0, 0, 0]
    
    nan_check=[]
    mbh_mstar_array =[]
    mass_range_array= np.empty((tot_input_variable_count, N_sampling, N_sampling))
    input_variable_count = -1
    for i in range(tot_input_variable_count):
        mass_range_array[i, :, :] = np.NaN
        mbh_mstar_array.append([])
        mbh_mstar_array[i].append([])
        mbh_mstar_array[i].append([])

    SOL_FOUND=0
    #First find ranges of mbh and mstar which reproduces Lobs and Tobs within their uncertainties
    input_variable_count += 1
    double_intersection, mbh_mstar_array, mass_range_array, nan, t0_range=module.find_mbh_mstar_from_input(Lpeak_array, sample,
    mbh_mstar_array, N_sampling, mbh, mstar, c1, del_omega, input_variable_count, mass_range_array, LPEAK, double_intersection, include_tcool, h_r)
    check_input[input_variable_count] = 1
    nan_check.append(nan)

    input_variable_count += 1
    double_intersection, mbh_mstar_array, mass_range_array, nan, t0_range=module.find_mbh_mstar_from_input(Tpeak_array, sample,
    mbh_mstar_array, N_sampling, mbh, mstar, c1, del_omega, input_variable_count, mass_range_array, TPEAK, double_intersection, include_tcool, h_r)
    check_input[input_variable_count] = 1
    nan_check.append(nan)


    #First find centroid of the ranges found and use it for the first guess in solver1_LT()
    centroid_bh, min_bh, max_bh, centroid_star, min_star, max_star, retv_centroid = module.find_centroid_range(N_sampling, mbh, mstar, double_intersection)
    print ("retv_centroid",retv_centroid, centroid_bh,centroid_star)
    if(retv_centroid == 0):
        #First try to find the solution
        retv, mbh_sol, mstar_sol = module.solver1_LT(Lpeak_array[0][sample], Tpeak_array[0][sample], centroid_bh, centroid_star, c1, del_omega)
        #Check if the solutions are correct. If not try the second solver
        error_L, error_T = module.relative_error_calc(LPEAK, TPEAK, Lpeak_array[0][sample], Tpeak_array[0][sample],  mbh_sol, mstar_sol, c1, del_omega, include_tcool, h_r)
        print ("Error 1", error_L, error_T)
        #if include_tcool=1, the final solution is found in solver2_LT()
        if(retv!=0 or error_L > TOL or error_T > TOL or include_tcool==1):
            retv, mbh_sol, mstar_sol = module.solver2_LT (Lpeak_array[0][sample], Tpeak_array[0][sample], centroid_bh, centroid_star, c1, del_omega, include_tcool, h_r)
            error_L, error_T = module.relative_error_calc(LPEAK, TPEAK, Lpeak_array[0][sample], Tpeak_array[0][sample],  mbh_sol, mstar_sol, c1, del_omega, include_tcool, h_r)
            print ("ERROR", error_L, error_T)
            if(include_tcool==1):
                mbh_sol_range = []
                mstar_sol_range = []
                retv, mbh_sol1, mstar_sol1 = module.solver2_LT (Lpeak_array[0][sample] - Lpeak_array[1][sample][0], Tpeak_array[0][sample] - Tpeak_array[1][sample][0], mbh_sol, mstar_sol, c1, del_omega, include_tcool, h_r)
                if(retv==0):
                    mbh_sol_range.append(mbh_sol1)
                    mstar_sol_range.append(mstar_sol1)
                retv, mbh_sol1, mstar_sol1 = module.solver2_LT (Lpeak_array[0][sample] - Lpeak_array[1][sample][0], Tpeak_array[0][sample] + Tpeak_array[1][sample][1], mbh_sol, mstar_sol, c1, del_omega, include_tcool, h_r)
                if(retv==0):
                    mbh_sol_range.append(mbh_sol1)
                    mstar_sol_range.append(mstar_sol1)
                retv, mbh_sol1, mstar_sol1 = module.solver2_LT (Lpeak_array[0][sample] + Lpeak_array[1][sample][1], Tpeak_array[0][sample] - Tpeak_array[1][sample][0], mbh_sol, mstar_sol, c1, del_omega, include_tcool, h_r)
                if(retv==0):
                    mbh_sol_range.append(mbh_sol1)
                    mstar_sol_range.append(mstar_sol1)
                retv, mbh_sol1, mstar_sol1 = module.solver2_LT (Lpeak_array[0][sample] + Lpeak_array[1][sample][1], Tpeak_array[0][sample] + Tpeak_array[1][sample][1], mbh_sol, mstar_sol, c1, del_omega, include_tcool, h_r)
                if(retv==0):
                    mbh_sol_range.append(mbh_sol1)
                    mstar_sol_range.append(mstar_sol1)
                if(len(mbh_sol_range)>1):
                    min_bh = min(np.amin(mbh_sol_range),min_bh)
                    max_bh = max(np.amax(mbh_sol_range),max_bh)
                    min_star = min(np.amin(mstar_sol_range),min_star)
                    max_star = max(np.amax(mstar_sol_range),max_star)
        if(error_L<TOL and error_T<TOL):
            error_bh_l = mbh_sol - min_bh
            error_bh_h = max_bh - mbh_sol
            error_star_l = mstar_sol - min_star
            error_star_h = max_star - mstar_sol
            mbh_sol_array.append([mbh_sol, error_bh_l, error_bh_h])
            mstar_sol_array.append([mstar_sol, error_star_l, error_star_h])
            if(error_star_l < 0.0 or error_star_h < 0.0):
                print(index_array[sample], ": Either one of the masses or both are out of range. Change the range accordingly")
                print("Mstar = %10.5e"%(mstar_sol))
                print("Mbh   = %10.5e"%(mbh_sol))
            if(max_bh>=mbh_range[sample][1]):
                print(index_array[sample], ": [WARNING] The upper limit of Mbh reaches the range given in model_info.txt. It is possible that the upper limit is not well-determined.")
            if(min_bh<=mbh_range[sample][0]):
                print(index_array[sample], ": [WARNING] The lower limit of Mbh reaches the range given in model_info.txt. It is possible that the lower limit is not well-determined.")
            if(max_star>=mstar_range[sample][1]):
                print(index_array[sample], ": [WARNING] The upper limit of Mstar reaches the range given in model_info.txt. It is possible that the upper limit is not well-determined.")
            if(min_star<=mstar_range[sample][0]):
                print(index_array[sample], ": [WARNING] The lower limit of Mstar reaches the range given in model_info.txt.  It is possible that the lower limit is not well-determined.")

            SOL_FOUND=1
        else:
            error_bh_l = centroid_bh - min_bh
            error_bh_h = max_bh - centroid_bh
            error_star_l = centroid_star - min_star
            error_star_h = max_star - centroid_star
            mbh_sol_array.append([centroid_bh, error_bh_l, error_bh_h])
            mstar_sol_array.append([centroid_star, error_star_l, error_star_h])
            if(include_tcool==0):
                SOL_FOUND=1


    
    if(SOL_FOUND==1):
        a0_t0_sol = module.get_t0_a0_error(mbh_sol_array[sample], mstar_sol_array[sample], c1)
        print ("{:^25}".format(index_array[sample]), " [   Solution] m_bh[10^6msol]= {0:.2g}".format(mbh_sol),
        "_{-","{0:.2g}".format(error_bh_l),"}^{+","{0:.2g}".format(error_bh_h),"{:^25}".format("},    m_star[msol] =")
        ,"{0:.2g}".format(mstar_sol),"_{-","{0:.2g}".format(error_star_l),"}^{+","{0:.2g}".format(error_star_h),"}")
    else:
        print ("{:^25}".format(index_array[sample])," [No Solution] within the given mass range: mbh = [",mbh_range[sample][0],
        "-",mbh_range[sample][1],"] 10^{6}msol, mstar = [", mstar_range[sample][0],"-",mstar_range[sample][1],"] msol")
        mbh_sol_array.append([-100,100,100])
        mstar_sol_array.append([-100,100,100])
        a0_t0_sol=[[-100,-100,-100],[-100,-100,-100]]

    module.write_output(output_file, retv_centroid, index_array[sample], mbh_sol_array[sample], mstar_sol_array[sample],
    Lpeak_array[0][sample], Lpeak_array[1][sample], Tpeak_array[0][sample], Tpeak_array[1][sample], a0_t0_sol)
    solution_exist.append(SOL_FOUND)
    
    if(SOL_FOUND==1):
        #Plot the solutions on a (M_BH - M_star) grid
        print ("PLOTTING")
        pm.plotting(index_array[sample], figure_storage, double_intersection, mbh, mstar, mass_range_array, mbh_mstar_array, mbh_sol, mstar_sol, plot_format, plot_quality, include_tcool, h_r, c1, del_omega)
    double_intersection_array.append(double_intersection)
    
    
    


output_file.close()
#Plot the solutions and their uncertainties on a (M_BH - M_star) grid for entire sample
pm.plot_double_intersection(figure_storage, index_array, double_intersection_array, mass_range_array, plot_format,
plot_quality, c1, del_omega, samplesize, mbh_sol_array, mstar_sol_array, solution_exist,include_tcool, h_r)
