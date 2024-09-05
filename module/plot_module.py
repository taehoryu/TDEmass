import numpy as np
import glob
import os
import os.path
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('Agg')
import matplotlib.colors as colors
import matplotlib.cm
#import matplotlib.cm as cm
#from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
import errno
import math
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import itertools
from matplotlib import rcParams
import module

def plotting(index_label, figure_storage, plot_array, mbh, mstar, mass_range_array, mbh_mstar_array, mean_bh, mean_star, plot_format, plot_quality, include_tcool, hr, f, c1, del_omega):
    matplotlib.rcParams.update({'font.size': 20})
    marker = itertools.cycle(('<', 'X', '^', 'o', '*','>','p','s'))
    colors=["blue", "red", "green", "magenta","orange","grey"]

    fig, ax = plt.subplots()
    fig.canvas.draw()
    
    index_label2=""
    for letter in index_label:
        if(letter=="/"):
            index_label2 = index_label2 + "_"
        else:
            index_label2 = index_label2 + letter
    
    if(include_tcool==1):
        index_label2 += "_with_cooling_hr_%3.3f"%(hr)
    filename = index_label2 + plot_format
    filename1 = os.path.join(figure_storage, filename)
    XBH, YSTAR = np.meshgrid(mbh, mstar)
    
    for count in range(2):
        if(len(mbh_mstar_array[count][0])>0):
            ax.contourf(XBH,YSTAR,mass_range_array[count],colors = colors[count],alpha=0.3)
    ax.scatter(mean_bh, mean_star,zorder=99, color="magenta", marker= "X", s= 200)
    plt.rcParams['hatch.color'] = "white"
    masked = np.ma.masked_where(plot_array!=2,plot_array)
    ax.contourf(XBH,YSTAR,masked,colors = "green",alpha=1.0,hatches=['x'])

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(mstar[0],mstar[-1])
    ax.set_xlim(mbh[0],mbh[-1])
      
    xtick_array= get_tick_array(mbh[0],mbh[-1],"x")
    ytick_array= get_tick_array(mstar[0],mstar[-1],"y")
    
    ax.text(0.65,0.23, r"$c_{1}$            = %4.1f"%(c1), transform = ax.transAxes, fontsize = 12)
    ax.text(0.65,0.18, r"$\Delta \Omega$           = %4.1f $\pi$"%(del_omega/np.pi), transform = ax.transAxes, fontsize = 12)

    if(include_tcool==1):
        ax.text(0.65,0.13, "Cooling?  = Yes", transform = ax.transAxes, fontsize = 12)
        ax.text(0.65,0.08, "$f$              = %4.2f"%(f), transform = ax.transAxes, fontsize = 12)
        ax.text(0.65,0.03, "$h/r$           = %4.2f"%(hr), transform = ax.transAxes, fontsize = 12)
    else:
        ax.text(0.65,0.13, "Cooling?  = No", transform = ax.transAxes, fontsize = 12)

    ax.set_xticks(xtick_array)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))

    ax.yaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(ytick_array)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func))
      
    ax.yaxis.set_label_coords(-0.1,0.5)
    ax.set_ylabel(r"$M_{\star}\/[\mathrm{M}_{\odot}]$")
    ax.set_xlabel(r"$M_{\mathrm{BH}}\/[10^{6}\/\mathrm{M}_{\odot}]$")
    
    ax.set_title(index_label, fontsize = 25)
    ax = plt.gca()
    plt.tight_layout()
    plt.gcf().subplots_adjust(left=0.16)
    plt.gcf().subplots_adjust(bottom=0.16)

    plt.savefig(filename1,dpi=plot_quality)
    plt.close()
    return



def plot_double_intersection(figure_storage, index_array, plot_array, mass_range_array, mbh_search_range, mstar_search_range, plot_format,plot_quality,c1,del_omega,samplesize,mbh_sol_array,mstar_sol_array,solution_exist, include_tcool, hr):
    marker = itertools.cycle(('<', 'X', '^', 'o', '*','>','p','s'))
    colors=["blue", "red", "green", "magenta","orange","grey"]
    matplotlib.rcParams.update({'font.size': 20})
    temp_array=[[[],[]],[[],[]]]
    for sample in range(samplesize):
        if(solution_exist[sample]==1):
            mbh_lo = mbh_sol_array[sample][0] - mbh_sol_array[sample][1]
            mbh_hi = mbh_sol_array[sample][2] + mbh_sol_array[sample][0]
            mstar_lo = mstar_sol_array[sample][0] - mstar_sol_array[sample][1]
            mstar_hi = mstar_sol_array[sample][2] + mstar_sol_array[sample][0]
            temp_array[0][0].append(mbh_lo)
            temp_array[0][1].append(mbh_hi)
            temp_array[1][0].append(mstar_lo)
            temp_array[1][1].append(mstar_hi)
        
    mbh = 10.0** np.linspace(np.log10(np.amin(temp_array[0][0]))*0.75,np.log10(np.amax(temp_array[0][1]))*1.25,100)
    mstar = 10.0** np.linspace(np.log10(np.amin(temp_array[1][0]))*0.75,np.log10(np.amax(temp_array[1][1]))*1.25,100)
        
    XBH, YSTAR = np.meshgrid(mbh, mstar)
    
    
    
    
    
    cmap = "jet"
    fig, ax = plt.subplots()
    int_count = 2
    filename = "c1_"+str(c1)[0:5]+"_del_omega_"+str(del_omega/math.pi)[0:3]+"_inferred_mass"
    if(include_tcool==1):
        filename +="_with_cooling_hr_%3.3f"%(hr)
    filename += plot_format
    filename1 = os.path.join(figure_storage, filename)
    good_sample_size = 0
    offset=[]
    for sample in range(samplesize):
        maxval = np.amax(plot_array[sample].flatten())
        if(maxval==int_count):
            good_sample_size = good_sample_size + 1
        offset.append(good_sample_size)
    
    vmin = 0
    vmax = max(good_sample_size,1)
    sample_count = 0
    min_bh = 1e100
    min_star = 1e100
    max_bh = -1e100
    max_star = -1e100
    for sample in range(samplesize):
      if(solution_exist[sample]==1):
        mean_bh = 0.0
        mean_star = 0.0
        area_total = 0.0

        bh_sol = mbh_sol_array[sample][0]
        error_bh_l = mbh_sol_array[sample][1]
        error_bh_h = mbh_sol_array[sample][2]
        star_sol = mstar_sol_array[sample][0]
        error_star_l = mstar_sol_array[sample][1]
        error_star_h = mstar_sol_array[sample][2]
        mean_bh = mbh_sol_array[sample][0]
        ccmap = plt.get_cmap(cmap)
        color = ccmap((offset[sample]-vmin)/(vmax-vmin))
        if (mean_bh>0.0):
            markers = next(marker)
            pp = ax.errorbar(bh_sol,star_sol,yerr=np.array([(error_star_l, error_star_h)]).T,xerr = np.array([(error_bh_l, error_bh_h)]).T,color=color, fmt='--',capsize=2,elinewidth= 1)
            pp = ax.scatter(bh_sol,star_sol,marker = markers,color=color,label=index_array[sample])
        
            min_bh = min(bh_sol - mbh_sol_array[sample][1], min_bh)
            max_bh = max(bh_sol + mbh_sol_array[sample][2], max_bh)
            min_star = min(star_sol - mstar_sol_array[sample][1], min_star)
            max_star = max(star_sol + mstar_sol_array[sample][2], max_star)
    ax.legend(loc=2,fontsize=5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    
#    ax.set_ylim(min_star*0.8,max_star*1.2)
#    ax.set_xlim(min_bh*0.8,max_bh*1.2)
    ax.set_xlim(mbh_search_range[0],mbh_search_range[1])
    ax.set_ylim(mstar_search_range[0],mstar_search_range[1])
    xtick_array= get_tick_array(min_star,max_star,"x")
    ytick_array= get_tick_array(min_bh,max_bh,"y")

    ax.xaxis.set_major_formatter(NullFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())

    ax.set_xticks(xtick_array)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(ytick_array)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func))
    
    ax.yaxis.set_label_coords(-0.11,0.5)
    ax.set_ylabel(r"$M_{\star}\/[\mathrm{M}_{\odot}]$")
    ax.set_xlabel(r"$M_{\mathrm{BH}}\/[10^{6}\/\mathrm{M}_{\odot}]$")
    ax = plt.gca()
    plt.tight_layout()
    plt.gcf().subplots_adjust(left=0.16)
    plt.gcf().subplots_adjust(bottom=0.16)

    plt.savefig(filename1,dpi=plot_quality)
    plt.close()

    return


def get_tick_array(min,max,axis):
  xarray=[0.1,0.5,1,3,5,10,50]
  yarray=[0.15,0.3,1,3,10]
  tick_array =[]
  if(axis=="x"):
    for i in xarray:
      if(i>=min and i<=max):
        tick_array.append(i)
  if(axis=="y"):
    for i in yarray:
      if(i>=min and i<=max):
        tick_array.append(i)

  return tick_array



def format_func(value, tick_number):
  # find number of multiples of pi/2

  if value ==0.01:
    return "0.01"
  elif value ==0.1:
    return "0.1"
  elif value ==0.15:
    return "0.15"
  elif value ==0.3:
    return "0.3"
  elif value == 0.5:
    return "0.5"
  elif value == 1.0:
    return "1"
  elif value == 2.0:
    return "2"
  elif value == 4.0:
    return "4"
  elif value == 6.0:
    return "6"
  elif value == 7.0:
    return "7"
  elif value == 8.0:
    return "8"
  elif value == 9.0:
    return "9"

  elif value == 0.8:
    return "0.8"#r"$6\times10^{-3}$"

  elif value == 3.0:
    return "3"#r"$4\times10^{-2}$"
  elif value == 5.0:
    return "5"#r"$4\times10^{-2}$"
  elif value == 10.0:
    return "10"#r"$3\times10^{-1}$"
  elif value == 50.0:
    return "50"#r"$3\times10^{-1}$"
  elif value == 100.0:
    return "100"#r"$3\times10^{-1}$"
