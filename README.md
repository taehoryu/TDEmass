# TDEmass


1. Description

    1) Model - slow circularization
    
        TDEmass is a tool for interpretation of Tidal Disruption Event (TDE) observations.   In TDEs, a supermassive black hole 
        at the center of a galaxy tears apart an ordinary star; the debris is placed on highly eccentric orbits and, in ways that remain 
        controversial, ultimately produces a very bright flare.   The spectrum of the optical/UV light in this flare is generally well 
        described by a black body with temperature ~ few x 10^4 K.

        Using this tool, one can infer the mass of the black hole (mbh) and the mass of the star (mstar) involved in a TDE 
        by solving two non-linear algebraic equations derived in Ryu+2020 (https://ui.adsabs.harvard.edu/abs/2020arXiv200713765R/abstract) 
        (Equations 9 and 10). These equations arise from a physical model for the optical/UV luminosity in which the energy 
        for this light is generated by dissipation in shocks near the debris' orbital apocenter. They may be written in the following form:

        mbh = 0.5 (Tobs/10^4)^{-8/3} (del_omega/2)^{-2/3} c_1^{-2} Xi(mbh,mstar)^3
        mstar = 5 (Lobs/10^{44})^{9/4} (Tobs/10^{4.5})^{-1} (del_omega/2)^{-1/4} c_1^{3/2} \Xi(mbh,msar)^{-9/2}.

        where Lobs refers to the observed peak UV/optical luminosity (in units of 10^44 erg/s) and Tobs are the observed temperature 
        at peak luminosity (in units of 10^4.5 K). mbh are mstar are in units of 10^6 Msol and 1 Msol, respectively. c_1 and del_omega are
        two unspecified parameters (see Section 2) below ). 

        \Xi(mbh,mstar) is the ratio of the width of the debris' specific energy distribution to its conventional estimate G mbh R_star/ r_t^{2} 
        where R_star is the stellar radius and r_t is the order of magnitude estimate for the tidal radius, i.e., (mbh/mstar)^(1/3) R_star 
        (Ryu+2020, Arxiv 2001:03501).  \Xi(mbh,mstar) is nonlinear, but analytic function of mbh and mstar.  Its nonlinearity is 
        what requires numerical solutions of these equations.

        The code has been updated recently to incorporate the cooling effect discussed in Krolik, Piran, Ryu (2024) arXiv:2409.02894.
        Two more parameters have been added (aspect ratio - h_r and the fraction of the bound mass within the semimajor axis of the most tightly
        bound debris at peak mass return time). If include_tcool_tdyn_ratio = 1, the code will instead solve for Equations 17 and 19 in Krolik, Piran, Ryu (2024). 
   
   3) Two unspecified parameters c_1 and Del_omega
   
        Our model includes two unspecified parameters, c_1 and Del_omega. c_1*a_0 is the distance from the black hole at which a significant 
        amount of energy is dissipated by shocks. Here, a_0 is the apocenter distance for the orbit of the most tightly bound debris. 
        Del_omega is the solid angle. So Del_omega * (c_1 * a_0)^2 corresponds to the surface area of the emitting region. 
        Note that c1a0 (Del Omega/ 4pi)^(1/2) is equivalent to what is called the "black body radius". For fixed Lobs and Tobs, 
        larger Del_omega leads to smaller mbh and mstar, but the dependence is weak. However, mbh and mstar are sensitive to c_1. 
        For fixed Lobs and Tobs, mbh \propto c_1^{a} with a = -(1.2-2) and mstar \propto c_1^{b} with b = 0.8-1.5. 

        The quantities c1 and del_omega are coefficients poorly-determined by current theory; we recommend setting 
        c1=1 and del_omega=2 (default values in the code).
  
  
    4) Code description
    
        The code solves the equations using the following algorithm:

        - Create a table of Lobs and Tobs values on a logarithmic grid in mbh and mstar within ranges specified by mbh_range and mstar_range. 
          There are 500 within each mass range.

        - Locate Lobs and Tobs within this table and define the uncertainty in mbh and mstar by the uncertainty in Lobs and Tobs.   
          These uncertainties are called dLobs-, dLobs+, dTobs-, and dTobs+.

        - Find the central values of mbh and mstar within the ranges determined by the uncertainty in mbh and mstar 
          (see find_centroid_range() in module.py for the definition of the central value)

        - With the first guess, find new values of mbh and mstar by using the two basic equations. Iterate to convergence. 
          The initial values of mbh and mstar are chosen to be the central values of mbh and mstar. The routine embodying this step is solver1_LT().

        - If the previous step fails to converge, use the 2-dimensional Newton-Raphson method to solve the equations.  
          The initial values of mbh and mstar  are the central values of mbh and mstar.  This routine is called solver2_LT().

    


        The units in the code are:
        - Luminosity(Lobs, Lobs_1, Lobs_2) = erg/s
        - Temperature (Tobs, Tobs_1, Tobs_2) = K
        - black hole mass (mbh) = 10^6msol (msol : solar mass)
        - stellar mass (mstar) = msol


1. Download and run

	To download, clone the git repository using HTTP access :

		https://github.com/taehoryu/TDE_mass_inference.git

	To run, you first enter the directory "TDE_mass_inference". In the directory, the code is run as

		Python3 main.py

	The code has three source dependencies : scipy, numpy and matplotlib


2. Inputs

    Before running, two input files must be placed in "input_file" directory:  
 
    1) model_info.txt : (This must be in plain text)
        - inputdata_file_name : the name of the input observational data file.
        
        - output_file_name : the name of the output file with the inferred black hole and stellar mass. If "-" is given, the output file name is given 
          with the case names: If there is one case in the input file, output_file_name = "case-name.out". If there are 
          more than two, output_file_name = "(first-case-name)_(last-case-name).out" where "first-case-name" ("last-case-name") 
          is the name of the first (last) case listed in the input data file.
        
        - c1 : the ratio between the lengthscale of the emitting region and the apocenter of the characteristic debris orbit; we recommend setting this to 1.
    
        - Del_omega : Solid angle (in units of pi) of radiation from the emission region;we recommend setting this number to 2.

        - mbh_range : the range of black hole mass, in units of 10^{6} solar mass, within which the solver tries to find the solution. 

        - mstar_range : These two define the range of stellar mass, in units of solar mass, within which the solver
                             tries to find the solution.

        - include_tcool_tdyn_ratio : integer value to decide whether the luminosity is limited by the Eddington luminosity if the cooling time is shorter than the dynamical time (1: yes, 2: no)

        - f                        : the fraction of the bound mass within the semimajor axis of the most tightly bound debris at peak mass return time (default value = 0.15). This is only used when include_tcool_tdyn_ratio=1                

        - aspect_ratio             : the aspect ratio of the debris disk used to calculate the cooling time (default value = 0.5). This is only used when include_tcool_tdyn_ratio=1

        
	2) input data file (the name is given in "model_info.txt")

		This can be either in plain text or in .csv format. The input data file should have 7 columns:

			Case_name   Lobs    dLobs-    dLobs+    Tobs    dTobs-    dTobs+

		- Case_name : you can put any name here. This name will be used as the title of the solution figures (see 4. output).

		- Lobs  : the peak luminosity in units of erg/s. 
    
		- dLobs-, dLobs+: These two values are determined by the uncertainties in Lobs. (IMPORTANT) Lobs - (dLobs-) and Lobs + (dLobs+) are the lower and upper limits 
                       for Lobs, respectively. For example, if log10(Lobs)=43_{-0.1}^{+0.1}, Lobs = 10^43, dLobs- = 10^{43}-10^{43-0.1} = 2.1 \times 10^42 and 
                       dLobs+ = 10^{43+0.1} - 10^{43} = 2.6 \times 10^42. This definition applies to the uncertainties for all quantities produced in the output file.

		- Tobs : the black body temperature at the peak luminosity in units of K.

		- dTobs-, dTobs+ : Same for dLobs- and dLobs+. See the definitions for dLobs- and dLobs+.

		The input.csv file in "example" directory has two example cases. you can take copy and paste the file and put it in "input_file" directory.



        NOTE that examples of the input file and input data file can be found in example directory (and input_file directory).

4. Outputs

	One text output file and a number of .png plot files are created in "output" directory.

	1) text file with the inferred black hole mass for the candidates

		This output file consists of 19 columns. The first seven columns are identical to those in the input data file. 
                The rest of the columns are the inferred black hole mass[mbh], its uncertainties (dmbh-,dmbh+), the inferred stellar mass [mstar], 
                its uncertainties (dmstar-,dmstar+), t0 and its uncertainties (dt0-, dt0+) and a0 and its uncertanties (da0-,da0+). 
                Here, t0 is the characteristic mass return time of the most tightly bound debris. 
                The uncertainties for the masses, t0 and a0 are defined the same way as those for Lobs (see Section 3-2).
                The name of the output file is given in the "model_info.txt" (see Section 3-1)

	2) solution figures on the mbh - mstar plane for each candidate

		The plots are created as many as the number of the candidates given in the input data file. There are two strips as long as the solutions 
               are found within the ranges of mbh and mstar set in the input data file. The blue strip demarcates the solutions for the given Lobs and the red strip 
               for the given Tobs. The X-hatched green region indicates the solutions for the given Lobs and Tobs. 

	3) One solution figure on the mbh - mstar plane for all cases
                This figure shows the inferred mbh and mstar found within the ranges of given mbh and mstar for the cases listed in the input data file. 
		The name of the figure is  "c1_(c1_value)_del_omega_(del Omega value)_inferred_mass.png"
		"c1_value" and "del Omega value" are the values of c_1 and Del_omega, respectively.
		The error bars are determined to enclose the extreme values of the solutions. 
        
        
        NOTE that the outputs for the example case in the example input file and input data file can be found in example directory.

5. Issues/requests

    - E-mail address for issues/requests : udraeo@gmail.com or Taeho.Ryu@colorado.edu

6. Acknowledgement

If this software contributes to your project, please acknowledge it by citing,

@ARTICLE{2020arXiv200713765R,
    author = {{Ryu}, Taeho and {Krolik}, Julian and {Piran}, Tsvi},
    title = "{Measuring stellar and black hole masses of tidal disruption events}",
    journal = {arXiv e-prints},
    keywords = {Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Astrophysics of Galaxies, Astrophysics - Solar and Stellar Astrophysics},
    year = 2020,
    month = jul,
    eid = {arXiv:2007.13765},
    pages = {arXiv:2007.13765},
    archivePrefix = {arXiv},
    eprint = {2007.13765},
    primaryClass = {astro-ph.HE},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv200713765R},
    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
