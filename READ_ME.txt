1. Description

	This code is developed to infer the black hole mass (mbh) and stellar mass (mstar) from two input data (peak UV/optical luminosity, denoted by Lobs, and black body temperature at peak, denoted by T) by solving the two non-linear equations given in Ryu+2020 () (Equations 9 and 10). For each candidate, the code tries to solve the equations using three different methods.

1) It repeatedly computes Lobs and T from logarithmically sampled mbh and mstar (within the range defined in the input file, bin size = N_sampling) and find the ranges of mbh and mstar yielding Lobs and T within their uncertainties. If the following two methods fail, the central value for the ranges are given as the solutions.

2) solver1_LT() : With the first guess, it computes mbh and mstar using the equations. It uses the mean value of the first guess and the first computed values to get the second values. It repeats this procedure until they converge simultaneously.

3) solver2_LT() : 2-dimensional Newtonian-Raphson method.


2. Install and run

To install, clone the git repository using HTTP access :

https://github.com/taehoryu/TDE_mass_inference.git

To run, you first enter the directory "TDE_mass_inference". In the directory, the code is run as

Python3 main.py

The code has two source dependencies : scipy and numpy


3. Units

The default units are:

-Luminosity(Lobs, Lobs_1, Lobs_2) = erg/s
-Temperature (Tobs, Tobs_1, Tobs_2) = K
-black hole mass (mbh) = 10^6msol (msol : solar mass)
-stellar mass (mstar) = msol
 

3. Input data format

Before running, you should prepare two input files in "input_file" directory:  one that contains the value of key parameters and  the other that contains the spectral data for the events that you are interested in. 
 

	1) model_info.txt : define key parameters. 
	
		- inputdata_file_name : the name of the input data file containing the spectral data of candidates
	        - output_file_name : the name of the output file with the inferred black hole and stellar mass. If "-" is given, the output file name is given with the candidate names: If there is one candidate in the input file, output_file_name = "candidate-name.out". If there are more than two, output_file_name = "(first-candidate-name)_(last-candidate-name).out" where "first-candidate-name" ("last-candidate_name") is the name of the first (last) candidate listed in the input data file
		- c1 : value of c_1 (the characteristic distance scale of the emission region in units of the apocenter distance of the most tightly bound debris). Defualt = 1
		- Del_omega : Solid angle (in units of pi) of radiation from the emission region. Default = 2.0
		- N_sampling : The number of bin in mbh (black hole mass) and mstar (stellar mass) to find the range of mbh and mstar for given Lobs and Tobs. It also determines the grid resolution of the contour plots produced as output. If N_sampling is small, it will produce output more quickly. N_sampling > 100 is suggested.


	2) input data file (the name is given in "model_info.txt")

		The code can handle two formats (".txt" or ".csv"). The input data file should consists of 11 columns:


			Candidate_name   mstar_min   mstar_max   mbh_min   mbh_max    Lobs    Lobs_1    Lobs_2    Tobs    Tobs_1    Tobs_2


		- Candidate_name : you can put any name here. This name has nothing to do with mass inference, but will be shown in resulting figures.

		- mstar_min, mstar_max : These two defines the range of stellar mass [mstar_min, mstar_max], in units of solar mass, within which the solver tries to find the solution. If the solutions are not found, you can increase the range and re-run the code. 

		- mbh_min, mbh_max : the range of black hole mass, in units of solar mass [mbh_min, mbh_max], within which the solver tries to find the solution. 

		- Lobs  : the peak luminosity in units of erg/s. 
		- Lobs_1, Lobs_2 : These two is determined by the uncertainties in L. (IMPORTANT) L_1 is the uncertainty for the lower bound and L_2 is the uncertainty for the upper bound. In other words, Lobs - Lobs_1 and Lobs + Lobs_2 are the lower and upper bounds for Lobs, respectively. For example, if log10(Lobs)=43_{-0.1}^{+0.1}, Lobs = 10^43, Lobs_1 = 10^{43}-10^{43-0.1} = 2.1*10^42 and Lobs_2 = 10^{43+0.1} - 10^{43} = 2.6*10^42. 

		- Tobs : the black body temperature at the peak luminosity in units of K

		- Tobs_1, Tobs_2 : Same for Lobs_1 and Lobs_2. See the definitions for Lobs_1 and Lobs_2.


		The input.csv file in "example" directory has two examples. you can take copy and paste the file to put it in the same level

4. Output data format

The code produces three different output files. 2) a plot for the solutions on the mbh - mstar plane for each candidate given in the input data file. 3) A plot for the solutions and the error bars on the mbh - mstar plane for all candidates given in the input data file. All these output files are saved in "output" directory.

	1) table with the inferred black hole mass for the candidates

		This output file consists of 13 columns. The first seven columns are identical to those in the input data file. The rest columns are the inferred black hole mass[mbh], its uncertainties [mbh_1,mbh_2], the inferred stellar mass [mstar] and its uncertainties[mstar_1,mstar_2]. The uncertainties for the inferred masses are defined the same way as those for Lobs and Tobs. The name of the output file is given in the "model_info.txt"

	2) solution figure for each candidate

		The plots are created as many as the number of the candidates given in the input data file. There are two strips as long as the solutions are within the ranges of mbh and mstar set in the input data file. The blue strip demarcates the solutions for the given Lobs and the red strip for the given Tobs. The X-hatched green region indicates the solutions for the given Lobs and Tobs. 

	3) solution figure for all candidates

		The name of the figure is  "c1_(c1_value)_del_omega_(del Omega value)_inferred_mass.png"
		"c1_value" and "del Omega value" are the values of c_1 and Del_omega, respectively.
		The error bars are determined to enclose the extreme values of the solutions. 

5. Remark
