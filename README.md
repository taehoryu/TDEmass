# TDEmass


1. Description

	TDEmass is developed to infer the black hole mass (mbh) and stellar mass (mstar) by solving the two non-linear equations derived in Ryu+2020 () 
  (Equations 9 and 10) assuming the main source of the observed bolometric luminosity is the heat dissipated by shocks near apocenter. 
  The input is the spectral data (peak UV/optical luminosity, denoted by Lobs, and black body temperature at peak, denoted by Tobs). 
  For each candidate, the code tries to solve the equations using three different methods.
	
	1) It repeatedly computes Lobs and Tobs from logarithmically sampled mbh and mstar (within the range defined in the input file, bin size = N_sampling) 
  and finds the ranges of mbh and mstar yielding Lobs and Tobs within their uncertainties. The central values for the ranges for mbh and mstar are used 
  as the first guess for the second method (solver1_LT). If the following two methods fail, the central values are given as the solutions.

	2) solver1_LT() : With the first guess, it computes mbh and mstar using the equations. It uses the mean value of the first guess and the first computed values 
  to get the second values. It repeats this procedure until they converge simultaneously.

	3) solver2_LT() : 2-dimensional Newtonian-Raphson method.


2. Two free parameters c_1 and Del_omega
Our model includes two free parameters, c_1 and Del_omega. c_1*a_0 is the distance from the black hole at which a significant amount of energy is dissipated by shocks. 
Here, a_0 is the apocenter distance for the orbit of the most tightly bound debris. Del_omega is the solid angle. 
So Del_omega * (c_1 * a_0)^2 corresponds to the surface area of the emitting region. Note that c1a0(Del Omega/ 4pi)^(1/2) is equivalent to what is called the "black body radius".
For fixed Lobs and Tobs, larger Del_omega leads to smaller mbh and mstar, but the dependence is weak.  
However, mbh and mstar are sensitive to c_1. For fixed Lobs and Tobs, mbh \propto c_1^{a} with a = -(1.2-2) and mstar \propto c_1^{b} with b = 0.8-1.5. 
The default values are c_{1}=1 and Del_omega=2pi. With these choices, it is likely that the solutions are correct within factors of a few (Ryu+2020). 


3. Download and run

	To download, clone the git repository using HTTP access :

		https://github.com/taehoryu/TDE_mass_inference.git

	To run, you first enter the directory "TDE_mass_inference". In the directory, the code is run as

		Python3 main.py

	The code has two source dependencies : scipy, numpy and matplotlib

4. Units

	The default units are:

	- Luminosity(Lobs, Lobs_1, Lobs_2) = erg/s
	- Temperature (Tobs, Tobs_1, Tobs_2) = K
	- black hole mass (mbh) = 10^6msol (msol : solar mass)
	- stellar mass (mstar) = msol
 

5. Input data format

	Before running, you should prepare two input files which should be placed in "input_file" directory:  one that contains the value of key parameters and 
  the other that contains the spectral data for the events that you are interested in. 
 
	1) model_info.txt : define key parameters. 
	
		- inputdata_file_name : the name of the input data file containing the spectral data of candidates.
    
	        - output_file_name : the name of the output file with the inferred black hole and stellar mass. If "-" is given, the output file name is given 
                         with the candidate names: If there is one candidate in the input file, output_file_name = "candidate-name.out". If there are 
                         more than two, output_file_name = "(first-candidate-name)_(last-candidate-name).out" where "first-candidate-name" ("last-candidate_name") 
                         is the name of the first (last) candidate listed in the input data file.
                         
		- c1 : value of c_1 (the characteristic distance scale of the emission region in units of the apocenter distance of the most tightly bound debris). Defualt = 1.
    
		- Del_omega : Solid angle (in units of pi) of radiation from the emission region. Default = 2.0.

		- mbh_range : the range of black hole mass, in units of 10^{6} solar mass, within which the solver tries to find the solution. 

		- mstar_range : These two define the range of stellar mass, in units of solar mass, within which the solver
                             tries to find the solution.

        
	2) input data file (the name is given in "model_info.txt")

		The code can handle two formats (".txt" or ".csv"). The input data file should have 11 columns:

			Candidate_name   Lobs    dLobs-    dLobs+    Tobs    dTobs-    dTobs+


		- Candidate_name : you can put any name here. This name does not affect the mass inference, but will simply be shown in solution figures (see 4. output).

		- Lobs  : the peak luminosity in units of erg/s. 
    
		- dLobs-, dLobs+: These two values are determined by the uncertainties in Lobs. (IMPORTANT) Lobs - (dLobs-) and Lobs + (dLobs+) are the lower and upper limits 
                       for Lobs, respectively. For example, if log10(Lobs)=43_{-0.1}^{+0.1}, Lobs = 10^43, dLobs- = 10^{43}-10^{43-0.1} = 2.1*10^42 and 
                       dLobs+ = 10^{43+0.1} - 10^{43} = 2.6*10^42. This definition applies to the uncertainties for all quantities produced in the output file.

		- Tobs : the black body temperature at the peak luminosity in units of K.

		- dTobs-, dTobs+ : Same for dLobs- and dLobs+. See the definitions for dLobs- and dLobs+.

		The input.csv file in "example" directory has two examples. you can take copy and paste the file and put it in "input_file" directory.

6. Output

	The code produces three different output files. All these output files are saved in "output" directory.

	1) text file with the inferred black hole mass for the candidates

		This output file consists of 19 columns. The first seven columns are identical to those in the input data file. 
                The rest columns are the inferred black hole mass[mbh], its uncertainties (dmbh-,dmbh+), the inferred stellar mass [mstar], 
                its uncertainties (dmstar-,dmstar+), t0 and its uncertainties (dt0-, dt0+) and a0 and its uncertanties (da0-,da0+). 
                Here, t0 is the characteristic mass return time of the most tightly bound debris. 
                The uncertainties for the masses, t0 and a0 are defined the same way as those for Lobs (see Section 5-2).
                The name of the output file is given in the "model_info.txt"

	2) solution figure on the mbh - mstar plane for each candidate

		The plots are created as many as the number of the candidates given in the input data file. There are two strips as long as the solutions 
               are within the ranges of mbh and mstar set in the input data file. The blue strip demarcates the solutions for the given Lobs and the red strip 
               for the given Tobs. The X-hatched green region indicates the solutions for the given Lobs and Tobs. 

	3) solution figure on the mbh - mstar plane for all candidates

		The name of the figure is  "c1_(c1_value)_del_omega_(del Omega value)_inferred_mass.png"
		"c1_value" and "del Omega value" are the values of c_1 and Del_omega, respectively.
		The error bars are determined to enclose the extreme values of the solutions. 

7. Remark

- E-mail address for inquiry : udraeo@gmail.com

8. Please acknowledge Ryu+2020()  by citing ....
