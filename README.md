Supplemental Materials for "The Macroeconomic Impact of NAFTA Termination”
Joseph Steinberg, University of Toronto

This document describes the source data and computer programs used to create the results in Canadian Journal of Economics Manuscript 201842, “The Macroeconomic Impact of NAFTA Termination.”  These supplementary materials contain all of the files necessary to replicate the analysis in the paper. Section 1 describes the source data. Section 2 describes the Python scripts used to process and analyze the data, and the scripts that create the tables and figures that contain empirical results. Section 3 describes the C program used to calibrate and solve the model. Section 4 describes the Python scripts that create the tables and figures that contain model results.
1.	Data
The raw data described below are contained in the data folder of the supplementary materials.
1.1	World Input Output Database
The 2016 release of the World Input Output Database is comprised of 15 Stata data files, each of which contains a world input-output matrix for a given year. Each file lists gross output, intermediate input usage, and final demand for 43 countries and 56 industries. For more information on the construction of the WIOD, see http://www.wiod.org/release16. I use only the 2014 data in this project, so I include only the input-output matrix for that year in this supplement. The file is WIOT2014_October16_ROW.dta.
1.2	MFN tariffs from 2014
The file wto_mfn_tariffs_nafta_2014_hs12.csv, which was downloaded from the WTO website, contains each NAFTA country’s MFN average applied tariff at the 6-digit HS 2012 level.
1.3	Bilateral tariffs from 1993
The file trains_data_extract.csv, which was downloaded from the TRAINS online database, contains bilateral average applied tariffs for each of the 6 NAFTA import relationships (USA-MEX, USA-CAN, MEX-USA, MEX-CAN, CAN-USA, CAN-MEX) at the 6-digit HS 2012 level.
1.4	Bilateral trade flows from 2014 (HS 2012 version)
The file comtrade_nafta_2014_hs12.csv, which was downloaded from the COMTRADE website, contains bilateral goods trade flows for each of the 6 NAFTA import relationships for the year 2014 at the 6-digit HS 2012 level.


1.5	Bilateral trade flows from 2014 (HS 1996 version)
The file comtrade_nafta_2014_hs96.csv, which was downloaded from the COMTRADE website, contains bilateral goods trade flows for each of the 6 NAFTA import relationships for the year 2014 at the 6-digit HS 1996 level. This second version of the COMTRADE data is required to compute sector-level pre-NAFTA tariffs because the TRAINS data from 1993 uses the 1996 version of the HS classification system.
2. Python scripts (empirical)
All of the Python scripts described in this section are contained in the folder programs/python in the supplementary materials. They require the Python version 2 codebase, along with the NumPy, SciPy, and Matplotlib libraries. The output of these scripts is stored in the programs/python/output folder. Some of these scripts (namely 2.1) produce intermediate files that are used by other scripts (these files are also contained in the supplement). Thus, script 2.1 must be run first before running the other scripts in this section.
2.1 wiod_preproc.py
This script processes the raw WIOD data (item 1.1) aggregating industries according to the specification in Table 1 and aggregating all non-NAFTA countries into a single rest of the world. The aggregated data are saved to three pickle files, wiod_m.pik, wiod_f.pik, wiod_vg.pik, that contain intermediate inputs, final demand, and gross output/value added, respectively.
2.2 table_wiod_key_facts.py
This script loads the pickle files created by 2.1 and uses them to create a LaTeX file that contains Table 3 (key_facts.tex).
2.3 plots_data.py
This script loads the pickle files created by 2.1 and uses them to create Figures 1 through 3 (fig1_bilateral_trade.pdf, fig2_sectoral_trade.pdf, fig3_trade_balances.pdf).
2.4 iomat.py
This script loads the pickle files created by 2.1 and uses them to create an aggregated input-output matrix that is used as an input to the C program that solves the model (see section 3). This input file is called iomat.txt.
2.5 table_tariffs.py
This script loads the HS-6 level tariff and trade data (items 1.3–1.5), merges them, and creates LaTeX files that contain Table 2 (tariffs.tex, which contains sector-level MFN tariffs on trade between NAFTA countries) and Table 7 (tariffs_old.tex, which contains sector-level tariffs on trade between NAFTA countries based on bilateral applied tariffs from 1993). This script also creates text files that are used as inputs to the C program: tariffs.txt (the baseline tariffs); tariffs_alt1.txt (sensitivity analysis with increased U.S. protectionism); and tariffs_old.txt (pre-NAFTA tariffs).
2.6 elasticities.py
This script loads the raw, disaggregated WIOD data (item 1.1), and uses them to compute sectoral the trade elasticities shown in panel (b) of Table 4. This script simply prints output to the screen; the user must manually type the output into the LaTeX file for the table (assigned_params.tex).
3. C program
The general equilibrium model described in section 3 and calibrated in section 4 is solved using a computer program written in C. This program is contained in the folder programs/c. The source code is contained in the folder programs/c/src. The binary executable, which is created by compiling the program, is contained in the folder programs/c/bin. The output of the program, which is created by running the executable, is contained in the folder programs/c/output.
In addition to the standard C codebase, this program requires the GNU gcc compiler (https://gcc.gnu.org/), the GNU GSL library (https://www.gnu.org/software/gsl/), OpenMP (https://www.openmp.org/), and the Intel Math Kernel Library (https://software.intel.com/en-us/mkl). I used the Ubuntu 16.04 Linux distribution to compile and run the program, and I cannot guarantee that the program will work without modifications in Windows or other operating systems.
The program uses OpenMP parallelization to evaluate the Jacobian matrix of the equilibrium system, and a computer with a large number of cores is required to run the program in a reasonable amount of time. The program also requires a large amount of memory to store this matrix and other program objects; the memory requirement depends on the number of threads used in the program. I used a 40-core Xeon workstation and the program required about 12 GB of RAM.
3.1 Compiling the program
The supplementary materials contain a makefile that compiles the source code into an executable. To compile the program, the user should run make from the command line in the programs/c folder. The executable created by the makefile is located in the folder programs/c/bin.
3.2 Running the program
Once the program has been compiled, the baseline quantitative exercise can be run by typing the command ./bin/nafta from the command line in the programs/c/ folder. The program has a number of command-line options that direct it to perform various sensitivity analyses instead of the baseline exercise. For example, to run the version of the analysis without capital adjustment costs, use the command ./bin/nafta -k. To see a full listing of the command line options, run ./bin/nafta -h (the “h” stands for help). To run the baseline analysis and all of the sensitivity analyses, you can use the bash script run (type ./run from the programs/c folder). Note that this will take a very long time!

3.3 Source code
The folder programs/c/src contains all of the header and source files required to compile the program (aside from the external libraries described above).
3.3.1 globals.h
This header file, which is included by all other files in the program, contains a number of macros, short inline functions, and global variables that are used throughout the program.
3.3.2 main.c
This source file contains the programs’s main() function, along with a driver for the quantitative exercise and a routine that processes command-line arguments.
3.3.3 calibrate.h and calibrate.c
This header and source file are responsible for calibrating the model’s parameters. The calibrate() function drives the entire calibration process, which loads the text files created by scripts 2.4 and 2.5, assigns values to the parameters listed in Table 4, and then performs the steps described in sections 4.2 and 4.3.
3.3.4 eqm.h and eqm.c
This header and source file are responsible for solving for the model’s equilibrium. The header file defines an eqm struct that contains all of the equilibrium variables, and the function solve_eqm() finds the values of these variables that satisfy the model’s equilibrium conditions, which are contained in the function eval_eqm_conds().
3.3.5 exporters.h and exporters.c
This header and source file contain the structures and functions that manage the heterogeneous-exporters portion of the model. The header file defines an exporter_vars struct that contains the value function, export participation rate, and threshold productivities, and the functions in the source file perform the steps needed to update these objects as the macroeconomy evolves (represented by equations 18 through 21 in the paper’s text).
3.3.6 solver.h, solver.c, gnewton.h, and gnewton.c
These header and source files contain a parallelized implementation of the gnewton method for solving square nonlinear systems from the GNU GSL library. The key difference between my version of this method and the one contained in the GSL is that mine uses OpenMP parallelization to evaluate the Jacobian matrix, and then calls a parallel linear solver from the Intel MKL to solve for the Newton step implied by this matrix.
3.4 Program output
The program creates several output files, all written to the programs/c/output folder. Note that the supplementary materials contain all of these output files, so there is no need for the user to re-run the C program if he or she wants to create different tables or figures using the output.

3.4.1 Model parameters
The file params.txt contains all of the calibrated values of the matrix parameters that are not reported in the text of the paper (e.g. the Armington shares, the fixed exporting costs, etc.)
3.4.2 Output of baseline quantitative analysis
The files named vars<x>_<y>.csv, where <x> = {0,1} and <y> = {usa,can,mex,row}, contain the output of the baseline quantitative exercise. The vars0_<y> files contain the equilibrium values for each country in the model in the benchmark equilibrium where NAFTA remains in force forever. The vars1_<y> files contain the termination equilibrium values.
3.4.3 Output of sensitivity analyses
The output of the various sensitivity analyses are contained in files with similar naming conventions as the files described in 3.4.2. These files, however, contain additional suffixes that denote the sensitivity analysis. The files are vars<x>_<y>_<z>.csv, where
•	<x> = {0,1}
•	<y> = {usa,can,mex}
•	<z> = {no_k_adj_cost, no_l_adj_cost, no_trd_adj_cost, eqkappa, nokappa, fix_tb, noio, sym_te, cobb_douglas, tariff_alt4, tariff_alt4_v2, tariff_alt5, tariff_alt2, tariff_alt1, tariff_alt7, cp_combo}
The <z> suffixes are listed in the order that the sensitivity analyses appear in Table 6.
4. Python scripts (model results)
All of the Python scripts described in this section are also contained in the folder programs/python in the supplementary materials. They require the same codebase as the scripts described in section 2. These scripts can be run in any order, but they musts be run after the C program because they rely on its output.
4.1 table_model_lr_changes.py
This script loads the output of the C program (namely, the files vars<x>_<y>.csv, where <x> = {0,1} and <y> = {usa,can,mex}), and produces Table 5, which shows the long-run effects of NAFTA termination in the baseline exercise (model_lr_changes.tex).
4.2 plots_model.py
This script loads the baseline C program output (the same files listed in 3.4.2) and produces Figures 4 through 6, which plot the baseline model’s long-run changes and transition dynamics (fig4_LR_trade.pdf, fig5_LR_realloc.pdf, fig6_trd_macro_dyn.pdf).
4.3 plots_sens.py
This script loads the baseline C program output (the files listed in 3.4.2) and produces Figure 7, which plots the transition dynamics in the model without trade adjustment costs (fig7_trd_macro_dyn_no_trd_adj_cost.pdf).
4.4 table_welfare.py
This script loads the C program output for the baseline model (files in 3.4.2) and all of the sensitivity analyses (files in 3.4.3), and creates Table 6 (welfare.tex), which contains all of the welfare calculati
