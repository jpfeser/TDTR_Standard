TDTR_Standard
=============

Tools for analyzing time-domain thermoreflectance experiments (standard refers to Concentric pump-probe and radially anisotropic stack of materials and interfaces).  Unidirectional Heat Flow

I.  In the "Main_Calls" folder, there are two files "TDTR_Main.m" and "TDTR_Fit.m" that can accomplish nearly everything by calling the other subroutines 
*THESE ARE THE PROGRAMS USERS NORMALLY EDIT/RUN FOR THEIR SPECIFIC NEEDS.  I DON'T RECOMMEND ALTERING ANY OF THE OTHER PROGRAMS.*

The "TDTR_Main.m" can by default:
	1)  Fit TDTR data (more on the details later...)
	2)  Calculate errorbars for real or simulated data fits
	3)  Make sensitivity plots
	4)  Calculate the Steady-state temperature rise in an experiment

II.  Here is the workflow for using the "TDTR_Main.m" program:
	Step 1: Type all of the properties of the system (materials properties, layer thicknesses, spot sizes, laser rep rate, etc. in the first section of the program labeled "TYPE THERMAL SYSTEM PARAMTERS HERE." The first element of each vector (lambda, h, C) is the top layer...i.e. where the laser hits.  The last element is the substrate, and the program always assumes the substrate is semi-infinite.
	Step 2: Choose the range of time delays you care about.  The smaller the absolute value of the time delay, the longer the program will take to run.
	Step 3: Choose your program options:  this is where you choose whether you will load/fit data, or calculate sensitivity plots.

To fit data to a thermal model (i.e. load a .txt file and fit data), set "importdata=1;" Otherwise type "importdata=0;" to skip loading real data.
To calculate sensitivity plots, based on the parameters you typed in during step 1, set "senseplot=1;"  To skip this step, set "senseplot=0;"  I recommend always looking at the sensivity plots.
	Step 4: Choose whether you would like to calculate the errorbars for the fitting.  If you want to calculate errorbars, type "ebar=1;" Otherwise, set "ebar=0;" to skip it.  Note:  if you load real data for fitting, the errorbar will be calculated using that data, but if you don't load real data, the program will simulate data using the parameters from step 1 and can still calculate the estimated errorbar.

	In order to calculate the errorbars, the program needs to know what the "percent error" is for all of the parameters you typed into step 1.  Type them in in the sections labeled "ERRORBAR OPTIONS."  For example, if the heat capacity of layer 3 has a 5% errorbar, then type "C(3)=0.05"  You can skip the calculation for layers which you know won't introduce large errors in the solution using the "_consider" variables.  For example, to not consider layer 3's thermal conductivity, type "lambda_consider(3)=0" otherwise, you should consider every layer except the layers that you are solving for (because changing your initial guess for that layer should not affect the solution).


	Step 5 (If you are fitting TDTR or calculating errorbars):  You need to alter the "TDTR_Main.m" program and the "TDTR_Fit.m" program so that the program knows which variables you are fitting for (and it can do simultanous/multivariable fits!).  Here's how you do it:

