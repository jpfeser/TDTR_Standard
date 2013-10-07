TDTR_Standard
=============

Tools for analyzing time-domain thermoreflectance experiments (standard refers to Concentric pump-probe and radially anisotropic stack of materials and interfaces).  Unidirectional Heat Flow

I.  In the "Main_Calls" folder, there are two files "TDTR_Main.m" and "TDTR_Fit.m" that can accomplish nearly everything by calling the other subroutines 
*THESE ARE THE PROGRAMS USERS NORMALLY EDIT/RUN FOR THEIR SPECIFIC NEEDS.  I DON'T RECOMMEND ALTERING ANY OF THE OTHER PROGRAMS.*

The "TDTR_Main.m" can by default:
	1)  Fit TDTR data (more on the details later...requires help of the program "TDTR_Fit.m")
	2)  Calculate errorbars for real or simulated data fits
	3)  Make sensitivity plots
	4)  Calculate the Steady-state temperature rise in an experiment
	5)  Save the results of all of the above, plus a plot of the model vs data fit.

II.  Workflow for using the "TDTR_Main.m" program:

	Step 1: Type all of the properties of the system (materials properties, layer thicknesses, spot sizes, laser rep rate, etc. in the first section of the program labeled "TYPE THERMAL SYSTEM PARAMTERS HERE." The first element of each vector (lambda, h, C) is the top layer...i.e. where the laser hits.  The last element is the substrate, and the program always assumes the substrate (last layer) is semi-infinite.

	Step 2: Choose the range of time delays you care about.  The smaller the absolute value of the time delay, the longer the program will take to run.

	Step 3: Choose your "program options":  this is where you choose whether you will load/fit data, or calculate sensitivity plots.

To fit data to a thermal model (i.e. load a .txt file and fit data), set "importdata=1;" Otherwise type "importdata=0;" to skip loading real data.
To calculate sensitivity plots, based on the parameters you typed in during step 1, set "senseplot=1;"  To skip this step, set "senseplot=0;"  I recommend always looking at the sensivity plots.
	Step 4: Choose whether you would like to calculate the errorbars for the fitting.  If you want to calculate errorbars, type "ebar=1;" Otherwise, set "ebar=0;" to skip it.  Note:  if you are loading/fitting data, the errorbar will be calculated using that data, but if you aren't fitting experimental data, the program will simulate data using the parameters from step 1 and can still calculate the estimated errorbar.  The way the program calculates errorbars is by redoing the fitting using different "fixed" parameters. For example, so suppose you aren't fitting the heat capacity because you "know" it already...if there is an "errorbar" associated with the heat capacity, then that will cause an error in the fitting.  The program asks for the "errorbar" associated with each "known" parameter and re-simulates the fitting using the worst case variation in the "known" parameters.  The difference between the "re-fit" and the "original fit" is the errorbar for the fit...add up the sum of the squares for all errorbars and you get the total errorbar.  

	In order to calculate the errorbars, the program needs to know what the "percent error" is for all of the parameters you typed into step 1.  Type them in in the sections labeled "ERRORBAR OPTIONS."  For example, if the heat capacity of layer 3 has a 5% errorbar, then type "C(3)=0.05"  You can skip the calculation for layers which you know won't introduce large errors in the solution using the "_consider" variables.  For example, to not consider layer 3's thermal conductivity, type "lambda_consider(3)=0" otherwise, you should consider every layer except the layers that you are solving for (because changing your initial guess for that layer should not affect the solution).


	Step 5: (Optional:  Only used if you are fitting TDTR data or calculating errorbars):  You need to alter the "TDTR_Main.m" program and the "TDTR_Fit.m" program so that the program knows which variables you are fitting for (it can do simultanous/multivariable fits!).  Here's how you do it (two steps):

A) In the "TDTR_Fit.m" program (the purpose of this program is to evaluate the sum-of-the-squares of the residuals, comparing experiment and simulation), you need to define which variables to "fit."  The vector of fitting variables is called "X" in th program.  So for example, to fit for the thermal conductivity of the 3rd layer and the thermal conductivity of the 2nd layer, you would type:

lambda(3) = X(1); %this might represent substrate thermal conductivity
lambda(2) = X(2); %this might represent thermal interface conductance

B) In the "PROGRAM OPTIONS" section of the "TDTR_Main.m" program, you will need to alter a variable called  "X_Fit_Variables" that is used to set the initial guesses for the "fitting parameters." The easiest way to do this, is to type something like (continuing the example from above):

X_Fit_Variables = [lambda(3) lambda(2)];

%OR

X_Fit_Variables(1) = lambda(3); %this represents initial guess for substrate thermal conductivity
X_Fit_Variables(2) = labmda(2); %this represents initial guess for thermal interface conductance


