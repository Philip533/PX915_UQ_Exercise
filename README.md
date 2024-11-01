This repository contains the code and data used for calculations in the corresponding report. Full details of the theory can be found there; this is a simple document for guiding the user through a basic calculation.

Dependencies:

Seaborn-0.13.2
TinyGP-0.3.0
SciPy 1.13.1
NumPy 2.0.2
Gfortran 11.4.1
Make sure that X-forwarding is enabled for your SSH connection if you want the final figure to show up automatically.

Available elements:
P, Cl, Mg, Al, In, Au, Ho

For this set of calculations, we are using a modified statistical distribution with a free parameter alpha to model the X-ray intensity of a 2→1 transition in a given element. We want to use the Akylas cascade code provided to both brute force, and perform first order variance estimation in order to compare the statistics.

First, a least squares regression must be performed using experimental data to find a suitable mean value of alpha. This can be done by doing:

        `python3 run_experiment.py <element symbol> <sensible guess for alpha> -1`

Sensible initial guesses of values of alpha can be taken from either Hartmann [1] or Vogel [2] for the corresponding element. (Hint: Tables I-VI of Vogel and Tables 3-7 of Hartmann)

This prints out an optimised value of alpha; make a note of this value, as it used for the next step.

Next, the code needs to be run again, but with different input parameters:

        `python3 run_experiment.py <element symbol> <fitted alpha from previous run> <number of samples>`

This will run the Akylas code <number of samples> times, using the fitted parameter from before. The intensities for the 2→1 transition are then calculated, and a histogram is plotted.

A plot containing both the first order estimation and the brute force samples will be created, along with the mean and standard deviations of both distributions. 
This file will be called "histogram.pdf" and will be located in the root of the repository.

[1] F. J. Hartmann, R. Bergmann, H. Daniel, H.-J. Pfeiffer, T. von Egidy, and W. Wilhelm, Measurement of the muonic x-ray cascade in Mg, AI, In, Ho, and Au, Z Physik A 305, 189 (1982).

[2] P. Vogel, Muonic cascade: General discussion and application to the third-row elements, Phys. Rev. A 22, 1600 (1980).

