#!/usr/env/python3
import sys
import seaborn as sns
import re
import math
import matplotlib.pyplot as plt
import shutil
import os
import numpy as np
import subprocess
import scipy.stats as st
from element_information import return_elements

# Modified statistical distribution that we need
def stat_dist(alpha,l):
    return (2 *l + 1)*np.exp(alpha*l)

# Routine to read the experimental data from custom 
# defined format as show in exp_intens.dat
# We also return the number of s transitions
def read_exp_data(element_name):

    # Open file
    f = open("exp_intens.dat", "r")
    lines = f.readlines()

    s_trans = 0
    intensity_dict = {}
    counter = 0
    for line in lines:

		# Split into each transition
        splitline = line.split()

        # Find the appropriate element and count how many starting
        # energy levels we have
        if(splitline[0] == element_name):
          counter += 1

		  # Start point is the energy level in which we start
          start_point = splitline[1]

		  # Put the entire line of transitions into the dictionary
          intensity_dict[start_point] = splitline[2:]
          if(str(splitline[2:][int(start_point)-2])) != "-2+":
              s_trans += 1
    return intensity_dict, s_trans

# Return an array containing a normal distribution for a given
# experimental mean and std
def generate_experimental_gaussian(mean, std, linspace):
    
    gaussian = []

    for i in linspace:
        gaussian.append(1/(math.sqrt(2*np.pi*std**2)) * np.exp((-(i - mean)**2)/(2*std**2)))

    return np.array(gaussian)

# Get all of the transitions and put them in a dictionary
def parse_akylas_transitions(out_name):
    akylas_intensities = {}
    # We want to loop over our entire set of transitions given by Akylas
    for i in range(2,21):

        # Make an empty dict
        akylas_intensities[i] = []

        # Now loop over the end points
        for j in range(1,i):

            # Here we need a space between the = and number
            if ((i <= 9) and (j <= 9)):
                line = subprocess.check_output(["grep", "N1= "+str(i)+", N2= "+str(j)+",", out_name])
                float_line = float(line.decode("utf-8").split()[7].split("(")[0])
            elif ((i > 9) and (j <= 9)):
                line = subprocess.check_output(["grep", "N1="+str(i)+", N2= "+str(j)+",", out_name])
                float_line = float(line.decode("utf-8").split()[6].split("(")[0])
            elif ((i > 9) and (j > 9)):
                line = subprocess.check_output(["grep", "N1="+str(i)+", N2="+str(j)+",", out_name])
                float_line = float(line.decode("utf-8").split()[5].split("(")[0])
            akylas_intensities[i].append(float_line)
    return akylas_intensities
    
if(len(sys.argv) == 1):
    print("Fatal error: Element must be provided at command line")
    sys.exit()

# For now we are keeping nmax at 20
nmax = 20

# Get the element from command line
element = sys.argv[1]
l_distribution = sys.argv[2]

# Get all the experimental data, but it's not split 
# into intensities and errors yet
exp_intens = {}
exp_intens, num_s_transitions = read_exp_data(element)

# First let's generate an array containing the mean population
distribution = np.zeros(nmax)

# Custom distribution used
# Unused behaviour for now
if(int(l_distribution) == -1):
    print ("Custom l-distribution chosen.")
    nop_line = "NOP   OR   "+str(l_distribution)
    experiment_name = element+"custom"+alpha

# Modified statistical distribution
elif(int(l_distribution) == 0):
    print ("Modified statistical distribution chosen")
    if(len(sys.argv) <= 3):
        print ("Fatal error: Exponent must be provided")
        sys.exit()
    else:
        alpha = sys.argv[3]
        print ("Exponent alpha = ", alpha)
        nop_line = "NOP   OR   "+str(l_distribution)
        experiment_name = element+"_stat_"+alpha
        alpha = float(alpha)
        for i in range(nmax):
            distribution[i] = stat_dist(alpha, i)

# Name of the folder which will contain the results
os.makedirs(experiment_name, exist_ok=True)

# This gets all our element specific parameters for Akylas
z_line, zs_line, a_line, be_line = (return_elements(element))

# Compile the code here
os.system("gfortran -g --std=legacy -fd-lines-as-comments -ffpe-trap=zero,invalid,denorm,underflow -fbacktrace -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -funsafe-math-optimizations -fno-align-commons -fmax-errors=1 akylas_cascade.f")

# More necessary lines for the input file
# which are required for every run
stop_line = "STO(P)NE"
xeq_line = "XEQ   NE"

# Move basic file to directory
shutil.copyfile("basic_input", experiment_name+"/basic_input")

# Move into the proper directory
os.chdir(experiment_name)

# Number of samples to use in each degree of freedom
# for the sensitivity analysis
# USER SPECIFY THIS NUMBER TO OBTAIN THE REPRODUCIBLE RESULT
iters = 100

# This dictionary will contain all possible n level transitions
# from Akylas
intensity_table = []

# This contains a dictionary for each sample of the parameters
list_of_dicts = []

# Counter for keeping track of the samples
counter = 0

# Loop over a set of bounds from the user specific exponent
if(int(l_distribution) == 0):

	# Set the range we are looking at
    sensitivity_range = 3.0
    upper_bound = alpha * ( 1 + sensitivity_range)
    lower_bound = alpha * ( 1 - sensitivity_range)

	# Zero the alpha array
    alpha_vals = np.zeros(iters)

	# Loop over our range of parameter space
    for i in np.arange(lower_bound, upper_bound, (upper_bound - lower_bound)/iters):

		# Round for naming the files
        j = np.round(i,3)
        name = "input_"+experiment_name+"_"+str(j)
        out_name = "output_"+experiment_name+"_"+str(j)
        nmx_line = "NMX   OR   "+str(nmax)+" "+str(i)

		# Copy the basic input and append all of the relevant
		# input lines to it
        shutil.copyfile("basic_input", name)
        f = open(name, "a")
        f.write(z_line+"\n")
        f.write(zs_line+"\n")
        f.write(a_line+"\n")
        f.write(be_line+"\n")
        f.write(nop_line+"\n")
        f.write(nmx_line+"\n")
        f.write(xeq_line+"\n")
        f.write(stop_line+"\n")
        f.close()
		
		# Run the code with our chosen distribution
        os.system("../a.out < "+name+" >"+out_name)
		
		# Collect all of the intensities
        intensity_table = parse_akylas_transitions(out_name)
		  
        # Keep track of the intensities and parameter value
        alpha_vals[counter] = j
        list_of_dicts.append(intensity_table)
        counter += 1
 
# Loss as a function of parameter space
loss_array = np.zeros(iters)

# Loop over each of our samples
for l in range(iters):
    loss_func = 0

    # Loop over each experimental start point
    for i,j in zip(exp_intens.values(), exp_intens.keys()):

        counter = 1

        # Loop over each transition 
        for k in i:

            # Split up into intensity and error
            intens = k.split("+")[0]
            err = k.split("+")[1]

            # If we reached a -2, then we can move on to the next line as
            # no experimental transition exists
            if(float(intens) < -1):
                continue

            # Add the contribution to the loss function, scaled with the 
            # experimental uncertainty
            loss_func += ((float(list_of_dicts[l][int(j)][counter-1]) - float(intens))/float(err))**2
            counter += 1
        loss_array[l] = loss_func

# Find the fitted parameter which corresponds to the minimum
# value of the loss function
fitted_alpha = alpha_vals[np.argmin(loss_array)]
print("Optimal parameter = ", fitted_alpha)

# Make a new directory to contain data we will
# want to plot with GNUplot
figure_data = "figure_data"
os.makedirs(figure_data, exist_ok=True)

# Name of the lossplot file
loss_data = "lossplot_"+experiment_name+".dat"

if(int(l_distribution) != 4):
    # Write the loss function to file
    f = open(loss_data, "w")
    if(int(l_distribution) == 0):
        for i in range(len(alpha_vals)):
            f.write(str(alpha_vals[i])+" "+ str(loss_array[i])+"\n")
    elif(int(l_distribution) == 2 or int(l_distribution) == 16):
        for i in range(len(a_vals)):
            # for j in range(len(b_vals)):
            f.write(str(a_vals[i])+" "+str(b_vals[i]) +" "+str(loss_array[i])+"\n")
    f.close()

    # Move the loss plot to a folder
    shutil.move(loss_data, figure_data+"/"+loss_data)
