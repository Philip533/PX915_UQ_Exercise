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
alpha = sys.argv[2]
iters = int(sys.argv[3])

nop_line = "NOP   OR   0"
experiment_name = element+"_"+alpha

# Get all the experimental data, but it's not split 
# into intensities and errors yet
exp_intens = {}
exp_intens, num_s_transitions = read_exp_data(element)

# Name of the folder which will contain the results
os.makedirs(experiment_name, exist_ok=True)

# This gets all our element specific parameters for Akylas
z_line, zs_line, a_line, be_line = (return_elements(element))

# Compile the code here
os.system("gfortran -g --std=legacy -fd-lines-as-comments -ffpe-trap=zero,invalid,denorm,underflow -fbacktrace -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -funsafe-math-optimizations -fno-align-commons -fmax-errors=1 akylas_cascade.f")

# More necessary lines for the input file
# which are required for every run
# XEQ starts the cascade, and stop stops the code
stop_line = "STO(P)NE"
xeq_line = "XEQ   NE"

# Move basic file to directory
shutil.copyfile("basic_input", experiment_name+"/basic_input")

# Move into the proper directory
os.chdir(experiment_name)

# This dictionary will contain all possible n level transitions
# from Akylas
intensity_table = []

# This contains a dictionary for each sample of the parameters
list_of_dicts = []

# Counter for keeping track of the samples
counter = 0

# Loop over a set of bounds from the user specific exponent

# Set the range we are looking at
sensitivity_range = 3.0
alpha = float(alpha)
upper_bound = alpha * ( 1 + sensitivity_range)
lower_bound = alpha * ( 1 - sensitivity_range)

# Zero the alpha array
alpha_vals = np.zeros(iters)

# Holds 2 values to calculate the gradient
grad_values = []
grad_step = 1e-6 * 0.16
opt_alpha = 0.16
alpha_std = 8e-4

# Loop over our range of parameter space
# for i in np.arange(lower_bound, upper_bound, (upper_bound - lower_bound)/iters):
for i in (opt_alpha, opt_alpha + grad_step):

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
    # alpha_vals[counter] = j
    # list_of_dicts.append(intensity_table)
    grad_values.append(intensity_table[2][0])
    counter += 1

first_order_std = (np.sqrt(((grad_values[1] - grad_values[0])/grad_step)**2 * alpha_std **2))
mean_val = grad_values[0]
 
# Loss as a function of parameter space
loss_array = np.zeros(iters)

# Loop over each of our samples
# for l in range(iters):
#     loss_func = 0

#     # Loop over each experimental start point
#     for i,j in zip(exp_intens.values(), exp_intens.keys()):

#         counter = 1

#         # Loop over each transition 
#         for k in i:

#             # Split up into intensity and error
#             intens = k.split("+")[0]
#             err = k.split("+")[1]

#             # If we reached a -2, then we can move on to the next line as
#             # no experimental transition exists
#             if(float(intens) < -1):
#                 continue

#             # Add the contribution to the loss function, scaled with the 
#             # experimental uncertainty
#             loss_func += ((float(list_of_dicts[l][int(j)][counter-1]) - float(intens))/float(err))**2
#             counter += 1
#         loss_array[l] = loss_func

# # Find the fitted parameter which corresponds to the minimum
# # value of the loss function
# fitted_alpha = alpha_vals[np.argmin(loss_array)]
fitted_alpha = 0.16
print("Optimal parameter = ", fitted_alpha)

# Now we must use the fitted parameter and a normal distribution to run the calculation
# multiple times
num_samples = 10

# List which holds all the intensities for each sample i.e num_samples amount of dictionaries
normal_dist_intensities = []

# Loop over number of samples we want from a normal distribution centred around
# the fitted alpha value
for j in range(num_samples):
    sampled_alpha = np.random.normal(fitted_alpha, alpha_std)

    # Each iteration gets a file
    name = "input_"+experiment_name+"_"+str(j)
    out_name = "output_"+experiment_name+"_"+str(j)
    nmx_line = "NMX   OR   "+str(nmax)+" "+str(sampled_alpha)

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
    normal_dist_intensities.append(intensity_table)
    
lyman_alpha_intensities = np.zeros(num_samples)
for j in range(len(normal_dist_intensities)):
    lyman_alpha_intensities[j] = (normal_dist_intensities[j][2][0])
# print(lyman_alpha_intensities)
sns.kdeplot(lyman_alpha_intensities)
print("Intensity of 2->1 using brute force sampling =  ", np.mean(lyman_alpha_intensities) ," +- ", np.std(lyman_alpha_intensities))
print("Intensity of 2->1 using first order variance = ", mean_val ," +- ", first_order_std)
x = np.linspace(mean_val - 3*first_order_std, mean_val + 3*first_order_std, 100)
plt.plot(x, st.norm.pdf(x, mean_val, first_order_std))
plt.show()
