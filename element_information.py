#!/usr/bin/python3
import sys
import re
import numpy as np

# This routine looks up the element of choice
# and returns the lines needed for the input file
# Z, ZS, A and BE, which are the element specific parameters
# BE are the binding energies for the Z-1 element
def return_elements(ele):

    zero_binding = "9.9E19"
    # Grab the binding energy file
    with open('element_information/binding.dat') as f:
        lines = f.read().splitlines()

    # List to contain the information for element
    element = []

    # Flag to check if an element has been found or not
    done = False

    # Get the appropriate element
    for i in range(len(lines)):
        split = []
        split = (lines[i].split())

        # Ignore comments
        if (split[0] == "#"):
            continue
        
        # If element matches input, then grab the previous line
        if (split[2] == ele):
            element = lines[i-1].split()
            done = True
            break

    # Catch incorrect name
    if (done == False):
        print("ERROR IN ELEMENT NAME", ele, "PROGRAM STOPPING")
        sys.exit()

    # All elements have a K binding energy
    k_energy = element[3]

    # Convert to floats
    element_float =  [float(i) for i in element[4:]]

    # Here we have the full set of L and M energies
    if(len(element[4:]) >= 8):
        l_energy = (sum(element_float[:3])/3)
        m_energy = (sum(element_float[3:8])/5)

    # Here we have full L but non full M
    elif(3 < len(element[4:]) < 8):
        l_energy = (sum(element_float[:3])/3)
        m_energy = (sum(element_float[3:])/len(element_float[3:]))

    # Here we catch the case of full L but empty M
    elif(len(element[4:]) == 3):
        l_energy = (sum(element_float[:3])/3)
        m_energy = zero_binding

    # Here we have empty M and non full L
    elif(0 < len(element[4:]) < 3):
        l_energy = (sum(element_float)/len(element_float))
        m_energy = zero_binding

    # No L energies at all
    elif(len(element[4:]) == 0):
        l_energy = zero_binding
        m_energy = zero_binding

    # Atomic charge
    charge = int(element[1])+1

    # Now get atomic mass table
    with open('element_information/atomic_masses.dat') as f:
        lines = f.read().splitlines()

    # Get the appropriate atomic mass
    for i in range(len(lines)):
        split_line = lines[i].split()
        if (ele == split_line[1]):
            mass = (lines[i].split()[3])

            # Remove the brackets that occur in some of the elements
            if(re.search("\(", str(mass))):
                mass = mass.split("(")[0]
            break

    # Charge of the element
    z_line = "Z     NE   "+str(charge)

    # Deal with negative charges
    charge_k = charge - 1 if (charge - 1) > 0 else 0
    charge_l = charge - 3 if (charge - 3) > 0 else 0
    charge_m = charge - 9 if (charge - 9) > 0 else 0

    # Effective charges for K, L, M shell
    zs_line = ("ZS    NE   "+str(charge_k) + " " + str(charge_l)+" "+ str(charge_m))

    # Atomic mass
    a_line  = ("A     OR   "+str(mass))

    # Binding energies for Z-1
    be_line = ("BE    NE  "+" "+ str(k_energy)+" "+str(l_energy)+" "+ str(m_energy))

    return z_line, zs_line, a_line, be_line

