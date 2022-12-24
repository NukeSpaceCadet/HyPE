# HyPE
A hybrid rocket simulation code (Hybrid Performance Estimator, HyPE)

This code uses lookup tables for propellant combos (varying O/F ratios and chamber pressures). These can be generated
using Propep or CEA. To add new propellants, a new table will have to be generated. A table for one propellant cobo is included
and can be used as a template.

The input to the code is a plain-text file. A template is included for a whip-cream nitrous engine, and can be run as-is to test the code.
The code generates a plain-text output file with temperatures, pressures, mass flow rates, and specific impulse vs time, and saves plots
as a png bitmap.

to run the code: ./python HyPE.py input_filename.in
