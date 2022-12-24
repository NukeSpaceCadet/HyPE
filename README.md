# HyPE
A hybrid rocket simulation code (Hybrid Performance Estimator, HyPE)

This code uses looku[ tables for propellant combos (varying O/F ratios and chamber pressures). 
To add new propellants, a new table will have to be generated. A table for one propellant cobo is included
and can be used as a template.

The input to the code is a plain-text xml file. A template is included, and should run as-is to test the code.

to run the code: $python HyPE.py input_filename.in
