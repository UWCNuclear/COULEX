# COULEX CODES: A collection of Coulomb exication codes to simplifly matters regarding low energy Coulomb excitation


G2CEtest.cpp: Code to calculate correlated error using GOSIA2
Code to determine correlated error using GOSIA2 using method outlined in GOSIA 2012 manual ( T. Czosnyka, D. Cline, C. Y. Wu, A. B. Hayes, \textit{Gosia user manual for simulation and analysis of Coulomb excitation experiments}, 109-114 (2012). )

Depedencies: GOSIA2 and c++11 or higher needs to be installed. 
Gosia2 input files: For projectile excitation the GOSIA2 minimazation input files is required. 
Too compile code: g++ -Wall G2CEtest.cpp -o CET.
Usage: ./CET A.MINI.inp 
Output: A.MINI.error.dat
