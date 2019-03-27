#! /usr/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os.path
from scipy import interpolate
import sys

# versions to compare
versions = []

if len(sys.argv) <= 1:
	
	print("Provide at least one argument, versions to plot (you provided " + str(len(sys.argv)-1) + ").")
	exit()
	
else:

	iterargs = iter(sys.argv)
	next(iterargs)
	for arg in iterargs:
		versions.append((arg))
	

# Benchmark Settings
runs = range(31)
max_fevals = 30000
objectives = [ "TwoOnOne","OmniTest", "SYMPART1", "SYMPART2", "SYMPART3","SSUF1", "SSUF3"] #,"genMEDconcave","genMEDmm"] #TwoOnOne #BD2

for v in range(len(versions)):

	print("")
	print("Version " + versions[v])
	print("IGD    Std    IGDX.  Std    Gens  Time    Sets  Runs Problem")

	for o in range(len(objectives)):

		generations = []
		time = []
		fevals = []
		igd = []
		igpd = []
		sets = []
		
		number_of_runs = 0;

		for ri in range(len(runs)):
			
			# read a single line
			random_seed = 100 + runs[ri]
			
			try:
				file = open("runlogs/runlogs_version" + str(versions[v]) + "/statistics_" + objectives[o] + "_" + str(random_seed) + ".txt","r")
			except:
				#print("File  -- runlogs_version" + str(versions[v]) + "/statistics_" + objectives[o] + "_" + str(random_seed) + ".txt -- Not Found.")
				break
				
			data = file.readlines()
			
			# note, first line is title, last line is empty
			lines = len(data)

			if lines > 0:
				# +1 to skip first line
				words = data[lines-1].split()

				generations.append(float(words[0]))
				time.append(float(words[2]))
				fevals.append(float(words[1]))
				igd.append(float(words[6]))
				igpd.append(float(words[7]))
				sets.append(float(words[9]))

				number_of_runs += 1

		if(len(time) > 0):
			print("%1.3f"% (sum(igd)/len(igd)) 
				+ " " + "%1.3f"% (np.std(igd))
				+ " " + "%1.3f"% (sum(igpd)/len(igpd))
				+ " " + "%1.3f"% (np.std(igpd))
				+ " " + "%4.f"% (sum(generations)/len(generations))  
				+ "  " + "%4.1f"% (sum(time)/len(time))
				+ "  " + "%6.2f"% (sum(sets)/len(time)) 
				+ "  " + "%2.f"% len(time)
				+ "  " + objectives[o]
				)

			

