#! /usr/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os.path
from scipy import interpolate
import sys

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'serif'

# versions to compare
versions = []
popsize = [-1] #todo 
#if len(sys.argv) <= 1:
#	
#	print("Provide at least one argument, versions to plot (you provided " + str(len(sys.argv)-1) + ").")
#	exit()
	
#else:

#	iterargs = iter(sys.argv)
#	next(iterargs)
#	for arg in iterargs:
#		versions.append((arg))
	

# Benchmark Settings
runs = range(30)
max_fevals = [100000,100000,100000,100000,1000000,1000000,1000000, 1000000] 
plot_limit_fevals_min = [ 1000, 1000, 1000, 1000, 10000, 10000,10000, 10000]
objectives = [ "SYMPART1","SYMPART3", "SSUF1", "SSUF3"]
objective_names = [ "SYM-PART1 (n=2)","SYM-PART3 (n=2)", "SSUF1 (n=2)", "SSUF3 (n=2)", "MinDist (m=2,n=10)", "MinDist (m=2,n=30)","MinDist (m=3,n=5)", "MinDist (m=3,n=10)"]
versions = ["0","100","101","110", "111"]

max_igd = [ 1.42e-3, 1.42e-3, 3.13e-4, 3.68e-4,7e-4, 7e-4, 5.6e-3, 5.6e-3]
max_igdx = [ 3.95e-3, 3.25e-3, 5.9e-3,7.58e-4, 9.8e-4, 9.8e-4, 6.85e-3, 6.85e-3]



fig, ax = plt.subplots(3,len(objectives),figsize=[4*len(objectives),9],sharex='col',sharey='row')
colors = ["r","b","c","y","g","k","m"]

for o in range(len(objectives)):

	plot_limit_log_fevals_max = (max_fevals[o])
	plot_limit_log_time_max = 5;

	interp1pts_fevals = np.arange(2.0, plot_limit_log_fevals_max, 100)
	interp1pts_time = np.arange(1.0, plot_limit_log_time_max, 1)

	print("Reading objective " + objectives[o])
	
	## IGD subplot

	if(o == 0):
		ax[0, o].set_ylabel('$\leftarrow$ IGD')
	
	ax[0, o].set_title(objective_names[o],fontweight='bold')
	ax[0, o].grid()
	ax[0, o].set_ylim([1e-4, 1e1])
	ax[0, o].set_xlim([plot_limit_fevals_min[o], plot_limit_log_fevals_max])
	ax[0, o].axhline(y=max_igd[o], color='k', linestyle=':',label='limit')
	
	## IGDX subplot

	if(o == 0):
		ax[1,o].set_ylabel('$\leftarrow$ IGDX')
	
	# plt.title(objectives[o])
	ax[1, o].grid()
	ax[1, o].set_ylim([1e-4, 1e2])
	#ax[1, o].set_xlim([1000, plot_limit_log_fevals_max])
	
	ax[1, o].axhline(y=max_igdx[o], color='k', linestyle=':',label='limit')
	
	## MR subplot

	if(o == 0):
		ax[2][o].set_ylabel('MR $\\rightarrow$')
	
	# plt.title(objectives[o])
	ax[2, o].grid()
	ax[2, o].set_ylim([0 ,1.1])
	#ax[2, o].set_xlim([1000, plot_limit_log_fevals_max])
	ax[2, o].set_xlabel('Function Evaluations')
	ax[2, o].axhline(y=1.0, color='k', linestyle=':',label='limit')
		
	for vi in range(len(versions)):
		
		f_fevals_igd = np.empty([len(runs),len(interp1pts_fevals)]);
		f_fevals_igpd = np.empty([len(runs),len(interp1pts_fevals)]);
		f_time_igd = np.empty([len(runs),len(interp1pts_time)]);
		f_time_igpd = np.empty([len(runs),len(interp1pts_time)]);
		f_fevals_pareto_sets = np.empty([len(runs),len(interp1pts_fevals)]);
		collected_data_lines = 0
		
		for ri in range(len(runs)):
			
			# read a single line
			random_seed = 100 + runs[ri]
			
			if(o == 5 or o == 7):
			
				try:
					file = open("runlogs/runlogs_version" + str(versions[vi]) + "/statistics_" + objectives[o] + "_" + str(random_seed) + ".txt","r")
				except:
					print("runlogs/runlogs_version" + str(versions[vi]) + "/statistics_" + objectives[o] + "_" + str(random_seed) + ".txt -- Not Found.")
					break
					
			else:
			
				try:
					file = open("../longrun/runlogs/runlogs_version" + str(versions[vi]) + "/statistics_" + objectives[o] + "_" + str(random_seed) + ".txt","r")
				except:
					print("../longrun/runlogs/runlogs_version" + str(versions[vi]) + "/statistics_" + objectives[o] + "_" + str(random_seed) + ".txt -- Not Found.")
					break
				
			data = file.readlines()
			
			# note, first line is title, last line is empty
			lines = len(data)

			# for each line, 
			time = np.empty(lines)
			fevals = np.empty(lines)
			igd = np.empty(lines)
			igpd = np.empty(lines)
			pareto_sets = np.empty(lines)
			
			for i in range(lines-1):

				# +1 to skip first line
				words = data[i+1].split()

				time[i+1] = (float(words[2]))
				fevals[i+1] = (float(words[1]))
				
				if(o == 6 or o == 7):
					igd[i+1] = (float(words[7])) # is the default front found?
					igpd[i+1] = (float(words[8])) # is the default front found?
					pareto_sets[i+1] = float(words[9])
				else:
					igd[i+1] = (float(words[6])) # is the default front found?
					igpd[i+1] = (float(words[7])) # is the default front found?
					pareto_sets[i+1] = float(words[8])

			# set the initial dpfs the same as the first one achieved.
			time[0] = 0
			fevals[0] = 2.0
			igd[0] = igd[1]
			igpd[0] = igpd[1]
			pareto_sets[0] = pareto_sets[1]
			
			
			if fevals[lines-1] < plot_limit_log_fevals_max:
				# break;
				fevals[lines-1] = plot_limit_log_fevals_max
			
			if fevals[lines-1] >= plot_limit_log_fevals_max:
				collected_data_lines += 1
			
				# here we have a vector of (fevals,igd) in the log10-domain
				# now interpolate it. 

				# these are interpolation objects
				fevals_igd = interpolate.interp1d(fevals,igd)
				fevals_igpd = interpolate.interp1d(fevals,igpd)
				fevals_pareto_sets = interpolate.interp1d(fevals,pareto_sets)
				time_igd = interpolate.interp1d(time,igd)
				time_igpd = interpolate.interp1d(time,igpd)

				f_fevals_igd[ri,:] = np.array(fevals_igd(interp1pts_fevals))
				f_fevals_igpd[ri,:] = np.array(fevals_igpd(interp1pts_fevals))
				f_fevals_pareto_sets[ri,:] = np.array(fevals_pareto_sets(interp1pts_fevals))
				# f_time_igd[ri,:] = np.array(time_igd(interp1pts_time))
				# f_time_igpd[ri,:] = np.array(time_igpd(interp1pts_time))
			
				# line = plt.plot(fevals,dpfs, '--', linewidth=1,color=colors[l],alpha=0.2)

		sym = '-'
		label = ''
		# collected all data, now make plot

		if(float(versions[vi]) == 0):
			sym = '--'
			label = "MAMaLGaM"

		if(float(versions[vi]) == 100):
			label = "MO-HillVallEA-MAM"
			sym = '-'

		if(float(versions[vi]) == 101):
			label = "MO-HillVallEA-MAMu"
			sym = '-'

		if(float(versions[vi]) == 110):
			label = "MO-HillVallEA-iMAM"
			sym = '-'

		if(float(versions[vi]) == 111):
			label = "MO-HillVallEA-iMAMu"
			sym = '-'
			
		## IGD PLOT
		ax[0, o].loglog(interp1pts_fevals,f_fevals_igd[0:collected_data_lines,:].mean(axis=0), sym, linewidth=1,color=colors[vi],label=label)

		if len(f_fevals_igd[0:collected_data_lines,:]) > 0 :
			ax[0, o].fill_between(interp1pts_fevals,f_fevals_igd[0:collected_data_lines,:].min(axis=0),f_fevals_igd[0:collected_data_lines,:].max(axis=0), facecolor=colors[vi], alpha=0.1, rasterized=True)

			
		## IGDX PLOT
		ax[1, o].loglog(interp1pts_fevals,f_fevals_igpd[0:collected_data_lines,:].mean(axis=0), sym, linewidth=1,color=colors[vi],label=label)

		if len(f_fevals_igpd[0:collected_data_lines,:]) > 0 :
			ax[1, o].fill_between(interp1pts_fevals,f_fevals_igpd[0:collected_data_lines,:].min(axis=0),f_fevals_igpd[0:collected_data_lines,:].max(axis=0), facecolor=colors[vi], alpha=0.1, rasterized=True)

			
		## MR PLOT
		ax[2, o].semilogx(interp1pts_fevals,f_fevals_pareto_sets[0:collected_data_lines,:].mean(axis=0), sym, linewidth=1,color=colors[vi],label=label)
	
		if len(f_fevals_pareto_sets[0:collected_data_lines,:]) > 0 :
			ax[2, o].fill_between(interp1pts_fevals,f_fevals_pareto_sets[0:collected_data_lines,:].min(axis=0),f_fevals_pareto_sets[0:collected_data_lines,:].max(axis=0), facecolor=colors[vi], alpha=0.1, rasterized=True)

				
		#if(o == 0):
		#	ax[2, o].legend(loc='upper left',ncol=6, bbox_to_anchor=(0.4,-0.25))
		
			
	plt.draw()

plt.savefig("longrun_part1.png",bbox_inches='tight')
plt.savefig("longrun_part1.pdf",bbox_inches='tight')
plt.close()
