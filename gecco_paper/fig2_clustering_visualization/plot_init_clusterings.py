#! /usr/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os.path
import seaborn
import sys

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'serif'

# Runs a single run of the given versions 
# plots objective space of each generation

if len(sys.argv) <= 1:
	
	print("Provide at least 1 argument: (1) init_popsize (you provided " + str(len(sys.argv)-1) + ").")
	exit()
	
else:

	version = sys.argv[1]
	random_seed = "100"


objectives = ["SYMPART1","SYMPART2","SYMPART3","OmniTest","TwoOnOne","SSUF1","SSUF3","genMEDmm","triangles2"] #[ "ZDT1", "ZDT2","ZDT3","ZDT4","ZDT6","BD1","BD2s","genMEDconvex","genMEDconcave","genMEDmm"]
objective_names = ["SYM-PART1 (m=2)", "SYM-PART2 (m=2)", "SYM-PART3 (m=2)", "Omni-Test (m=2)",'Two-On-One (m=2)', "SSUF1 (m=2)", "SSUF3 (m=2)", "MinDist (m=2)","MinDist (m=3)"]
xlimmin = [-20, -20, -20, 0, -3, 1, 0,-4,-5]
xlimmax = [20, 20, 20, 6, 3, 3, 1,4,5]
ylimmin = [-20, -20, -20, 0, -3, -1, 0,-4,-5]
ylimmax = [20, 20, 20, 6, 3, 1, 2,4,5]

obj1min = [0,0,0,-5,0.6,0,0,0,0]
obj1max = [4,4,4,0,2,1,1,2,3]
obj2min = [0,0,0,-5,0,0,0,0,0]
obj2max = [4,4,4,0,0.6,1,1,2,3]



#objectives = ["OmniTest"]
#xlimmin = [-4]
#xlimmax = [4]
#ylimmin = [-4]
#ylimmax = [4]
#obj1min = [0]
#obj1max = [2]
#obj2min = [0]
#obj2max = [2]

objectives_order = [4, 2, 5,6,8] #[4,3,0,1,2,5,6]


write_directory = "initial_plots/"

if not os.path.exists(write_directory):
	os.makedirs(write_directory)
#else:
#	for the_file in os.listdir(write_directory):
#		file_path = os.path.join(write_directory, the_file)
#		try:
#			if os.path.isfile(file_path):
#				os.unlink(file_path)
#			#elif os.path.isdir(file_path): shutil.rmtree(file_path)
#		except Exception as e:
#			print(e)


fig, ax = plt.subplots()
fig.set_tight_layout(True)
fig.set_figheight(3)
fig.set_figwidth(15)

plot_nr = 1; 


# for each objective create a new figure:
for o in objectives_order:

	objective = objectives[o]
	print("Problem " + objective)


	colors = ["r","b","g","c","m","y","aqua","aquamarine","fuchsia","gold","coral","plum","chocolate","indigo","maroon","silver","tan","wheat","sienna","olive","pink","greenyellow"]


	if(o == 8):
		colors = ["r","b","g","c","m","y","aqua","aquamarine","fuchsia","gold","coral","plum","chocolate","indigo","maroon","silver","tan","wheat","sienna","olive","pink","greenyellow"]

	if(o == 2):
		colors = ["r","b","g","c","m","y","aqua","aquamarine","fuchsia","gold","coral","plum","chocolate","indigo","maroon","silver","tan","wheat","sienna","olive"]
	
	# Parameter space subplot
	plt.subplot(1,len(objectives_order),plot_nr)
	plot_nr += 1
	plt.title(objective_names[o])
	plt.xlabel('$x_0$')
	plt.ylabel('$x_1$')
	plt.xlim([xlimmin[o], xlimmax[o]])
	plt.ylim([ylimmin[o], ylimmax[o]])
	ax.grid(linestyle='-', linewidth=1)
				


	# read default front
	filename = "./runlogs/runlogs_single_run_version"+str(version)+"/pareto_set_" + objective + "_" + str(random_seed)  + ".txt"
			
	pareto_set_available = False

	try: 
		file = open(filename,"r")
		pareto_set_available = True

	except:
		print("No default front. ")
		print("File " + filename + " -- Not Found.")
	

	if(pareto_set_available):
		
		data = file.readlines()

		lines = len(data)

		pareto_set = np.empty([lines,2]);
		pareto_front = np.empty([lines,2]);
		
		for i in range(lines):
		
			words = data[i].split()
			columns = len(words)
			number_of_parameters = columns - 4

			pareto_set[i,0] = float(words[0])
			pareto_set[i,1] = float(words[1])

			pareto_front[i,0] = float(words[number_of_parameters])
			pareto_front[i,1] = float(words[number_of_parameters+1])


	gen = 0;

	# read the population
	population_available = False
	filename = "./runlogs/runlogs_single_run_version"+version+"/population_00000_generation_" +  str(gen).zfill(5) + "_" + objective + "_" + str(random_seed)  + ".txt"
	try:
		file = open(filename,"r")
		population_available = True;

	except:
		print("File " + filename + " -- Not Found.")
		continue


	if(population_available):

		data = file.readlines()
		
		lines = len(data)
		
		population_obj = np.empty([lines,2])
		population_param = np.empty([lines,2])
		population_cluster = []

		for i in range(lines):
		
			words = data[i].split()
			columns = len(words)
			
			number_of_parameters = columns - 4

			population_param[i,0] = float(words[0])
			population_param[i,1] = float(words[1])

			population_obj[i,0] = float(words[number_of_parameters])
			population_obj[i,1] = float(words[number_of_parameters+1])
			
			population_cluster.append(colors[int(float(words[number_of_parameters + 3])%len(colors))])


	# plot stuff


	if(pareto_set_available):
		plt.plot(pareto_set[:,0],pareto_set[:,1],'k.',markersize=3,rasterized=True)

	if(population_available):
		plt.scatter(population_param[:,0],population_param[:,1], facecolors=population_cluster, marker='.',rasterized=True)



plt.legend()
plt.draw()
plt.savefig(write_directory + "version" + str(version) + "_" + str(random_seed) + ".png",bbox_inches='tight')
plt.savefig(write_directory + "version" + str(version) + "_" + str(random_seed) + ".pdf",bbox_inches='tight')
plt.close()






