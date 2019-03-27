import subprocess
import time
import os
import shutil
import sys

# Runs a single run of the given versions 
# plots objective space of each generation
# versions

if len(sys.argv) != 3:
	
	print("Provide at 2 arguments: (1) init_popsize, (1) maximum_number_of_processes, (you provided " + str(len(sys.argv)-1) + " arguments).")
	exit()
	
else:

	init_popsize = float(sys.argv[1])
	version = sys.argv[1];
	lopt = "0"
	random_seed = "100"
	maximum_number_of_processes = int(sys.argv[2])
		
		

# Benchmark Settings
max_fevals = 10  * init_popsize
# problems_list = [10, 11, 12, 13, 14, 15, 16]
#problem_variables = [2,5,2,2,2,2,2]
problems_list = [10, 14, 15, 16, 19]
problem_variables = [2,2,2,2,2]

# standard parameter settings
lower_user_range = str(-20)
upper_user_range = str(20)
maximum_number_of_populations = str(1)
popsize = str(init_popsize)
elitist_archive_size_target = str(1000)
approximation_size_target = str(100)
maximum_number_of_evaluations = str(max_fevals)
vtr = str(1e-8)
maximum_number_of_seconds = str(0)
number_of_subgenerations_per_population_factor = str(2)


# check the directory -> empty it if it is already full


# parallel running
processes = set()
versions = [version]

for vi in range(len(versions)):
	
	# version = versions[vi]
	
	write_directory = "./runlogs/runlogs_single_run_version"+ str(version)+"/"
	
	if not os.path.exists(write_directory):
		os.makedirs(write_directory)
	else:
		for the_file in os.listdir(write_directory):
			file_path = os.path.join(write_directory, the_file)
			try:
				if os.path.isfile(file_path):
					os.unlink(file_path)
				#elif os.path.isdir(file_path): shutil.rmtree(file_path)
			except Exception as e:
				print(e)


	# copy executable to runlog directory
	src_exe = "../../MOHillVallEA/mo-hillvallea.app"
	#dst_exe = write_directory + "mmmo"
	#shutil.copy2(src_exe, dst_exe)


	print("Version "+ str(version))
	
	# run = 3
	# random_seed_str = str(100+run)
	
	for pi in range(len(problems_list)): #
	
		problem_index = str(problems_list[pi])
		number_of_parameters = str(problem_variables[pi])
		
		processes.add(subprocess.Popen([src_exe,"-s","-w","-v","-V",str(version), problem_index,number_of_parameters,lopt, lower_user_range, upper_user_range, maximum_number_of_populations, popsize, elitist_archive_size_target, approximation_size_target, maximum_number_of_evaluations, vtr, maximum_number_of_seconds, number_of_subgenerations_per_population_factor, str(random_seed), write_directory], shell=False))
		
		while len(processes) >= maximum_number_of_processes:
			time.sleep(10)
			processes.difference_update([
				p for p in processes if p.poll() is not None])
