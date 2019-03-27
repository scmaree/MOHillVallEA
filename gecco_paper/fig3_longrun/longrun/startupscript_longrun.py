import subprocess
import time
import os
import shutil
import sys

# Runs the given version of MMMO on the entire benchmark
# the benchmark consists of 10 functions

if len(sys.argv) != 2:
	
	print("Provide 1 argument CPU Cores to use (you provided " + str(len(sys.argv)-1) + ").")
	exit()
	
else:
	maximum_number_of_processes = int(sys.argv[1])
	
	print("Running version on " + str(maximum_number_of_processes) + " cores.")

# Benchmark Settings
number_of_runs = 30
max_fevals_list = [100000, 100000, 100000, 100000, 1000000, 10000000]
problems_list = [12,14, 15, 16, 9, 19]
problem_variables = [2, 2, 2, 2, 10, 5]
elitist_archive_size_target = [1000, 1000, 1000, 1000, 1000, 2500]

# standard parameter settings
lower_user_range = str(-20)
upper_user_range = str(20)
maximum_number_of_populations = str(8)
popsize = str(-1)

vtr = str(1)
maximum_number_of_seconds = str(0)
number_of_subgenerations_per_population_factor = str(8)


# check the directory -> empty it if it is already full

versions = [0, 100, 101, 110, 111]
lopts = [0, 0, 1, 10, 11]

for i in range(len(versions)):


	version = versions[i]
	write_directory = "./runlogs/runlogs_version" + str(version) + "/"
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



# copy executable to runlog directory
src_exe = "../../../HillVallEA/mo-hillvallea.app"


# parallel running

processes = set()

for run in range(number_of_runs):

	print("Run "+ str(run))
	random_seed_str = str(100+run)
	
	for i in range(len(versions)):

		version = versions[i]
		lopt = lopts[i]
		write_directory = "./runlogs/runlogs_version" + str(version) + "/"

		for pi in range(len(problems_list)): #

			problem_index = str(problems_list[pi])
			number_of_parameters = str(problem_variables[pi])
			maximum_number_of_evaluations = str(max_fevals_list[pi])
			
			processes.add(subprocess.Popen([src_exe,"-s","-V",str(version), problem_index, number_of_parameters, str(lopt), lower_user_range, upper_user_range, maximum_number_of_populations, popsize, str(elitist_archive_size_target[pi]), str(elitist_archive_size_target[pi]), maximum_number_of_evaluations, vtr, maximum_number_of_seconds, number_of_subgenerations_per_population_factor, random_seed_str, write_directory], shell=False))
			
			while len(processes) >= maximum_number_of_processes:
				time.sleep(1)
				processes.difference_update([
					p for p in processes if p.poll() is not None])
