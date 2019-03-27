import subprocess
import time
import os
import shutil
import sys

# Runs the given version of MMMO on the entire benchmark
# the benchmark consists of 10 functions

if len(sys.argv) != 4:
	
	print("Provide three arguments, version number, local optimizer, and CPU Cores to use (you provided " + str(len(sys.argv)-1) + ").")
	exit()
	
else:
	version = sys.argv[1]
	lopt = sys.argv[2]
	maximum_number_of_processes = int(sys.argv[3])
	
	print("Running version " + str(version) + " on " + str(maximum_number_of_processes) + " cores.")

# Benchmark Settings
number_of_runs = 31
max_fevals = 30000
problems_list = [10,12,13,14,15,16]
problem_variables = [2,2,2,2,2,2]
popsize = str(-1);
# standard parameter settings
lower_user_range = str(-20)
upper_user_range = str(20)
maximum_number_of_populations = str(1)
elitist_archive_size_target = str(1000)
appr_size_target = str(100)
maximum_number_of_evaluations = str(max_fevals)
vtr = str(1e-8)
maximum_number_of_seconds = str(0)
number_of_subgenerations_per_population_factor = str(2)
write_directory = "runlogs/runlogs_version"+ str(version)+"/"


# check the directory -> empty it if it is already full

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

# parallel running
processes = set()

for run in range(number_of_runs):

	print("Run "+ str(run))

	random_seed_str = str(100+run)
	
	for pi in range(len(problems_list)): #

		problem_index = str(problems_list[pi])
		number_of_parameters = str(problem_variables[pi])
		
		processes.add(subprocess.Popen([src_exe,"-V",str(version), problem_index, number_of_parameters,str(lopt),lower_user_range, upper_user_range, maximum_number_of_populations, popsize,  elitist_archive_size_target, appr_size_target, maximum_number_of_evaluations, vtr, maximum_number_of_seconds, number_of_subgenerations_per_population_factor, random_seed_str, write_directory], shell=False))
		
		while len(processes) >= maximum_number_of_processes:
			time.sleep(0.1)
			processes.difference_update([
				p for p in processes if p.poll() is not None])
