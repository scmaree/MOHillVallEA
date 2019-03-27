------------ Step 1: Compile MO-HillVallEA --------------------

cd ./MOHillVallEA
make

------------ Step 2: Start experiments ------------------------
This reproduces Figure 3, the long run. It is a bit of a lazy implementation in two parts
because my script writes the folder based on the problem name, and the MinDist problem is
repeated multiple times with different settings. 

Go to the experiments folder,

cd ../gecco_paper/fig3_longrun/longrun
python startupscript_longrun.py <# of threads to be used>

and for the second part

cd ../longrun2
python startupscript_longrun.py <# of threads to be used>

------------- Step 3: Generate Figures ------------------------

python plot_longrun_results_part1.py
python plot_longrun_results_part2.py

------------- Bonus step: Re-use raw GECCO data ---------------
The runlogs used in the GECCO paper are stored in runlogs_gecco. Rename folder to runlogs to generate figures based on this data.
