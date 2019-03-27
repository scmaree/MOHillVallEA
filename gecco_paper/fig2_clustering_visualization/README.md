------------ Step 1: Compile MO-HillVallEA --------------------

cd ./MOHillVallEA
make

------------ Step 2: Start experiments ------------------------
This reproduces Table 1 Go to the experiments folder,

cd ../gecco_paper/fig2_clustering_visualization

To generate the results for N = 250 and N = 10000 run respectively

python startupscript_single_run_analysis.py 250 <# of threads>
python startupscript_single_run_analysis.py 10000 <# of threads>


------------- Step 3: Generate Figures ------------------------

python plot_init_clusterings.py 250
python plot_init_clusterings.py 10000 

------------- Bonus step: Re-use raw GECCO data ---------------
The runlogs used in the GECCO paper are stored in runlogs_gecco. Rename folder to runlogs to generate figures based on this data.
