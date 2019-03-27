------------ Step 1: Compile MO-HillVallEA --------------------

cd ./MOHillVallEA
make

------------ Step 2: Start experiments ------------------------
This reproduces Table 1

Go to the experiments folder,

cd ../gecco_paper/table1_igdx

To generate the results for MAMaLGaM, run

python startupscript_tablescript.py 0 0 <# of threads to be used>

To generate the results for MO-HillVallEA, run

python startupscript_tablescript.py 100 0 <# of threads to be used>


------------- Step 3: Generate Figures ------------------------

python print_tables.py 0 100

------------- Bonus step: Re-use raw GECCO data ---------------
The runlogs used in the GECCO paper are stored in runlogs_gecco. Rename folder to runlogs to generate figures based on this data.
