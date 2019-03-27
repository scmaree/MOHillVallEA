# MO-HillVallEA (GECCO 2019)
S.C. Maree, T. Alderliesten, P.A.N. Bosman


The Multi-Objective Hill-Valley Evolutionary Algorithm (MO-HillVallEA) is a real-valued multi-objective evolutionary algorithm specifically aimed for multi-modal optimization. It has been described initially in 

>**Real-valued Evolutionary Multi-modal Multi-objective Optimization by Hill-valley Clustering**
> S.C. Maree, T. Alderliesten, D. Thierens, and P.A.N. Bosman. 
> GECCO-2019, ACM Press, New York, New York, 2019.

---

# Getting started
Start by making a clone of the repository,

``` git clone https://github.com/SCMaree/HillVallEA```


 Call `make` to build using your favorite compiler in the directory `./MOHillVallEA`. This builds MO-HillVallEA with a command line interface and a number of benchmark functions pre-configured. Call `./mo-hillvallea.app` for a description. 
 
 For example, the following runs MO-HillVallEA-MAM on Mindist(2,2):
 
 `./mo-hillvallea.app -s -v -V 100 9 2 0 -20 20 1 -1 1000 100 30000 0 0 2 135 "./"`

This outputs the runlog, and writes it to `./statistics_genMEDmm_135.txt`
 
Furthermore, three python scripts have been provided to replicate the results and figures presented in the above mentioned paper. 


