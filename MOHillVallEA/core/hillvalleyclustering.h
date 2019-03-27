#pragma once

/*
 
 MO-HillVallEA
 
 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "hillvallea_internal.h"
#include "fitness.h"
#include "elitist_archive.h"
namespace hillvallea
{

  class hvc_t
  {

  public:

    // (de)initialize
    //----------------------------------------------------------------------------------
    hvc_t(fitness_pt fitness_function);
    ~hvc_t();

    void cluster(const population_t & pop, std::vector<population_pt> & subpopulations,  unsigned int  & number_of_evaluations, double & average_edge_length, bool add_test_solutions, bool recheck_elites, int optimizer_number, rng_pt & rng);
    void cluster_mo_distances(population_t & pop, std::vector<population_pt> & subpopulations, unsigned int  & number_of_evaluations, double & average_edge_length, bool add_test_solutions, bool recheck_elites, int optimizer_number, rng_pt & rng);
    
    bool hill_valley_test(const solution_t & sol1, const solution_t & sol2, const unsigned int max_trials, unsigned int & number_of_evaluations, std::vector<solution_pt> & test_solutions);
    bool hill_valley_test_single_objective(const solution_t & sol1, const solution_t & sol2, size_t objective_number, unsigned int max_trials, unsigned int & number_of_evaluations, std::vector<solution_pt> & test_solutions);
    bool hill_valley_test_point_single_objective(const hillvallea::solution_t &sol1, const hillvallea::solution_t &sol2, const hillvallea::solution_t &test, size_t objective_number);
    bool hill_valley_test_point_single_objective(double worst_objective_value, const hillvallea::solution_t &test, size_t objective_number);
    bool hill_valley_test_single_objective_reuse(const solution_t & sol1, const solution_t & sol2, size_t objective_number, unsigned int max_trials, unsigned int & number_of_evaluations, std::vector<solution_pt> & test_solutions);
    
    
  private:

    // input parameters
    //---------------------------------------------------------------------------------
    fitness_pt fitness_function;


  };


}
