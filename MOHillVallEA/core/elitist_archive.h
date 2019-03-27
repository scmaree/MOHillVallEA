#pragma once

/*
 
 MO-HillVallEA
 
 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "hillvallea_internal.h"
#include "solution.h"
#include "population.h"

namespace hillvallea
{


  /* 
    Elitist archive
  
    The elitist archive keeps a copy of all non-dominated solutions.
    It is a population with some extra functionality.

  */

  class elitist_archive_t : public population_t
  {

  public:

    // Rule of three
    //-----------------------------------------------------------------------------------------------
    elitist_archive_t( size_t target_size, rng_pt rng); // default objective discretization
    elitist_archive_t(const elitist_archive_t & other);
    ~elitist_archive_t();

    // data members
    //-----------------------------------------------------------------------------------------------
    size_t target_size;
    bool use_niching;
    
    // Update the archive
    //-----------------------------------------------------------------------------------------------
    int updateArchive(const solution_pt sol, bool make_copy_of_sol);
    int updateArchive(const solution_pt sol);
    
    // discretization of the objective space to reduce the number of solutions in the archive
    //-----------------------------------------------------------------------------------------------
    void adaptArchiveSize();

    // restart
    //-----------------------------------------------------------------------------------------------
    void clear();

    // check if a solution dominates
    //-----------------------------------------------------------------------------------------------
    bool solutionHasImproved(const hillvallea::solution_t & sol) const;
    double get_best_objective_values_in_elitist_archive(size_t obj_index);

    // Multi-Modal Archive
    //-----------------------------------------------------------------------------------------------
    elitist_archive_pt initNewArchive() const;
    void addArchive(elitist_archive_pt elitist_archive);
    void removeEmptyClusters();
    
    size_t size()  const;
    size_t number_of_clusters() const;
    size_t actualSize() const;
    void copyElitesToClusters(std::vector<cluster_pt> & clusters, size_t max_number_of_elites, const vec_t & objective_ranges, rng_pt rng) const;
    size_t addElitesToPopulation(population_t & population, int max_number_of_elites);
    
    void computeApproximationSet(size_t approximation_set_size, std::vector<optimizer_pt> & optimizers, elitist_archive_pt & elitist_archive, bool use_parameter_space_diversity, bool terminate_pops);
    
    void set_use_parameter_distances(bool value);
    void set_use_greedy_selection(bool value);
    size_t collectSubarchives(std::vector<population_pt> & subpopulations);

    void initArchiveForEachSubpop(std::vector<population_pt> & subpopulations) const;
    
    void objectiveRanges(vec_t & objective_ranges); // overload from population_t
    
  private:

    // private data members
    //-----------------------------------------------------------------------------------------------
    bool use_parameter_distances;
    bool use_greedy_selection;
    bool objective_discretization_in_effect;
    bool parameter_discretization_in_effect;
    vec_t objective_discretization;
    vec_t parameter_discretization;
    vec_t best_objective_values_in_elitist_archive;
    rng_pt rng;
    
    std::vector<elitist_archive_pt> clustered_archive;

    // add a solution to the archive without checking (adds solutions in the first 'nullptr' spot)
    //-----------------------------------------------------------------------------------------------
    void getAllSols(std::vector<solution_pt> & all_sols);
    void getAllSols(std::vector<solution_pt> & all_sols, std::vector<elitist_archive_pt> & origin);

    void addToArchive(const solution_pt sol, size_t & insert_index, bool make_copy_of_sol);
    void removeFromArchive(size_t sol_index);

    // To reduce the number of solutions in the archive
    //-----------------------------------------------------------------------------------------------
    bool sameObjectiveBox(const solution_t & sol1, const solution_t & sol2) const;

    
    void adaptObjectiveDiscretization();
    void adaptSizeBygreedyScatteredSubsetSelection();
    void adaptObjectiveDiscretization_mm();
    void adaptSizeBygreedyScatteredSubsetSelection_mm();
    
    void removeSolutionNullptrs();

    


  };

}
