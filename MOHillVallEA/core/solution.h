#pragma once

/*
 
 MO-HillVallEA
 
 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "hillvallea_internal.h"
#include "param.h"

namespace hillvallea
{

  class solution_t {

  public:

    // essential data members
    //-----------------------------------------
    vec_t param;             // position of the solution, i.e. the coordinate vector
    vec_t obj;               // fitness values
    std::vector<vec_t> gradients; // one gradient-vector per objective
    double penalty;          // penalty value ( set to > 0 if the solution is infeasible )
    size_t rank;             // some ranking of solutions
    elitist_archive_pt elite_origin; // register from which elitist achrive the sol came
    int cluster_number;
    int population_number;
    
    // constructor & destructor
    //-----------------------------------------
    solution_t();
    solution_t(size_t number_of_parameters);
    solution_t(size_t number_of_parameters, size_t number_of_objectives);
    solution_t(vec_t & param);
    solution_t(vec_t & param, vec_t & obj);
    solution_t(const solution_t & other);
    ~solution_t();

    // compare two solutions to see which is best
    //-----------------------------------------
    static bool better_solution_via_pointers(const solution_pt & sol1, const solution_pt & sol2);
    static bool better_solution_per_objective_via_pointers(const solution_pt & sol1, const solution_pt & sol2, size_t objective_number);
    
    static bool better_solution(const solution_t & sol1, const solution_t & sol2);
    static bool better_solution_per_objective(const solution_t & sol1, const solution_t & sol2, size_t objective_number);
    static bool better_solution_unconstraint(const solution_t & sol1, const solution_t & sol2);
    
    bool better_than(const solution_t & sol) const;
    bool better_than_unconstraint(const vec_t & other_obj) const;
    
    bool strictly_better_than(const solution_t & other) const;
    bool strictly_better_than_unconstraint(const vec_t & other_obj) const;

    bool better_than_unconstraint_per_objective(size_t objective_number, double other_obj) const;
    bool strictly_better_than_unconstraint_per_objective(size_t objective_number, double other_obj) const;
    bool better_than_per_objective(const solution_t & other, size_t objective_number) const;
    
    static bool better_rank_via_pointers(const solution_pt & sol1, const solution_pt & sol2);
    
    bool same_objectives(const solution_t & other) const;

    double minimal_objective_distance(const std::vector<solution_pt> & sols, const vec_t & obj_ranges) const;
    double minimal_param_distance(const std::vector<solution_pt> & sols) const;
    double transformed_objective_distance(const solution_t & other, const vec_t & obj_ranges) const;
    double transformed_objective_distance(const vec_t & other_obj, const vec_t & obj_ranges) const;

    // problem dimensions
    //-----------------------------------------
    size_t number_of_parameters() const;
    size_t number_of_objectives() const;
    void set_number_of_parameters(size_t number);
    void set_number_of_objectives(size_t number);
    
    // Euclidean distance to another solution in parameter space
    //-----------------------------------------
    double param_distance(solution_t & sol2) const;
    double param_distance(vec_t & param2) const;

  };



}
