

/*
 
 MO-HillVallEA
 
 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "solution.h"

namespace hillvallea
{

  // initialize solution
  //----------------------------------------------
  solution_t::solution_t()
  {
    this->param.resize(0);
    this->obj.resize(0);
    this->penalty = 0.0;
    this->rank = -1;
    this->elite_origin = nullptr;
    this->cluster_number = -1;
    this->population_number  = -1;
  }

  solution_t::solution_t(size_t number_of_parameters)
  {
    this->param.resize(number_of_parameters, 0.0);
    this->obj.resize(0);
    this->penalty = 0.0;
    this->rank = -1;
    this->elite_origin = nullptr;
    this->cluster_number = -1;
    this->population_number  = -1;
  }

  solution_t::solution_t(size_t number_of_parameters, size_t number_of_objectives )
  {
    this->param.resize(number_of_parameters, 0.0);
    this->obj.resize(number_of_objectives, 1e308);
    this->penalty = 0.0;
    this->rank = -1;
    this->elite_origin = nullptr;
    this->cluster_number = -1;
    this->population_number  = -1;
  }

  solution_t::solution_t(vec_t & param)
  {
    this->param = param;
    this->obj.resize(0);
    this->penalty = 0.0;
    this->rank = -1;
    this->elite_origin = nullptr;
    this->cluster_number = -1;
    this->population_number  = -1;
  }

  solution_t::solution_t(vec_t & param, vec_t & obj)
  {
    this->param = param;
    this->obj = obj;
    this->penalty = 0.0;
    this->rank = -1;
    this->elite_origin = nullptr;
    this->cluster_number = -1;
    this->population_number  = -1;
  }

  solution_t::solution_t(const solution_t & other)
  {
    this->param = other.param;
    this->obj = other.obj;
    this->penalty = other.penalty;
    this->rank = other.rank;
    this->elite_origin = other.elite_origin;
    this->cluster_number = other.cluster_number;
    this->population_number  = other.population_number;
  }

  // delete solution
  // nothing to be done
  //----------------------------------------------
  solution_t::~solution_t() {}

  // problem dimensions
  //-----------------------------------------
  size_t solution_t::number_of_parameters() const {
    return param.size();
  }
  size_t solution_t::number_of_objectives() const {
    return obj.size();
  }

  void solution_t::set_number_of_parameters(size_t number) {
    this->param.resize(number);
  }
  void solution_t::set_number_of_objectives(size_t number) {
    this->obj.resize(number);
  }


  // comparison for solution_t pointers
  // returns true if sol1 is strict better than sol2
  //-----------------------------------------------
  // defined as static!

  bool solution_t::better_rank_via_pointers(const solution_pt & sol1, const solution_pt & sol2)
  {
    return (sol1->rank < sol2->rank);
  }

  bool solution_t::better_solution_via_pointers(const solution_pt & sol1, const solution_pt & sol2)
  {
    return better_solution(*sol1, *sol2);
  }

  // defined as static!
  // returns true of the first solution is strict better than the second
  bool solution_t::better_solution(const solution_t & sol1, const solution_t & sol2)
  {
    return sol1.better_than(sol2);
  }

  bool solution_t::better_than(const solution_t & other) const
  {
    
    if (this->penalty > 0) // this is infeasible
    {
      if (other.penalty > 0) // Both are infeasible
      {
        if (this->penalty < other.penalty)
          return true;
      }
    }
    else // this is feasible
    {
      if (other.penalty > 0) { // this is feasible and other is not
        return true;
      }
      else { // Both are feasible */
        return this->better_than_unconstraint(other.obj);
      }
    }
    
    return false;
  }
  
  
  bool solution_t::better_solution_per_objective_via_pointers(const solution_pt &sol1, const solution_pt &sol2, size_t objective_number)
  {
    return better_solution_per_objective(*sol1, *sol2, objective_number);
  }
  
  bool solution_t::better_solution_per_objective(const solution_t & sol1, const solution_t & sol2, size_t objective_number)
  {
    return sol1.better_than_per_objective(sol2, objective_number);
  }
  
  // defined as static!
  // returns true of the first solution is strict better than the second
  bool solution_t::better_solution_unconstraint(const solution_t & sol1, const solution_t & sol2)
  {
    return sol1.better_than_unconstraint(sol2.obj);
  }

  bool solution_t::better_than_unconstraint(const vec_t & other_obj) const
  {

    assert(other_obj.size() == this->obj.size());

    bool strict = false;

    for (size_t i = 0; i < this->obj.size(); ++i)
    {
      // if (std::abs(this->obj[i] - other_obj[i]) >= 0.00001) // not 'equal'
      {
        if (!this->better_than_unconstraint_per_objective(i, other_obj[i]))
        {
          return false;
          // break;
        }
        
        if (this->strictly_better_than_unconstraint_per_objective(i,other_obj[i])) {
          strict = true;
        }

      }
    }

    if (strict == false) {
      return false;
    }

    return true;

  }
  
  // strictly better if the pentaly is strictly better, or,
  // when all objectives are strictly better
  bool solution_t::strictly_better_than_unconstraint(const vec_t & other_obj) const
  {

    for (size_t i = 0; i < this->obj.size(); ++i)
    {
      if (this->strictly_better_than_unconstraint_per_objective(i,other_obj[i]))
      {
        return false;
        break;
      }
    }

    return true;

  }

  bool solution_t::strictly_better_than(const solution_t & other) const
  {

    if (this->penalty > 0) // this is infeasible 
    {
      if (other.penalty > 0) // Both are infeasible
      {
        if (this->penalty < other.penalty)
          return true;
      }
    }
    else // this is feasible
    {
      if (other.penalty > 0) { // this is feasible and other is not
        return true;
      }
      else { // Both are feasible */
        return this->strictly_better_than_unconstraint(other.obj);
      }
    }

    return false;

  }
  
  
  bool solution_t::better_than_per_objective(const solution_t & other, size_t objective_number) const
  {
    if (this->penalty > 0) // this is infeasible
    {
      if (other.penalty > 0) // Both are infeasible
      {
        if (this->penalty < other.penalty)
          return true;
      }
    }
    else // this is feasible
    {
      if (other.penalty > 0) { // this is feasible and other is not
        return true;
      }
      else { // Both are feasible */
        return this->strictly_better_than_unconstraint_per_objective(objective_number, other.obj[objective_number]);
      }
    }
    
    return false;

  }

  bool solution_t::better_than_unconstraint_per_objective(size_t objective_number, double other_obj) const
  {
    assert(this->obj.size() > objective_number);
    return (this->obj[objective_number] <= other_obj);
  }

  bool solution_t::strictly_better_than_unconstraint_per_objective(size_t objective_number, double other_obj) const
  {
    assert(this->obj.size() > objective_number);
    return (this->obj[objective_number] < other_obj);
  }
  
  bool solution_t::same_objectives(const solution_t & other) const
  {
    assert(this->obj.size() == other.obj.size());
    
    bool same_objectives = true;
    for (size_t i = 0; i < number_of_objectives(); i++)
    {
      if (this->obj[i] != other.obj[i])
      {
        same_objectives = false;
        break;
      }
    }

    if (same_objectives && (this->penalty != other.penalty))
    {
      same_objectives = false;
    }

    return same_objectives;
  }


  // computes the distance to another solution
  //------------------------------------------------------------------------------------
  double solution_t::param_distance(solution_t & sol2) const {
    return param_distance(sol2.param);
  }

  double solution_t::param_distance(vec_t & param2) const
  {
    assert(this->param.size() == param2.size());
    return (this->param - param2).norm();
  }


  double solution_t::minimal_objective_distance(const std::vector<solution_pt> & sols, const vec_t & obj_ranges) const
  {
    
    double distance, distance_smallest = -1.0;
    for (size_t j = 0; j < sols.size(); j++)
    {
      if(sols[j] != nullptr)
      {
        distance = transformed_objective_distance(*sols[j], obj_ranges);
        if ((distance_smallest < 0) || (distance < distance_smallest)) {
          distance_smallest = distance;
        }
      }
    }

    if (distance_smallest == -1.0) {
      distance_smallest = 1e308;
    }

    return distance_smallest;

  }

  double solution_t::minimal_param_distance(const std::vector<solution_pt> & sols) const
  {

    double distance, distance_smallest = -1.0;
    for (size_t j = 0; j < sols.size(); j++)
    {
      if(sols[j] != nullptr)
      {
        distance = param_distance(*sols[j]);
        if ((distance_smallest < 0) || (distance < distance_smallest)) {
          distance_smallest = distance;
        }
      }
    }

    if (distance_smallest == -1.0) {
      distance_smallest = 1e308;
    }

    return distance_smallest;

  }

  double solution_t::transformed_objective_distance(const solution_t & other, const vec_t & obj_ranges) const {
    return transformed_objective_distance(other.obj, obj_ranges);
  }

  double solution_t::transformed_objective_distance(const vec_t & other_obj, const vec_t & obj_ranges) const
  {
    
    return this->obj.scaled_euclidean_distance(other_obj, obj_ranges);

  }

}







