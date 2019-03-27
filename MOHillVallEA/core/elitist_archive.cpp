/*
 
 MO-HillVallEA
 
 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "elitist_archive.h"
#include "optimizer.h"
#include "mathfunctions.h"
#include "hillvalleyclustering.h"

hillvallea::elitist_archive_t::elitist_archive_t(size_t target_size, rng_pt rng) : population_t()
{
  this->target_size = target_size;
  this->use_parameter_distances = false;
  this->use_greedy_selection = false;
  this->objective_discretization_in_effect = false;
  this->parameter_discretization_in_effect = false;
  this->rng = rng;
  this->use_niching = false;
  
  this->sols.reserve(1.5*target_size);
}

hillvallea::elitist_archive_t::~elitist_archive_t() {}



hillvallea::elitist_archive_t::elitist_archive_t(const elitist_archive_t & other) : population_t(other)
{
  this->target_size = other.target_size;
  this->use_parameter_distances = other.use_parameter_distances;
  this->use_greedy_selection = other.use_greedy_selection;
  this->use_niching = other.use_niching;
  
  this->objective_discretization_in_effect = other.objective_discretization_in_effect;
  this->parameter_discretization_in_effect = other.parameter_discretization_in_effect;
  this->objective_discretization = other.objective_discretization;
  this->parameter_discretization = other.parameter_discretization;
  this->best_objective_values_in_elitist_archive = other.best_objective_values_in_elitist_archive;
  this->rng = other.rng;
  
  // Copy each of hte clustered archives.
  this->clustered_archive.resize(other.clustered_archive.size());
  
  for(size_t i = 0; i < other.clustered_archive.size(); ++i) {
    this->clustered_archive[i] = std::make_shared<elitist_archive_t>(*other.clustered_archive[i]);
  }

}


/**
* Adds a solution to the elitist archive.
*/
void hillvallea::elitist_archive_t::addToArchive(const solution_pt sol, size_t & insert_index, bool make_copy_of_sol)
{

  assert(sol != nullptr);
  
  // to prevent gaps please.
  if (insert_index > sols.size()) {
    insert_index = sols.size();
  }

  if (insert_index == sols.size())
  {
    
    if(make_copy_of_sol) {
      sols.push_back(std::make_shared<solution_t>(*sol));
      sols.back()->elite_origin = nullptr;
    }
    else {
      sols.push_back(sol);
    }
    
  }
  else
  {
    
    if(make_copy_of_sol) {
      sols[insert_index] = std::make_shared<solution_t>(*sol);
      sols[insert_index]->elite_origin = nullptr;
    }
    else {
      sols[insert_index] = sol;
    }
    
  }
  

  // note, if all solutions have a penalty, the archive consists of the solution with the lowest penalty,
  // and is always size 1.
  if (sols.size() == 1)
  {
    best_objective_values_in_elitist_archive = sol->obj;
  } 
  else
  {
    
    // in this case, there are solutions with zero penalty.
    for (size_t j = 0; j < sol->number_of_objectives(); j++)
    {
      if (sol->strictly_better_than_unconstraint_per_objective(j, best_objective_values_in_elitist_archive[j])) {
        best_objective_values_in_elitist_archive[j] = sol->obj[j];
      }
    }
  }
  
}

void hillvallea::elitist_archive_t::removeFromArchive(size_t sol_index)
{
  sols[sol_index] = nullptr;
}

// Update the Archive
// Updates the elitist archive by offering a new solution to possibly be added to the archive.
// If there are no solutions in the archive yet, the solution is added.
// Otherwise, the number of times the solution is dominated is computed.
//  Solution A is always dominated by solution B that is in the same domination-box if B dominates A or A and B do not dominate each other.

int hillvallea::elitist_archive_t::updateArchive(const solution_pt sol)
{
  return updateArchive(sol, false);
}

int hillvallea::elitist_archive_t::updateArchive(const solution_pt sol, bool make_copy_of_sol)
{
  
  int added_index = -1;
  
  if(sol == nullptr) {
    return added_index;
  }
  
  bool is_extreme_compared_to_archive = false;
  
  size_t insert_index = sols.size();

  if (sol->penalty == 0)
  {
    if (sols.size() == 0) {
      is_extreme_compared_to_archive = true;
    }
    else
    {
      for (size_t i = 0; i < sol->number_of_objectives(); ++i)
      {
        if (sol->strictly_better_than_unconstraint_per_objective(i, best_objective_values_in_elitist_archive[i])) {
          is_extreme_compared_to_archive = true;
          break;
        }
      }
    }
  }

  // if the archive is empty, just insert the solution.
  if (sols.size() == 0)
  {
    addToArchive(sol, insert_index, make_copy_of_sol);
    added_index = (int) insert_index;
  }
  else
  {

    bool is_dominated_itself = false;
    bool all_to_be_removed = true;

    for (size_t i = 0; i < sols.size(); i++)
    {
      // we replace the first 'empty' spot in the archive.
      if (sols[i] == nullptr)
      {
        if (i < insert_index) {
          insert_index = i;
        }

        continue;
      }

      all_to_be_removed = false;
      if (sols[i]->better_than(*sol))
      {
        is_dominated_itself = true;
      }
      else
      {
        if (!sol->better_than(*sols[i]))
        {
          // we want only one solution per objective box
          // so if it is not better than the solution in the box,
          // we still discard it.
          if (!is_extreme_compared_to_archive) {
            if (sameObjectiveBox(*sols[i], *sol))
            {
              is_dominated_itself = true;
            }
          }
        }
      }

      if (is_dominated_itself) {
        break;
      }
    }

    // if all in the elitist archive are nullpointers, insert_index is set to 0 above, and
    // we can add the solution to the archive.
    if (all_to_be_removed) {
      addToArchive(sol, insert_index, make_copy_of_sol);
      added_index = (int) insert_index;
    }
    else
    {
      // don't add dominated solutions.
      if (!is_dominated_itself)
      {
        for (size_t i = 0; i < sols.size(); i++)
        {
          if (sols[i] == nullptr) {
            continue;
          }
          if (sol->better_than(*sols[i]) || (sameObjectiveBox(*sols[i], *sol)))
          {
            // if the to-be-removed solution is the best in the archive, update the best. 
            for (size_t j = 0; j < sol->number_of_objectives(); j++)
            {
              if (sols[i]->obj[j] == best_objective_values_in_elitist_archive[j])
              {
                best_objective_values_in_elitist_archive[j] = sol->obj[j];
              }
            }

            removeFromArchive(i);

          }
        }

        addToArchive(sol, insert_index, make_copy_of_sol);
        added_index = (int) insert_index;
      }

    }
  }
  
  return added_index;
  
}

// returns true if two solutions are in the same discretization box
bool hillvallea::elitist_archive_t::sameObjectiveBox(const solution_t & sol1, const solution_t & sol2) const
{

  // when using parameter distances, the objective discretization is ignored
  if(!use_parameter_distances)
  {
    if (!objective_discretization_in_effect)
    {
      // If the solutions are identical, they are still in the (infinitely small) same objective box.
      for (size_t i = 0; i < sol1.number_of_objectives(); i++)
      {
        if (sol1.obj[i] != sol2.obj[i])
          return false;
      }

      return true;
    }

    for (size_t i = 0; i < sol1.number_of_objectives(); i++)
    {

      int box1 = (int)(sol1.obj[i] / objective_discretization[i]);
      int box2 = (int)(sol2.obj[i] / objective_discretization[i]);

      if (box1 != box2) {
        return false;
      }
    }

    return true;
  }
  else
  {
    return false;
  }
  
}

void hillvallea::elitist_archive_t::clear()
{
  if(use_niching)
  {
    for(size_t c = 0; c < clustered_archive.size(); ++c) {
      if(clustered_archive[c] != nullptr) {
        clustered_archive[c]->clear();
      }
    }
    clustered_archive.clear();
  }
  else
  {
    sols.clear();
  }
}


/**
 * Discard similar solutions when the archive  archive needs to be filtered for
 *
 */
void hillvallea::elitist_archive_t::adaptArchiveSize()
{
  
  if(use_niching)
  {
    // niching
    if(use_greedy_selection) {
      adaptSizeBygreedyScatteredSubsetSelection_mm();
    }
    else {
      adaptObjectiveDiscretization_mm();
    }
  }
  else
  {
    // not niching
    if(use_greedy_selection) {
      adaptSizeBygreedyScatteredSubsetSelection();
    }
    else {
      adaptObjectiveDiscretization();
    }
  }
  
}


/**
* Adapts the objective box discretization. If the numbre of solutions in the elitist archive is too high or too low
* compared to the population size, the objective box discretization is adjusted accordingly. In doing so, the
* entire elitist archive is first emptied and then refilled.
*/
void hillvallea::elitist_archive_t::adaptObjectiveDiscretization()
{
  size_t elitist_archive_size_target_lower_bound = (size_t)(0.75*target_size);
  size_t elitist_archive_size_target_upper_bound = (size_t)(1.25*target_size); // must be at least 1

  // disable the discretization if the archive is too small
  if (sols.size() < elitist_archive_size_target_lower_bound)
  {
    objective_discretization_in_effect = false;
    
    for(size_t c = 0; c < clustered_archive.size(); ++c) {
      clustered_archive[c]->objective_discretization_in_effect = false;
    }
  }

  if(sols.size() > elitist_archive_size_target_upper_bound){
    removeSolutionNullptrs();
  }
  
  // if the archive size crosses the upperbound, adapt the discretization and dump solutions
  if (sols.size() > elitist_archive_size_target_upper_bound)
  {

    // std::cout << "Achive size reduced by objective space discretization." << std::endl;
    vec_t objective_ranges;
    objectiveRanges(objective_ranges);
    
    objective_discretization_in_effect = true;

    int na = 1;
    int nb = (int)pow(2.0, 25.0);
    
    std::vector<solution_pt> archive_copy = sols;
    
    for (size_t k = 0; k < 25; k++)
    {
      int nc = (na + nb) / 2;

      objective_discretization = objective_ranges / ((double)nc);

      // Clear the elitist archive
      sols.clear();

      // Rebuild the elitist archive 
      for (size_t i = 0; i < archive_copy.size(); i++) {
        updateArchive(archive_copy[i]);
      }

      if (sols.size() <= elitist_archive_size_target_lower_bound) {
        na = nc;
      }
      else {
        nb = nc;
      }

    }

  }
}


bool hillvallea::elitist_archive_t::solutionHasImproved(const hillvallea::solution_t & sol) const
{
  
  if (sols.size() == 0)
    return true;
  
  bool result = false;

  // if it is better than the best in a single objective
  if (sol.penalty == 0)
  {
    for (size_t j = 0; j < sol.number_of_objectives(); j++)
    {
      if (sol.strictly_better_than_unconstraint_per_objective(j, best_objective_values_in_elitist_archive[j]))
      {
        result = true;
        break;
      }
    }
  }

  // if not, it must dominate a solution
  if (result != true)
  {
    result = true;
    // for each sol in the archive,
    for (size_t i = 0; i < sols.size(); i++)
    {

      if (sols[i] == nullptr) {
        continue;
      }

      // if the solution in the archive is better than the solution at hand, it has not improved the archive.
      if (sols[i]->better_than(sol))
      {
        result = false;
        break;
      }
      else if (!sol.better_than(*sols[i]))
      {
        if (sameObjectiveBox(*sols[i],sol))
        {
          result = false;
          break;
        }
      }
    }
  }

  return result;
}


void hillvallea::elitist_archive_t::adaptSizeBygreedyScatteredSubsetSelection()
{
  size_t target_upper_bound = (size_t)(1.25*target_size); // must be at least 1
  
  if(size() > target_upper_bound) {
    removeSolutionNullptrs();
  }
  
  if(size() > target_upper_bound)
  {
    
    std::vector<solution_pt> new_sols;
    
    for(size_t i =0; i < sols.size(); ++i)
    {
      if(sols[i] != nullptr) {
        new_sols.push_back(sols[i]);
      }
    }
    
    sols = new_sols; // removes all nullptrs
    
    std::vector<solution_pt> selected_sols;
  
    if(use_parameter_distances) {
      
      // std::cout << "Achive size reduced by greedy parameter subset selection" << std::endl;
      selectSolutionsBasedOnParameterDiversity(sols, target_size, selected_sols, rng);
    }
    else {
      
      // std::cout << "Achive size reduced by greedy objective subset selection" << std::endl;
      vec_t objective_ranges;
      objectiveRanges(objective_ranges);
      
      selectSolutionsBasedOnObjectiveDiversity(sols, target_size, selected_sols, objective_ranges, rng);
    }
    
    sols = selected_sols;

  }
  
}

void hillvallea::elitist_archive_t::removeSolutionNullptrs()
{
  std::vector<solution_pt> new_sols;
  new_sols.reserve(sols.size());
  
  for(size_t i = 0; i < sols.size(); ++i)
  {
    if(sols[i] != nullptr) {
      new_sols.push_back(sols[i]);
    }
  }
  
  sols = new_sols;
}

double hillvallea::elitist_archive_t::get_best_objective_values_in_elitist_archive(size_t obj_index)
{
  
  if(this->size() == 0) {
    return NAN;
  }
  
  // compute best objectives.
  
  if(use_niching)
  {
    solution_t best_objective;
    best_objective.obj = clustered_archive[0]->best_objective_values_in_elitist_archive;
    
    // niching, collect best objectives
    for(size_t i = 0; i < clustered_archive.size(); ++i)
    {
      
      for (size_t j = 0; j < best_objective.number_of_objectives(); j++)
      {
        // I do !strictly better because best_obj is a vector, not a solution, so i cannot call the comparator the other way around
        if (!best_objective.strictly_better_than_unconstraint_per_objective(j, clustered_archive[i]->best_objective_values_in_elitist_archive[j])) {
          best_objective.obj[j] = clustered_archive[i]->best_objective_values_in_elitist_archive[j];
        }
      }
    }
    
    best_objective_values_in_elitist_archive = best_objective.obj;
  }
  else
  {
      solution_t best_objective;
      
    // niching, collect best objectives
    for(size_t i = 0; i < sols.size(); ++i)
    {
      if(sols[i] != nullptr)
      {
        if(best_objective.number_of_objectives() == 0)
        {
          best_objective.obj = sols[i]->obj;
        }
        
        for (size_t j = 0; j < best_objective.number_of_objectives(); j++)
        {
          if (!best_objective.strictly_better_than_unconstraint_per_objective(j, sols[i]->obj[j])) {
            best_objective.obj[j] = sols[i]->obj[j];
          }
        }
      }
 
      best_objective_values_in_elitist_archive = best_objective.obj;
      
    }
  }

  return best_objective_values_in_elitist_archive[obj_index];;

}

hillvallea::elitist_archive_pt hillvallea::elitist_archive_t::initNewArchive() const
{
  elitist_archive_pt new_archive = std::make_shared<elitist_archive_t>(target_size, rng);
  new_archive->target_size = target_size;
  new_archive->use_parameter_distances = use_parameter_distances;
  new_archive->use_greedy_selection = use_greedy_selection;
  new_archive->parameter_discretization_in_effect = parameter_discretization_in_effect;
  new_archive->objective_discretization_in_effect = objective_discretization_in_effect;
  new_archive->parameter_discretization = parameter_discretization;
  new_archive->objective_discretization = objective_discretization;
  new_archive->use_niching = false;
  new_archive->best_objective_values_in_elitist_archive = best_objective_values_in_elitist_archive;
  
  return new_archive;
}



void hillvallea::elitist_archive_t::addArchive(elitist_archive_pt elitist_archive)
{
  clustered_archive.push_back(elitist_archive);
}

void hillvallea::elitist_archive_t::removeEmptyClusters()
{
  
  std::vector<elitist_archive_pt> temp;
  
  for(size_t i =0 ; i < clustered_archive.size(); ++i)
  {
    if(clustered_archive[i] != nullptr && clustered_archive[i]->sols.size() != 0)
    {
      temp.push_back(clustered_archive[i]);
    }
  }
  
  clustered_archive = temp;
  
}


size_t hillvallea::elitist_archive_t::size() const {
  
  if(!use_niching)
  {
    return sols.size();
  }
  else
  {
    size_t total_size = 0;
    
    for(size_t i =0 ; i < clustered_archive.size(); ++i)
      total_size += clustered_archive[i]->size();
    
    return total_size;
  }
  
  
}

size_t hillvallea::elitist_archive_t::actualSize() const {
  
  size_t n = 0;
  if(!use_niching)
  {
    for(size_t i = 0 ; i < sols.size(); ++i) {
      if(sols[i] != nullptr) {
        n++;
      }
    }
  }
  else
  {

    for(size_t i =0 ; i < clustered_archive.size(); ++i) {
      for(size_t j = 0; j < clustered_archive[i]->size(); ++j) {
        if(clustered_archive[i]->sols[j] != nullptr) {
          n++;
        }
      }
    }
  }
  
  return n;
  
  
}

size_t hillvallea::elitist_archive_t::number_of_clusters() const
{
  return clustered_archive.size();
}





// multi-modal version of the objectivespace discretization
// we use the same discretization for all archives
void hillvallea::elitist_archive_t::adaptObjectiveDiscretization_mm()
{
  size_t elitist_archive_size_target_lower_bound = (size_t)(0.75*target_size);
  size_t elitist_archive_size_target_upper_bound = (size_t)(1.25*target_size); // must be at least 1
  
  // Get All Sols
  size_t archive_size = actualSize();
  
  // disable the discretization if the archive is too small
  if (archive_size < elitist_archive_size_target_lower_bound)
  {
    objective_discretization_in_effect = false;
    
    for(size_t i = 0; i < clustered_archive.size(); ++i) {
      clustered_archive[i]->objective_discretization_in_effect = false;
    }
  }
  
  // if the archive size crosses the upperbound, adapt the discretization and dump solutions
  if (archive_size > elitist_archive_size_target_upper_bound)
  {

    vec_t objective_ranges;
    objectiveRanges(objective_ranges);

    // per cluster, set the discretization
    for(size_t i = 0; i < clustered_archive.size(); ++i) {
      clustered_archive[i]->objective_discretization_in_effect = true;
    }
    
    objective_discretization_in_effect = true;
    
    int na = 1;
    int nb = (int)pow(2.0, 25.0);
    
    std::vector<std::vector<solution_pt>> archive_copy(clustered_archive.size());
    
    for(size_t c = 0; c < clustered_archive.size(); ++c) {
      archive_copy[c] = clustered_archive[c]->sols;
    };
    
    for (size_t k = 0; k < 25; k++)
    {
      int nc = (na + nb) / 2;
      
      objective_discretization = objective_ranges / ((double)nc);
      
      for(size_t c = 0; c < clustered_archive.size(); ++c)
      {

        // set discretization
        clustered_archive[c]->objective_discretization = objective_discretization;
        clustered_archive[c]->sols.clear();
        
        // Rebuild the elitist archive
        for (size_t i = 0; i < archive_copy[c].size(); i++) {
          clustered_archive[c]->updateArchive(archive_copy[c][i]);
        }
      }
      
      archive_size = actualSize();
      
      // std::cout << archive_size << ",";
      
      if (archive_size <= elitist_archive_size_target_lower_bound) {
        na = nc;
      }
      else {
        nb = nc;
      }
      
      
    }
    
    
    // std::cout << std::endl;

    
  }
  
}

void hillvallea::elitist_archive_t::adaptSizeBygreedyScatteredSubsetSelection_mm()
{
  size_t target_upper_bound = (size_t)(1.25*target_size);
  
  if(size() > target_upper_bound){
    removeSolutionNullptrs();
  }
  
  if(size() > target_upper_bound)
  {
    
    population_t new_sols;
    std::vector<elitist_archive_pt> sol_origin;
    getAllSols(new_sols.sols, sol_origin);
    
    if(new_sols.sols.size() > target_upper_bound)
    {
      
      std::vector<size_t> selected_indices;
    
      if(use_parameter_distances)
      {
        // convert to the right format
        std::vector<vec_t> parameters(new_sols.sols.size());
        
        for (size_t i = 0; i < new_sols.sols.size(); ++i) {
          parameters[i] = new_sols.sols[i]->param;
        }
        
        greedyScatteredSubsetSelection(parameters, (int)(0.75*target_size), selected_indices, rng);

      }
      else
      {
        vec_t objective_ranges;
        objectiveRanges(objective_ranges);
        
        // we scale the objectives to the objective ranges
        // before performing subset selection
        std::vector<vec_t> scaled_objectives(new_sols.sols.size());

        for (size_t i = 0; i < new_sols.sols.size(); ++i)
        {
          scaled_objectives[i].resize(new_sols.sols[i]->number_of_objectives());
          
          for(size_t j =0; j < new_sols.sols[i]->number_of_objectives(); ++j) {
            scaled_objectives[i][j] = new_sols.sols[i]->obj[j] / objective_ranges[j];
          }
          
        }
        
        greedyScatteredSubsetSelection(scaled_objectives, (int)target_size, selected_indices, rng);

      }
      

      for(size_t j = 0; j < clustered_archive.size(); ++j)
      {
        clustered_archive[j]->sols.clear(); // withing all nullptr's.
      }
        
      
      for(size_t i = 0; i < selected_indices.size(); ++i) {
        sol_origin[selected_indices[i]]->sols.push_back(new_sols.sols[selected_indices[i]]);
      }
    }
  }
  
}



/**
 * Elitism: copies at most 1/k*tau*n solutions per cluster
 * from the elitist archive.
 */
void hillvallea::elitist_archive_t::copyElitesToClusters(std::vector<cluster_pt> & clusters, size_t max_number_of_elites, const vec_t & objective_ranges, rng_pt rng) const
{
  
  for (size_t i = 0; i < clusters.size(); i++) {
    clusters[i]->average_fitness(clusters[i]->objective_mean);
    clusters[i]->elites.clear();
    clusters[i]->elites.reserve(2* max_number_of_elites);
  }
  
  double distance, distance_smallest;
  
  // divide elites over the clusters
  for (size_t i = 0; i < sols.size(); ++i)
  {
    
    if (sols[i] == nullptr) {
      continue;
    }
    
    distance_smallest = -1;
    size_t j_min = 0;
    for (size_t j = 0; j < clusters.size(); ++j)
    {
      
      distance = sols[i]->transformed_objective_distance(clusters[j]->objective_mean, objective_ranges);
      if ((distance_smallest < 0) || (distance < distance_smallest))
      {
        j_min = j;
        distance_smallest = distance;
      }
      
    }
    
    clusters[j_min]->elites.push_back(sols[i]);
    
    
  }
  
  // if there are more than 'max' elites, do a diversity based selection.
  for (size_t i = 0; i < clusters.size(); i++)
  {
    if (clusters[i]->elites.size() > max_number_of_elites)
    {
      std::vector<solution_pt> selected_elites;
      selectSolutionsBasedOnObjectiveDiversity(clusters[i]->elites, max_number_of_elites, selected_elites,  objective_ranges, rng);
      
      clusters[i]->elites = selected_elites;
      
    }
  }
  
}



/**
 * Elitism: copies at most 1/k*tau*n solutions per cluster
 * from the elitist archive.
 */
size_t hillvallea::elitist_archive_t::addElitesToPopulation(population_t & population, int max_number_of_elites)
{
  
  size_t popsize = population.size();
  
  if(max_number_of_elites < 0)
  {
    if(use_niching)
    {
      for (size_t i = 0; i < clustered_archive.size(); i++)
      {
        size_t old_popsize = population.size();
        population.addCopyOfSolutions(clustered_archive[i]->sols);
        
        for(size_t j = old_popsize; j < population.size(); ++j) {
          population.sols[j]->elite_origin = clustered_archive[i];
        }
      }
    }
    else
    {
      population.addCopyOfSolutions(sols);
    }
  }
  else
  {
    if(use_niching)
    {
      // if there are more than 'max' elites, do a diversity based selection.
      size_t max_number_of_elites_per_archive = std::max((size_t) 1, max_number_of_elites / clustered_archive.size());
      
      // Add elites to the population
      vec_t objective_ranges;
      objectiveRanges(objective_ranges);
      
      for (size_t i = 0; i < clustered_archive.size(); i++)
      {
        std::vector<solution_pt> selected_elites;
        selectSolutionsBasedOnObjectiveDiversity(clustered_archive[i]->sols, max_number_of_elites_per_archive, selected_elites, objective_ranges, rng);
        
        size_t init_popsize = population.size();
        population.addCopyOfSolutions(selected_elites);
        
        for(size_t j = init_popsize; j < init_popsize + selected_elites.size(); ++j) {
          population.sols[j]->elite_origin = clustered_archive[i];
        }
        
      }
    } // end use_niching
    else
    {
      std::vector<solution_pt> selected_elites;
      // computeObjectiveRanges();
      // selectSolutionsBasedParameterDiversity(sols, max_number_of_elites, selected_elites, rng);
      selectSolutionsBasedOnParameterDiversity(sols, max_number_of_elites, selected_elites, rng);
      
      population.addCopyOfSolutions(selected_elites);
    }
  }
  
  return (population.size() - popsize);
  
}

void hillvallea::elitist_archive_t::getAllSols(std::vector<solution_pt> & all_sols)
{
  
  all_sols.clear();
  
  if(use_niching)
  {
  
  for(size_t j = 0; j < clustered_archive.size(); ++j)
  {
    all_sols.reserve(all_sols.size() + clustered_archive.size());
    
    std::vector<solution_pt> temp_new_sols;
    temp_new_sols.reserve(clustered_archive[j]->sols.size());
    
    for(size_t i =0; i < clustered_archive[j]->sols.size(); ++i)
    {
      if(clustered_archive[j]->sols[i] != nullptr)
      {
        temp_new_sols.push_back(clustered_archive[j]->sols[i]);
        clustered_archive[j]->sols[i]->cluster_number = (int) j;
      }
    }
    
    clustered_archive[j]->sols = temp_new_sols; // removing all nullptr's.
    
    for(size_t i = 0; i < temp_new_sols.size(); ++i)  {
      all_sols.push_back(temp_new_sols[i]);
    }
  }
  }
  else
  {
    all_sols.reserve(sols.size());
    
    for(size_t i =0; i < sols.size(); ++i)
    {
      if(sols[i] != nullptr)
      {
        all_sols.push_back(sols[i]);
        sols[i]->cluster_number = 0;
      }
    }
    
    sols = all_sols; // removing all nullptr's.

  }

}


void hillvallea::elitist_archive_t::getAllSols(std::vector<solution_pt> & all_sols, std::vector<elitist_archive_pt> & origin)
{
  
  all_sols.clear();
  origin.clear();
  
  if(use_niching)
  {
    
    for(size_t j = 0; j < clustered_archive.size(); ++j)
    {
      all_sols.reserve(all_sols.size() + clustered_archive.size());
      origin.reserve(origin.size() + clustered_archive.size());
      
      std::vector<solution_pt> temp_new_sols;
      temp_new_sols.reserve(clustered_archive[j]->sols.size());
      
      for(size_t i =0; i < clustered_archive[j]->sols.size(); ++i)
      {
        if(clustered_archive[j]->sols[i] != nullptr)
        {
          temp_new_sols.push_back(clustered_archive[j]->sols[i]);
          origin.push_back(clustered_archive[j]);
        }
      }
      
      clustered_archive[j]->sols = temp_new_sols; // removing all nullptr's.
      
      for(size_t i = 0; i < temp_new_sols.size(); ++i)  {
        all_sols.push_back(temp_new_sols[i]);
      }
    }
  }
  else
  {
    all_sols.reserve(sols.size());
    origin.reserve(sols.size());
    
    for(size_t i =0; i < sols.size(); ++i)
    {
      if(sols[i] != nullptr)
      {
        all_sols.push_back(sols[i]);
        //origin.push_back(*this); // todo?
        assert(false);
      }
    }
    
    sols = all_sols; // removing all nullptr's.
    
  }
  
}

void hillvallea::elitist_archive_t::set_use_parameter_distances(bool value)
{
  use_parameter_distances = value;
  
  if(use_niching)
  {
    for(size_t i =0 ; i < clustered_archive.size(); ++i) {
      clustered_archive[i]->set_use_parameter_distances(value);
    }
  }
}

void hillvallea::elitist_archive_t::set_use_greedy_selection(bool value)
{
  use_greedy_selection = value;
  
  if(use_niching)
  {
    for(size_t i =0 ; i < clustered_archive.size(); ++i) {
      clustered_archive[i]->set_use_greedy_selection(value);
    }
  }
}


void hillvallea::elitist_archive_t::computeApproximationSet(size_t approximation_set_size, std::vector<optimizer_pt> & optimizers, elitist_archive_pt & elitist_archive, bool use_parameter_space_diversity, bool terminate_pops)
{

  if(optimizers.size() == 0) {
    return;
  }
  
  if (elitist_archive->use_niching == false || elitist_archive->target_size <= approximation_set_size)
  {
    this->sols.clear();
    
    elitist_archive->getAllSols(this->sols);

    // see if adaptation is required.
    this->target_size = approximation_set_size;
    this->set_use_greedy_selection(true);
    this->set_use_parameter_distances(use_parameter_space_diversity);
    this->adaptArchiveSize();
    
  }
  else
  {
    // use niching and approximation set is smaller than elitist archive.
    // first, we get rid of the clusters that have no rank-0 solution.
    this->sols.clear();
    
    //std::vector<solution_pt> all_sols;
    // elitist_archive->getAllSols(all_sols);
    //for(size_t i =0; i < all_sols.size(); ++i) {
    //  this->updateArchive(all_sols[i]);
    //}
    
    for(size_t i = 0; i < elitist_archive->clustered_archive.size(); ++i)
    {
      for(size_t j = 0; j < elitist_archive->clustered_archive[i]->sols.size(); ++j)
      {
        if(elitist_archive->clustered_archive[i]->sols[j] != nullptr)
        {
          elitist_archive->clustered_archive[i]->sols[j]->cluster_number = (int) i;
          this->updateArchive(elitist_archive->clustered_archive[i]->sols[j]);
        }
      }
    }
    
    std::vector<bool> subpop_with_rank0(elitist_archive->clustered_archive.size(),false);
    
    for(size_t i = 0; i < this->sols.size(); ++i)
    {
      if(this->sols[i] != nullptr) {
      subpop_with_rank0[this->sols[i]->cluster_number] = true; // not sure if cluster number is set properly..
      }
    }
    
    this->sols.clear();
    
    for(size_t i = 0; i < elitist_archive->clustered_archive.size(); ++i) {
      
      if(subpop_with_rank0[i]) {
        this->addSolutions(*elitist_archive->clustered_archive[i]);
      }
    }
    
    // see if adaptation is required.
    this->target_size = approximation_set_size;
    this->set_use_greedy_selection(true);
    this->set_use_parameter_distances(use_parameter_space_diversity);
    this->adaptArchiveSize();
    
  }
  
  
  
  // Terminate populations!
  // size_t archive_size = this->actualSize();
  
  if(terminate_pops && optimizers.size() > 1)
  {
    
    size_t elites_added = 0;
    
    for(size_t i =0; i < optimizers.size(); ++i)
    {
      if(!optimizers[i]->terminated) {
        elites_added += optimizers[i]->new_elites_added;
      }
    }
    
    for(size_t i =0; i < optimizers.size() -1; ++i)
    {
      // if(optimizers[i]->new_elites_added < 0.5 * optimizers[i+1]->new_elites_added) {
      if(!optimizers[i]->terminated && optimizers[i]->new_elites_added < 0.25 * elites_added) {
        optimizers[i]->terminated = true; std::cout << "-";
      }
    }
  }
    
}

size_t hillvallea::elitist_archive_t::collectSubarchives(std::vector<population_pt> & subpopulations)
{
 
  initArchiveForEachSubpop(subpopulations);
  
  // Create an Elitist archive from all subarchives
  //---------------------------------------------------------------------
  size_t number_of_solutions_added = 0;
  clear();

  for(size_t i = 0; i < subpopulations.size(); ++i)
  {
    addArchive(subpopulations[i]->elitist_archive);
    number_of_solutions_added += subpopulations[i]->new_elites_added;
  }

  return number_of_solutions_added;
}

// creates an archive per subpop
void hillvallea::elitist_archive_t::initArchiveForEachSubpop(std::vector<population_pt> & subpopulations) const
{

  for(size_t i = 0; i < subpopulations.size(); ++i)
  {
    // create new archive
    subpopulations[i]->elitist_archive = initNewArchive();
    subpopulations[i]->new_elites_added = 0;
    
    // first, add all previously sampled elites (to reconstruct the archive).
    for(size_t j = 0; j < subpopulations[i]->size(); ++j)
    {
      // add the solution to the archive, and count how many non-elites are added (to measure performance of this optimizer)
      if(subpopulations[i]->sols[j]->elite_origin != nullptr) {
        subpopulations[i]->elitist_archive->updateArchive(subpopulations[i]->sols[j], true);
      }
    }
    
    // then, add all novel solutions.
    for(size_t j = 0; j < subpopulations[i]->size(); ++j)
    {
      // add the solution to the archive, and count how many non-elites are added (to measure performance of this optimizer)
      if(subpopulations[i]->sols[j]->elite_origin == nullptr) {
        if (subpopulations[i]->elitist_archive->updateArchive(subpopulations[i]->sols[j], true) >= 0) {
          subpopulations[i]->new_elites_added++;
        }
      }
    }
  }
}

void hillvallea::elitist_archive_t::objectiveRanges(vec_t & objective_ranges)
{
  
  if(use_niching)
  {
    std::vector<solution_pt> backup_sols = this->sols;
    
    getAllSols(this->sols);
    population_t::objectiveRanges(objective_ranges);
    
    this->sols = backup_sols;
  }
  else
  {
    population_t::objectiveRanges(objective_ranges);
  }
  
}















