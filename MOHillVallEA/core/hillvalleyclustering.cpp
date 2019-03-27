

/*
 
 MO-HillVallEA

 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "hillvalleyclustering.h"
#include "mathfunctions.h"
#include "cluster.h"

// c++ rule of three (two)
//----------------------------------------------------------------------------------
hillvallea::hvc_t::hvc_t(fitness_pt fitness_function)
{
  this->fitness_function = fitness_function;
}

hillvallea::hvc_t::~hvc_t() {}

// returns true if the solutions belong to the same niche
bool hillvallea::hvc_t::hill_valley_test(const hillvallea::solution_t &sol1, const hillvallea::solution_t &sol2, const unsigned int max_trials, unsigned int & number_of_evaluations, std::vector<solution_pt> & test_solutions)
{
 
  for (size_t k = 0; k < max_trials; k++)
  {

    solution_pt test = std::make_shared<solution_t>(sol1.number_of_parameters(), sol1.number_of_objectives());

    test->param = sol1.param + ((k + 1.0) / (max_trials + 1.0)) * (sol2.param - sol1.param);
    number_of_evaluations++;
    fitness_function->evaluate(test);

    test_solutions.push_back(test);
    
    // if the solution belongs to a different valley in any objective, 
    // the solutions do not belong to the same niche
    for (size_t objective_number = 0; objective_number < sol1.number_of_objectives(); ++objective_number) 
    {
      if (hill_valley_test_point_single_objective(sol1, sol2, *test, objective_number) == false) {
        return false;
      }
    }

  }

  return true;

}

// hill-valley test in a single objective dimension
// returns true if two solutions belong to the same basin
bool hillvallea::hvc_t::hill_valley_test_single_objective(const hillvallea::solution_t &sol1, const hillvallea::solution_t &sol2, size_t objective_number, unsigned int max_trials, unsigned int & number_of_evaluations, std::vector<solution_pt> & test_solutions)
{
  
  // sample max_trials test points along the line segment between sol1 and sol2
  for (size_t k = 0; k < max_trials; k++)
  {
    solution_pt test = std::make_shared<solution_t>(sol1.number_of_parameters(), sol1.number_of_objectives());
    test->param = sol1.param + ((k + 1.0) / (max_trials + 1.0)) * (sol2.param - sol1.param);
    number_of_evaluations++;
    fitness_function->evaluate(test);
    test_solutions.push_back(test);

    // if the solution belongs to a different valley in the given objective_number, 
    // the solutions do not belong to the same niche
    if (hill_valley_test_point_single_objective(sol1,sol2,*test,objective_number) == false) {
      return false;
    }

  }

  return true;

}



// hill-valley test in a single objective dimension
// returns true if two solutions belong to the same basin
bool hillvallea::hvc_t::hill_valley_test_single_objective_reuse(const hillvallea::solution_t &sol1, const hillvallea::solution_t &sol2, size_t objective_number, unsigned int max_trials, unsigned int & number_of_evaluations, std::vector<solution_pt> & test_solutions)
{
  
  // if 1 is the better sol, then 2 is the worst, so take that value
  double worst_objective;
  double best_objective;
  bool reverse_direction = false; // default, sample from sol1 to sol2.
  
  if (sol1.better_than_unconstraint_per_objective(objective_number, sol2.obj[objective_number])) {
    worst_objective = sol2.obj[objective_number];
    best_objective = sol1.obj[objective_number];
    reverse_direction = true;
  }
  else {
    worst_objective = sol1.obj[objective_number];
    best_objective = sol2.obj[objective_number];
  }
  
  
  // sample all solutions
  // sample max_trials test points along the line segment between sol1 and sol2
  if(test_solutions.size() == 0)
  {
    test_solutions.resize(max_trials);
    
    for (size_t k = 0; k < max_trials; k++)
    {
      test_solutions[k] = std::make_shared<solution_t>();
      test_solutions[k]->param = sol1.param + ((k + 1.0) / (max_trials + 1.0)) * (sol2.param - sol1.param);
    }
  }
  
  // try the points, and evaluate only if the objective is not set yet.
  size_t current_index;
  for (size_t k = 0; k < max_trials; k++)
  {
    
    if(reverse_direction) {
      current_index = max_trials - k - 1;
    } else {
      current_index = k;
    }
      
    
    if(test_solutions[current_index]->number_of_objectives() == 0) {
      number_of_evaluations++;
      fitness_function->evaluate(test_solutions[current_index]);
    }
    
    // if the solution belongs to a different valley in the given objective_number,
    // the solutions do not belong to the same niche
    if (hill_valley_test_point_single_objective(worst_objective,*test_solutions[current_index],objective_number) == false) {
      return false;
    }
    
    // if the solution is better, update the worst. But not if it is better than the best of sol1,sol2.
    //if(!test_solutions[current_index]->better_than_unconstraint_per_objective(objective_number, best_objective)) {
    //  worst_objective = test_solutions[current_index]->obj[objective_number];
    //}
  }
  
  return true;
  
}

bool hillvallea::hvc_t::hill_valley_test_point_single_objective(double worst_objective_value, const hillvallea::solution_t &test, size_t objective_number)
{
  
  // if 'test' is better than the worst_objective, accept the edge,
  // if 'test' is worse, reject the edge, there is a hill in between two valleys (minimization)
  if (test.better_than_unconstraint_per_objective(objective_number, worst_objective_value)) {
    return true;
  }
  else {
    return false;
  }
}


bool hillvallea::hvc_t::hill_valley_test_point_single_objective(const hillvallea::solution_t &sol1, const hillvallea::solution_t &sol2, const hillvallea::solution_t &test, size_t objective_number)
{
  solution_t worst(sol1.number_of_parameters(), sol1.number_of_objectives());
  double worst_objective;

  // if 1 is the better sol, then 2 is the worst, so take that value
  if (sol1.better_than_unconstraint_per_objective(objective_number, sol2.obj[objective_number])) {
    worst_objective = sol2.obj[objective_number];
  }
  else {
    worst_objective = sol1.obj[objective_number];
  }

  // if 'test' is better than the worst_objective, accept the edge, 
  // if 'test' is worse, reject the edge, there is a hill in between two valleys (minimization)
  if (test.better_than_unconstraint_per_objective(objective_number, worst_objective)) {
    return true;
  }
  else {
    return false;
  }
}

// this is the SO code, which assumes that the pop is sorted on fitness.
void hillvallea::hvc_t::cluster(const population_t & pop, std::vector<population_pt> & subpopulations, unsigned int  & number_of_evaluations, double & average_edge_length, bool add_test_solutions, bool recheck_elites, int optimizer_number, rng_pt & rng)
{

  // nothing to cluster
  if (pop.size() == 0) {
    return;
  }
  
  // set problem parameters
  size_t initial_popsize = pop.size();
  size_t number_of_parameters = pop.problem_size();
  size_t number_of_objectives = pop.sols[0]->number_of_objectives();
  
  // Hill-Valley Parameters
  vec_t param_ranges;
  pop.computeParameterRanges(param_ranges);
  double search_volume = param_ranges.prod();
  
  average_edge_length = pow(search_volume / initial_popsize, 1.0 / number_of_parameters);

  // such that for d={1,2} still Nn = d+1
  size_t max_number_of_neighbours;
  
  if( number_of_parameters <= 2 ) {
    max_number_of_neighbours = number_of_parameters + 1;
  }
  else {
    max_number_of_neighbours = 2 + log(number_of_parameters);
  }
  
  
  vec_t dist(initial_popsize,0.0);
  
  bool edge_added = false;
  unsigned int max_number_of_trial_solutions;
  
  // number clusters
  vec_t number_of_clusters(number_of_objectives, 1);
  std::vector<vec_t> cluster_index(initial_popsize);
  
  for(size_t i =0 ; i < initial_popsize; ++i) {
    cluster_index[i].resize(number_of_objectives);
  }
  
  // adding test solutions to the clusters
  // we store all test solutions for sol i and sol j in test_solutions[i][j] with i < j.
  size_t test_solution_hash, mini, maxi; // i * initial_popsize + j for (i > j)
  std::map<size_t,std::vector<solution_pt>> test_solution_hashes;
  /* std::vector<std::vector<std::vector<solution_pt>>> test_solutions(initial_popsize);
  for(size_t i = 0 ; i < initial_popsize; ++i) {
    test_solutions[i].resize(initial_popsize);
  } */

  for(size_t objective_number = 0; objective_number < number_of_objectives; ++objective_number)
  {
    
    // compute population ranks for current objective_number
    std::vector<size_t> sols_sorted_on_fitness(initial_popsize);
    vec_t collected_fitness(initial_popsize);
    
    for(size_t i = 0; i < initial_popsize; ++i)
    {
      collected_fitness[i] = pop.sols[i]->obj[objective_number] + pop.sols[i]->penalty * 1000;
    }
    
    compute_ranks_asc(collected_fitness, sols_sorted_on_fitness); // oops, i hardcoded that we do minimization. assumed
    
    // pop.getSingleObjectiveRanks(fitness_ranks);
    
    // The best is the first cluster.
    cluster_index[sols_sorted_on_fitness[0]][objective_number] = 0;
    
    for (size_t i = 1; i < initial_popsize; i++)
    {
      
      // compute the distance to all better solutions.
      dist[i] = 0.0;
      size_t nearest_better_index = 0;
      size_t furthest_better_index = 0;
      size_t old_nearest_better_index = 0;
      
      for (size_t j = 0; j < i; j++) {
        dist[j] = pop.sols[sols_sorted_on_fitness[i]]->param_distance(*pop.sols[sols_sorted_on_fitness[j]]);
        
        if (dist[j] < dist[nearest_better_index]) {
          nearest_better_index = j;
        }
        
        if (dist[j] > dist[furthest_better_index]) {
          furthest_better_index = j;
        }
      }
      
      edge_added = false;
      
      // check the few nearest solutions with better fitness
      for (size_t j = 0; j < std::min(i, max_number_of_neighbours); j++)
      {
        
        // find the next-to nearest index
        if (j > 0)
        {
          old_nearest_better_index = nearest_better_index;
          nearest_better_index = furthest_better_index;
          
          for (size_t k = 0; k < i; k++)
          {
            if (dist[k] > dist[old_nearest_better_index] && dist[k] < dist[nearest_better_index]) {
              nearest_better_index = k;
            }
          }
          
        }
        
        if (!recheck_elites)
        {
          // if both are elites that belong to the same archive, accept edge.
          // if both are elites that belong to different archives, reject edge.
          
          if( pop.sols[sols_sorted_on_fitness[i]]->elite_origin != nullptr && pop.sols[sols_sorted_on_fitness[nearest_better_index]]->elite_origin != nullptr)
          {
            
            if(pop.sols[sols_sorted_on_fitness[i]]->elite_origin ==  pop.sols[sols_sorted_on_fitness[nearest_better_index]]->elite_origin)
            {
              // accept edge
              cluster_index[sols_sorted_on_fitness[i]][objective_number] = cluster_index[sols_sorted_on_fitness[nearest_better_index]][objective_number];
              edge_added = true;
              break;
            }
            else
            {
              // reject edge
              // double check them, so that we combine more clusters.
              // continue;
            }
          }
        }
        
        // set budget
        max_number_of_trial_solutions = 1 + ((unsigned int)(dist[nearest_better_index] / average_edge_length));
        
        // such that mini < maxi
        mini = std::min(sols_sorted_on_fitness[nearest_better_index], sols_sorted_on_fitness[i]);
        maxi = std::max(sols_sorted_on_fitness[nearest_better_index], sols_sorted_on_fitness[i]);
        test_solution_hash = maxi * initial_popsize + mini;
        
        // there are no pre-saved test_solutions yet
        if(test_solution_hashes.find(test_solution_hash) == test_solution_hashes.end()) {
          test_solution_hashes[test_solution_hash] = {};
        }
        
        if (hill_valley_test_single_objective_reuse(*pop.sols[mini], *pop.sols[maxi], objective_number, max_number_of_trial_solutions, number_of_evaluations, test_solution_hashes[test_solution_hash]))
        {
          cluster_index[sols_sorted_on_fitness[i]][objective_number] = cluster_index[sols_sorted_on_fitness[nearest_better_index]][objective_number];
          edge_added = true;
          break;
        }
      }
      
      
      // its a new clusters, label it like that.
      if (!edge_added)
      {
        cluster_index[sols_sorted_on_fitness[i]][objective_number] = number_of_clusters[objective_number];
        number_of_clusters[objective_number]++;
      }
    }
  }
  
  // create & fill the clusters
  //---------------------------------------------------------------------------
  // map cluster indices:
  // n0 = 0,...,N0 - 1
  // n1 = 0,...,N1 - 1
  // n2 = 0,...,N2 - 1
  // current_cluster_index = N0*N1*n2 + N0*n1 + n0
  
  vec_t terms(number_of_objectives,1.0);
  
  for(size_t j = 0; j < number_of_objectives; ++j)
  {
    for(size_t subj = 0; subj < j; ++subj) {
      terms[j] *= number_of_clusters[subj];
    }
    
  }
  
  std::vector<population_pt> new_clusters;
  std::map<int, size_t> cluster_hashes;
  
  int current_hash;
  for (size_t i = 0; i < initial_popsize; ++i)
  {
    current_hash = (int) terms.dot(cluster_index[i]);
    
    // if the hash is not found yet, add it and create corresponding new cluster.
    if(cluster_hashes.find(current_hash) == cluster_hashes.end())
    {
      cluster_hashes[current_hash] = new_clusters.size();
      new_clusters.push_back(std::make_shared<population_t>());
    }
    
    new_clusters[cluster_hashes[current_hash]]->sols.push_back(pop.sols[i]);
  }
  
  // now i want to add the test solutions to the clusters.
  if( add_test_solutions )
  {
    int clusteri, clusterj;
    for(size_t i = 0; i < initial_popsize; ++i)
    {
      clusteri = terms.dot(cluster_index[i]);
      
      for(size_t j = i + 1; j < initial_popsize; ++j)
      {
        clusterj = terms.dot(cluster_index[j]);
        
        if(clusteri == clusterj) // only if the two solutions belong to the same cluster
        {
          test_solution_hash = j * initial_popsize + i;
          for(size_t k = 0; k < test_solution_hashes[test_solution_hash].size(); ++k)
          {
            if(test_solution_hashes[test_solution_hash][k]->number_of_objectives() > 0) // only if it is actually evaluated
            {
              new_clusters[cluster_hashes[clusteri]]->sols.push_back(test_solution_hashes[test_solution_hash][k]);
              new_clusters[cluster_hashes[clusteri]]->sols.back()->population_number = optimizer_number;
            }
          }
        }
      }
    }
  }
  
  
  // add only the clusters with sols
  for(size_t i = 0; i < new_clusters.size(); ++i)
  {
    subpopulations.push_back(new_clusters[i]);
    int cluster_number = (int) (subpopulations.size() - 1);
    
    for(size_t j = 0; j < subpopulations[cluster_number]->sols.size(); ++j) {
      subpopulations.back()->sols[j]->cluster_number = cluster_number;
    }
  }
  
  
}

void hillvallea::hvc_t::cluster_mo_distances(population_t & pop, std::vector<population_pt> & subpopulations, unsigned int  & number_of_evaluations, double & average_edge_length, bool add_test_solutions, bool recheck_elites, int optimizer_number, rng_pt & rng)
{
  
  if (pop.size() == 0) {
    return;
  }
  
  // Get problem parameters
  size_t initial_popsize = pop.size();
  size_t number_of_parameters = pop.problem_size();
  size_t number_of_objectives = pop.sols[0]->number_of_objectives();

  
  // Compute parameter range and expected edge length
  vec_t param_ranges;
  pop.computeParameterRanges(param_ranges);
  double search_volume = param_ranges.prod();
  average_edge_length = pow(search_volume / initial_popsize, 1.0 / number_of_parameters);
  
  vec_t objective_ranges;
  pop.objectiveRanges(objective_ranges);
  
  // Set clustering parameters
  size_t max_number_of_neighbours = (size_t) (2 + log((double) number_of_objectives));
  
  // cluster index, one index per objective for each solution
  std::vector<vec_t> cluster_index(initial_popsize);

  for(size_t i =0 ; i < initial_popsize; ++i) {
    cluster_index[i].resize(number_of_objectives);
  }
  
  // Start Clustering
  // vec_t param_dist(initial_popsize,0.0);
  vec_t objective_dist(initial_popsize, 0.0);
  bool edge_added = false;
  unsigned int max_number_of_trial_solutions;
  vec_t number_of_clusters(number_of_objectives, 1);
  
  // store the test solutions by a hash.
  size_t test_solution_hash, mini, maxi; // i * initial_popsize + j for (i > j)
  std::map<size_t,std::vector<solution_pt>> test_solution_hashes;

  
  for(size_t objective_number = 0; objective_number < number_of_objectives; ++objective_number)
  {
    
    //--------------------------------------------------------------------
    // compute population ranks for current objective_number
    
    std::vector<size_t> sols_sorted_on_fitness(initial_popsize);
    vec_t collected_fitness(initial_popsize);
    
    // this is one ugly implementation, how the penalties are added and stuff.
    // but i guess it works.
    // oops, i hardcoded that we do minimization. assumed
    for(size_t i = 0; i < initial_popsize; ++i) {
      collected_fitness[i] = pop.sols[i]->obj[objective_number] + pop.sols[i]->penalty * 1000;
    }
    
    compute_ranks_asc(collected_fitness, sols_sorted_on_fitness);
    
    //--------------------------------------------------------------------
    // The best is the first cluster. + iterate for the rest.
    cluster_index[sols_sorted_on_fitness[0]][objective_number] = 0;
    
    for (size_t i = 1; i < initial_popsize; i++)
    {
      edge_added = false;
      
      // compute the distance to all better solutions.
      // param_dist[i] = 0.0;
      objective_dist[i] = 0.0;
      size_t nearest_better_index = 0;
      size_t furthest_better_index = 0;
      size_t old_nearest_better_index = 0;
      
      for (size_t j = 0; j < i; j++) {
        // param_dist[j] = pop.sols[sols_sorted_on_fitness[i]]->param_distance(*pop.sols[sols_sorted_on_fitness[j]]);
        objective_dist[j] = pop.sols[sols_sorted_on_fitness[i]]->transformed_objective_distance(*pop.sols[sols_sorted_on_fitness[j]], objective_ranges);
        
        // if (param_dist[j] < param_dist[nearest_better_index]) {
        if (objective_dist[j] < objective_dist[nearest_better_index]) {
          nearest_better_index = j;
        }
        
        // if (param_dist[j] > param_dist[furthest_better_index]) {
        if (objective_dist[j] > objective_dist[furthest_better_index]) {
          furthest_better_index = j;
        }
      }
      
      // check the few nearest solutions with better fitness
      for (size_t j = 0; j < std::min(i, max_number_of_neighbours); j++)
      {
        
        // find the next-to nearest index
        if (j > 0)
        {
          old_nearest_better_index = nearest_better_index;
          nearest_better_index = furthest_better_index;
          
          for (size_t k = 0; k < i; k++)
          {
            //if (param_dist[k] > param_dist[old_nearest_better_index] && param_dist[k] < param_dist[nearest_better_index]) {
            if (objective_dist[k] > objective_dist[old_nearest_better_index] && objective_dist[k] < objective_dist[nearest_better_index]) {
              nearest_better_index = k;
            }
          }
          
        }
        
        if (!recheck_elites)
        {
          // if both are elites that belong to the same archive, accept edge.
          // if both are elites that belong to different archives, reject edge. If you dont re-check, the number of clusters explodes as they can (almost) never be merged again.
          
          if( pop.sols[sols_sorted_on_fitness[i]]->elite_origin != nullptr && pop.sols[sols_sorted_on_fitness[nearest_better_index]]->elite_origin != nullptr)
          {
            
            if(pop.sols[sols_sorted_on_fitness[i]]->elite_origin ==  pop.sols[sols_sorted_on_fitness[nearest_better_index]]->elite_origin)
            {
              // accept edge
              cluster_index[sols_sorted_on_fitness[i]][objective_number] = cluster_index[sols_sorted_on_fitness[nearest_better_index]][objective_number];
              edge_added = true;
              break;
            }
            else
            {
              // reject edge
              // double check them, so that we combine more clusters.
              // continue;
            }
          }
        }
        
        // set budget
        double param_dist = pop.sols[sols_sorted_on_fitness[i]]->param_distance(*pop.sols[sols_sorted_on_fitness[nearest_better_index]]);
        max_number_of_trial_solutions = 1 + ((unsigned int)(param_dist / average_edge_length));
        
        // such that mini < maxi
        mini = std::min(sols_sorted_on_fitness[nearest_better_index], sols_sorted_on_fitness[i]);
        maxi = std::max(sols_sorted_on_fitness[nearest_better_index], sols_sorted_on_fitness[i]);
        test_solution_hash = maxi * initial_popsize + mini;
        
        // there are no pre-saved test_solutions yet
        if(test_solution_hashes.find(test_solution_hash) == test_solution_hashes.end()) {
          test_solution_hashes[test_solution_hash] = {};
        }
        
        if (hill_valley_test_single_objective_reuse(*pop.sols[mini], *pop.sols[maxi], objective_number, max_number_of_trial_solutions, number_of_evaluations, test_solution_hashes[test_solution_hash]))
        {
          cluster_index[sols_sorted_on_fitness[i]][objective_number] = cluster_index[sols_sorted_on_fitness[nearest_better_index]][objective_number];
          edge_added = true;
          break;
        }
      }
      
      
      // its a new clusters, label it like that.
      if (!edge_added)
      {
        cluster_index[sols_sorted_on_fitness[i]][objective_number] = number_of_clusters[objective_number];
        number_of_clusters[objective_number]++;
      }
    }
  }
  
  // create & fill the clusters
  //---------------------------------------------------------------------------
  // map cluster indices:
  // n0 = 0,...,N0 - 1
  // n1 = 0,...,N1 - 1
  // n2 = 0,...,N2 - 1
  // current_cluster_index = N0*N1*n2 + N0*n1 + n0
  
  vec_t terms(number_of_objectives,1.0);
  
  for(size_t j = 0; j < number_of_objectives; ++j)
  {
    for(size_t subj = 0; subj < j; ++subj) {
      terms[j] *= number_of_clusters[subj];
    }
    
  }
  
  std::vector<population_pt> new_clusters;
  std::map<int, size_t> cluster_hashes;
  
  int current_hash;
  for (size_t i = 0; i < initial_popsize; ++i)
  {
    current_hash = (int) terms.dot(cluster_index[i]);
    
    // if the hash is not found yet, add it and create corresponding new cluster.
    if(cluster_hashes.find(current_hash) == cluster_hashes.end())
    {
      cluster_hashes[current_hash] = new_clusters.size();
      new_clusters.push_back(std::make_shared<population_t>());
    }
    
    new_clusters[cluster_hashes[current_hash]]->sols.push_back(pop.sols[i]);
  }
  
  // now i want to add the test solutions to the clusters.
  if( add_test_solutions )
  {
    int clusteri, clusterj;
    for(size_t i = 0; i < initial_popsize; ++i)
    {
      clusteri = terms.dot(cluster_index[i]);
      
      for(size_t j = i + 1; j < initial_popsize; ++j)
      {
        clusterj = terms.dot(cluster_index[j]);
        
        if(clusteri == clusterj) // only if the two solutions belong to the same cluster
        {
          test_solution_hash = j * initial_popsize + i;
          for(size_t k = 0; k < test_solution_hashes[test_solution_hash].size(); ++k)
          {
            if(test_solution_hashes[test_solution_hash][k]->number_of_objectives() > 0) // only if it is actually evaluated
            {
              new_clusters[cluster_hashes[clusteri]]->sols.push_back(test_solution_hashes[test_solution_hash][k]);
              new_clusters[cluster_hashes[clusteri]]->sols.back()->population_number = optimizer_number;
            }
          }
        }
      }
    }
  }
  
  
  // add only the clusters with sols
  for(size_t i = 0; i < new_clusters.size(); ++i)
  {
    subpopulations.push_back(new_clusters[i]);
    int cluster_number = (int) (subpopulations.size() - 1);
    
    for(size_t j = 0; j < subpopulations[cluster_number]->sols.size(); ++j) {
      subpopulations.back()->sols[j]->cluster_number = cluster_number;
    }
  }
  
  
}
