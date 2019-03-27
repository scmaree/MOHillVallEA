#define _CRT_SECURE_NO_WARNINGS

/*
 
 MO-HillVallEA
 
 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "recursion_scheme.h"
#include "optimizer.h"
#include "elitist_archive.h"
#include "population.h"
#include "solution.h"
#include "mathfunctions.h"
#include "fitness.h"
#include "hillvalleyclustering.h"

hillvallea::recursion_scheme_t::recursion_scheme_t
(
  fitness_pt fitness_function,
  const vec_t & lower_init_ranges,
  const vec_t & upper_init_ranges,
  double vtr,
  bool use_vtr,
  int version,
  int local_optimizer_index,
  size_t elitist_archive_size_target,
  size_t approximation_set_size,
  size_t maximum_number_of_populations,
  int base_population_size,
  unsigned int number_of_subgenerations_per_population_factor,
  unsigned int maximum_number_of_evaluations,
  unsigned int maximum_number_of_seconds,
  int random_seed,
  bool write_generational_solutions,
  bool write_generational_statistics,
  bool print_generational_statistics,
  const std::string & write_directory,
  const std::string & file_appendix,
  bool print_verbose_overview
)
{
  this->fitness_function = fitness_function;
  this->lower_init_ranges = lower_init_ranges;
  this->upper_init_ranges = upper_init_ranges;
  this->maximum_number_of_populations = maximum_number_of_populations;
  this->base_population_size = base_population_size;
  base_number_of_mixing_components = -1;
  this->maximum_number_of_evaluations = maximum_number_of_evaluations;
  this->maximum_number_of_seconds = maximum_number_of_seconds;
  this->vtr = vtr;
  this->use_vtr = use_vtr;
  this->version = version;
  this->local_optimizer_index = local_optimizer_index;
  this->random_seed = random_seed;
  this->write_generational_solutions = write_generational_solutions;
  this->write_generational_statistics = write_generational_statistics;
  this->print_generational_statistics = print_generational_statistics;
  this->write_directory = write_directory;
  this->file_appendix = file_appendix;
  this->number_of_subgenerations_per_population_factor = number_of_subgenerations_per_population_factor;
  this->print_verbose_overview = print_verbose_overview;
  
  this->elitist_archive_size_target = elitist_archive_size_target;
  this->approximation_set_size = approximation_set_size;

  // truncate initialization ranges to allowable ranges in of the objective function
  vec_t lower_range_bounds, upper_range_bounds;
  fitness_function->get_param_bounds(lower_range_bounds, upper_range_bounds);

  for (size_t i = 0; i < fitness_function->get_number_of_parameters(); ++i)
  {
    this->lower_init_ranges[i] = std::max(lower_init_ranges[i], lower_range_bounds[i]);
    this->upper_init_ranges[i] = std::min(upper_init_ranges[i], upper_range_bounds[i]);
  }

}

hillvallea::recursion_scheme_t::~recursion_scheme_t() {}



bool hillvallea::recursion_scheme_t::checkNumberOfEvaluationsTerminationCondition() const
{
  // or use : fitness_function->number_of_evaluations
  // could be used as a check as well

  assert(fitness_function->number_of_evaluations == number_of_evaluations);
  //todo?: update from fevals_start.
  if (maximum_number_of_evaluations > 0 && number_of_evaluations >= maximum_number_of_evaluations) {
    return true;
  }

  return false;
}

/*
 Returns 1 if the value-to-reach has been reached for the multi-objective case. This means that
 the D_Pf->S metric has reached the value-to-reach. If no D_Pf->S can be computed, 0 is returned.
*/
// non-const because value_to_reach_reached can store pre-computed stuff. 
bool hillvallea::recursion_scheme_t::checkVTRTerminationCondition(const elitist_archive_t & approximation_set)
{

  if (use_vtr)
  {
    // return (approximation_set.computeIGD(fitness_function->pareto_set) <= vtr);
    double threshold = 0;
    
    if(fitness_function->number_of_objectives == 2) {
      threshold = 5e-2;
    }
    else {
      threshold = 1e-1;
    }
    
    vec_t ones(fitness_function->pareto_sets.size(),1.0);
    return (approximation_set.computeSR(fitness_function->pareto_sets, threshold, ones) == 1.0);
  }

  return false;

}

double hillvallea::recursion_scheme_t::getTimer() const 
{
  clock_t end_time = clock();
  return (((double)(end_time - start_time)) / CLOCKS_PER_SEC);
}

bool hillvallea::recursion_scheme_t::checkTimeLimitTerminationCondition() 
{

  if (maximum_number_of_seconds > 0 && getTimer() > maximum_number_of_seconds) {
    return true;
  }

  return false;
}

bool hillvallea::recursion_scheme_t::checkGlobalTerminationCondition(const elitist_archive_t & approximation_set)
{
  if (populations.size() == 0) {
    return false;
  }

  if (checkNumberOfEvaluationsTerminationCondition()) {
    // std::cout << "Terminated because of evaluations" << std::endl;
    return true;
  }

  if (checkTimeLimitTerminationCondition()) {
    // std::cout << "Terminated because of time" << std::endl;
    return true;
  }

  return false;
}


bool hillvallea::recursion_scheme_t::checkTerminationConditionAllPopulations(const elitist_archive_t & approximation_set)
{

  if (checkGlobalTerminationCondition(approximation_set)) {
    return true;
  }

  if (checkVTRTerminationCondition(approximation_set)) {
    // std::cout << "Terminated because of VTR" << std::endl;
    return true;
  }

  for (size_t i = 0; i < populations.size(); ++i)
  {
    if (populations[i]->checkTerminationCondition()) {
      populations[i]->terminated = true;
    }
  }

  if (populations.size() < maximum_number_of_populations) {
    return false;
  }

  bool all_populations_terminated = true;
  for (size_t i = 0; i < populations.size(); ++i)
  {
    if (!populations[i]->terminated) {
      all_populations_terminated = false;
      break;
    }
  }

  return all_populations_terminated;
}

bool hillvallea::recursion_scheme_t::checkTerminationConditionOnePopulation(optimizer_t & population)
{

  if (checkGlobalTerminationCondition(*approximation_set)) {
    return true;
  }

  if (population.checkTerminationCondition()) {
    population.terminated = true;
  }

  return false;
}


void hillvallea::recursion_scheme_t::generationalStepAllPopulationsRecursiveFold(size_t population_index_smallest, size_t population_index_biggest)
{
 
  for (unsigned int i = 0; i < number_of_subgenerations_per_population_factor - 1; i++)
  {
    for (size_t population_index = population_index_smallest; population_index <= population_index_biggest; population_index++)
    {
      if (!populations[population_index]->terminated)
      {

        switch (version)
        {

        case 0: 
          populations[population_index]->generation(*elitist_archive, number_of_evaluations);
          break;

        default:
          bool largest_population = (population_index == populations.size() - 1);
          populations[population_index]->generation_mm(*elitist_archive, number_of_evaluations, largest_population);
          break;
        }

        // if any of the global termininations is hit (fevals, time, vtr), stop. 
        if (checkTerminationConditionOnePopulation(*populations[population_index]))
        {
          for (size_t j = 0; j < populations.size(); j++) {
            populations[j]->terminated = true;
          }
          return;
        }
        
        // if(maximum_number_of_populations == 1)
        {
          // if(populations[population_index]->number_of_generations % 10 == 0) // reduce computation time
          {
            if(write_generational_solutions || write_generational_statistics)
            {
              approximation_set->computeApproximationSet(approximation_set_size, populations, elitist_archive, (version != 0), false);
            }
            
            if (write_generational_statistics) {
              writeGenerationalStatisticsForOnePopulation(*populations[population_index], population_index, use_vtr, reached_vtr, print_generational_statistics,false);
            }
            if (write_generational_solutions)
            {
              for(size_t i = 0; i < populations.size(); ++i)
              {
                if(!populations[i]->terminated) {
                  writePopulation(i, *populations[i]);
                }
              }
              
              writeGenerationalSolutions(*approximation_set, false, number_of_solution_sets_written);
              
            }
          }
        }
      }
    }

    for (size_t population_index = population_index_smallest; population_index < population_index_biggest; population_index++) {
      generationalStepAllPopulationsRecursiveFold(population_index_smallest, population_index);
    }
  }
}

void hillvallea::recursion_scheme_t::generationalStepAllPopulations()
{

  size_t population_index_biggest = populations.size() - 1;
  size_t population_index_smallest = 0;

  while (population_index_smallest <= population_index_biggest)
  {
    if (!populations[population_index_smallest]->terminated) {
      break;
    }

    population_index_smallest++;
  }

  generationalStepAllPopulationsRecursiveFold(population_index_smallest, population_index_biggest);
  
  approximation_set->computeApproximationSet(approximation_set_size, populations, elitist_archive, (version != 0), true);
}


void hillvallea::recursion_scheme_t::run()
{

  initialize();

  if (print_verbose_overview) {
    printVerboseOverview();
  }

  while (!checkTerminationConditionAllPopulations(*approximation_set))
  {

    size_t population_size, number_of_mixing_components;
    if (populations.size() < maximum_number_of_populations)
    {

      if (populations.size() == 0)
      {
        population_size = base_population_size;
        number_of_mixing_components = base_number_of_mixing_components;
      }
      else
      {
        population_size = 2 * populations.back()->population_size;
        number_of_mixing_components = (int) (populations.back()->number_of_mixing_components * 1.5); // fitness_function->get_number_of_objectives();
      }

      int population_index = (int) populations.size();
      
      optimizer_pt population = std::make_shared<optimizer_t>(
        fitness_function, lower_init_ranges, upper_init_ranges,
        local_optimizer_index, population_size, number_of_mixing_components, population_index,
        maximum_number_of_evaluations, maximum_number_of_seconds, 
        vtr, use_vtr, rng);

      switch (version)
      {
        case 0:
          population->initialize(*elitist_archive, number_of_evaluations);
          break;
          
        default:
          population->initialize_mm(*elitist_archive, number_of_evaluations);
          break;
      }

      populations.push_back(population);

      if(write_generational_solutions || write_generational_statistics || maximum_number_of_populations > 0)
      {
        approximation_set->computeApproximationSet(approximation_set_size, populations, elitist_archive, (version != 0), false);
      }
      
      if (write_generational_statistics) {
        writeGenerationalStatisticsForOnePopulation(*populations.back(), populations.size()-1, use_vtr, reached_vtr, print_generational_statistics,( total_number_of_generations == 0 && populations.size() == 1));
      }
      if (write_generational_solutions)
      {
        for(size_t i = 0; i < populations.size(); ++i)
        {
          if(!populations[i]->terminated) {
            writePopulation(i, *populations[i]);
          }
        }
        
        writeGenerationalSolutions(*approximation_set, false, number_of_solution_sets_written);
        
      }
      
    }

    generationalStepAllPopulations();

    total_number_of_generations++;
  }

  // finish up
  approximation_set->computeApproximationSet(approximation_set_size, populations, elitist_archive, (version != 0), true);
  
  
  writeGenerationalStatisticsForOnePopulation(*populations[populations.size() - 1], populations.size() - 1, use_vtr, reached_vtr, print_generational_statistics, ( total_number_of_generations == 0));
  if (write_generational_solutions)
  {
    
    for(size_t i = 0; i < populations.size(); ++i)
    {
      if(!populations[i]->terminated) {
        writePopulation(i, *populations[i]);
      }
    }
    
    // write also the pareto set
    char  string[1000];
    sprintf(string, "%spareto_set%s.txt", write_directory.c_str(), file_appendix.c_str());

    fitness_function->pareto_set.writeToFile(string);
    
    writeGenerationalSolutions(*approximation_set, true, number_of_solution_sets_written);
  }
}

void hillvallea::recursion_scheme_t::initialize()
{

  // i don't really know why its defined like this, but there's an infinite loop when number_of_subgenerations_per_population_factor == 1
  assert(number_of_subgenerations_per_population_factor > 1);
  
  number_of_evaluations = 0;
  total_number_of_generations = 0;
  populations.clear();
  start_time = clock();
  number_of_solution_sets_written = 0;
  
  rng = std::make_shared<std::mt19937>((unsigned long)(random_seed));

  elitist_archive = std::make_shared<elitist_archive_t>(elitist_archive_size_target, rng);
  approximation_set = std::make_shared<elitist_archive_t>(approximation_set_size, rng);

  if(version != 0) {
    elitist_archive->use_niching = true;
  }
  else
  {
    elitist_archive->use_niching = false;
  }
  
  elitist_archive->set_use_parameter_distances(false);
  elitist_archive->set_use_greedy_selection(false);
  //}
  

  // try to generate the pareto set (and check if the fitness function has the VTR available)
  if(!fitness_function->get_pareto_set()) {
    use_vtr = false;
  }
  
  // if the SR is not yet activated, try to activate it
  if(!fitness_function->sr_available)
  {
    
    // but for that, we need the IGDX
    if(fitness_function->igdx_available && fitness_function->igd_available)
    {
      hvc_t hvc(fitness_function);
      double average_edge_length = 0.01; // irrelevant?
      rng_pt rng2 = std::make_shared<rng_t>(2001);
      unsigned int fevals = 0;
      
      // find pareto sets from the Pareto set.
      fitness_function->pareto_sets.clear();
      hvc.cluster(fitness_function->pareto_set, fitness_function->pareto_sets, fevals, average_edge_length,false, false, 0, rng2);
      fitness_function->number_of_evaluations = 0;
      fitness_function->pareto_sets_max_igdx.reset(fitness_function->pareto_sets.size(), 0.0);
      
      size_t local_approximation_set_size = std::max((size_t) 1, approximation_set_size / fitness_function->pareto_sets.size());
      population_t collected_sols;
      for(size_t i = 0 ; i < fitness_function->pareto_sets.size(); ++i)
      {
        population_t skinny_pareto_set, non_selected_sols;
        selectSolutionsBasedOnParameterDiversity(fitness_function->pareto_sets[i]->sols, local_approximation_set_size, skinny_pareto_set.sols, non_selected_sols.sols, rng2);
        fitness_function->pareto_sets_max_igdx[i] = skinny_pareto_set.computeIGDX(*fitness_function->pareto_sets[i]);
        collected_sols.addSolutions(skinny_pareto_set);
      }
      
      // std::cout << collected_sols.computeIGDX(fitness_function->pareto_set);
      // std::cout << fitness_function->pareto_sets_max_igdx;
      fitness_function->sr_available = true;
    }
  } // end computing SR
  
  if (maximum_number_of_populations == 1)
  {
    if (base_number_of_mixing_components <= 0) {
      base_number_of_mixing_components = 20;
    }

    if(base_population_size < 0)
    {
      size_t cluster_size = 0;
      
      hillvallea::rng_pt rng = std::make_shared<hillvallea::rng_t>(100); // idk.
      hillvallea::cluster_pt dummy = hillvallea::init_cluster(local_optimizer_index, fitness_function, rng);
      
      cluster_size = dummy->recommended_popsize(fitness_function->number_of_parameters);
      
      base_population_size = (int)(0.5*base_number_of_mixing_components*cluster_size);
    }
  }
  else
  {
    
    if(base_population_size < 0)
    {
      base_number_of_mixing_components = 1 + (int) fitness_function->get_number_of_objectives();
      base_population_size = 10 * base_number_of_mixing_components  * (int) (1 + log( (double) fitness_function->get_number_of_parameters()));
    }
    else
    {
      base_number_of_mixing_components = base_population_size / 10;
      // base_population_size = 10 * base_number_of_mixing_components;
    }
  }


}

void hillvallea::recursion_scheme_t::printVerboseOverview() const
{

  cluster_pt dummy_cluster = init_cluster(local_optimizer_index, fitness_function, rng);

  printf("### Settings ######################################\n");
  printf("#\n");
  printf("# Statistics writing every generation: %s\n", write_generational_statistics ? "enabled" : "disabled");
  printf("# Population file writing            : %s\n", write_generational_solutions ? "enabled" : "disabled");
  std::cout << "# Write Directory                    : " << write_directory << std::endl;
  std::cout << "# File Appendix                      : " << file_appendix << std::endl;
  printf("# Use of value-to-reach (vtr)        : %s\n", use_vtr ? "enabled" : "disabled");
  std::cout << "# Version                            : " << version << std::endl;
  std::cout << "# Local optimizer                    : " << dummy_cluster->name() << std::endl;
  printf("#\n");
  printf("###################################################\n");
  printf("#\n");
  printf("# Problem                  = %s\n", fitness_function->name().c_str());
  printf("# Number of parameters     = %zu\n", fitness_function->get_number_of_parameters());
  printf("# Number of objectives     = %zu\n", fitness_function->get_number_of_objectives());
  
  rng_pt rng = std::make_shared<rng_t>(110);
  population_t skinny_pareto_set, non_selected_sols;
  
  if(fitness_function->igd_available)
  {
    vec_t objective_ranges;
    fitness_function->pareto_set.objectiveRanges(objective_ranges);
    selectSolutionsBasedOnObjectiveDiversity(fitness_function->pareto_set.sols, approximation_set_size, skinny_pareto_set.sols, objective_ranges, non_selected_sols.sols, rng);
    printf("# Max IGD                  = %e\n", skinny_pareto_set.computeIGD(fitness_function->pareto_set));
  }
  
  if(fitness_function->igdx_available)
  {
    skinny_pareto_set.sols.clear();
    selectSolutionsBasedOnParameterDiversity(fitness_function->pareto_set.sols, approximation_set_size,skinny_pareto_set.sols, non_selected_sols.sols, rng);
    printf("# Max IGDX                 = %e\n", skinny_pareto_set.computeIGDX(fitness_function->pareto_set));
  }
  
  if(fitness_function->sr_available)
  {
    printf("# Number of Pareto Sets    = %zu\n", fitness_function->pareto_sets.size());
  }
  
  printf("# Initialization ranges    = ");
  // for (size_t i = 0; i < fitness_function->get_number_of_parameters(); i++)
  {
    size_t i = fitness_function->get_number_of_parameters() - 1;
    printf("x_%zu: [%5.2f; %5.2f]", i, lower_init_ranges[i], upper_init_ranges[i]);
    if (i < fitness_function->get_number_of_parameters() - 1)
      printf("\n#                            ");
  }
  printf("\n");
  printf("# Boundary ranges          = ");
  // for (size_t i = 0; i < fitness_function->get_number_of_parameters(); i++)
  {
    size_t i = fitness_function->get_number_of_parameters() - 1;
    vec_t lower_range_bounds, upper_range_bounds;
    fitness_function->get_param_bounds(lower_range_bounds, upper_range_bounds);

    printf("x_%zu: [%5.2f;%5.2f]", i, lower_range_bounds[i], upper_range_bounds[i]);
    if (i < fitness_function->get_number_of_parameters() - 1) {
      printf("\n#                            ");
    }
  }
  
  
  printf("\n");
  // printf("# Rotation angle           = %e\n", rho);
  // printf("# Tau                      = %e\n", tau);
  printf("# Population size          = %i\n", base_population_size);
  printf("# Number of populations    = %zu\n", maximum_number_of_populations);
  printf("# Number of mix. com. (k)  = %i\n", base_number_of_mixing_components);
  printf("# IMS Subgenerations       = %u\n", number_of_subgenerations_per_population_factor);
  // printf("# Dis. mult. decreaser     = %e\n", distribution_multiplier_decrease);
  // printf("# St. dev. rat. threshold  = %e\n", st_dev_ratio_threshold);
  printf("# Elitist ar. size target  = %zu\n", elitist_archive_size_target);
  printf("# Approximation set size   = %zu\n", approximation_set_size);
  printf("# Maximum numb. of eval.   = %d\n", maximum_number_of_evaluations);
  printf("# Value to reach (vtr)     = %e\n", vtr);
  printf("# Time limit (s)           = %d\n", maximum_number_of_seconds);
  printf("# Random seed              = %ld\n", (long)random_seed);
  printf("#\n");
  printf("###################################################\n");
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


void hillvallea::recursion_scheme_t::writeGenerationalStatisticsForOnePopulation(const optimizer_t & population, size_t population_number, bool & use_vtr, double & reached_vtr, bool print_output, bool write_header) const
{

  char    string[1000];
  bool   enable_hyper_volume;
  FILE   *file;

  enable_hyper_volume = (fitness_function->get_number_of_objectives() == 2);
  file = NULL;
  if (write_header)
  {

    sprintf(string, "%sstatistics%s.txt", write_directory.c_str(), file_appendix.c_str());
    file = fopen(string, "w");

    sprintf(string, "# Generation  Evaluations  Time (s)");
    fputs(string, file); if (print_output) std::cout << string;
    for (size_t i = 0; i < fitness_function->get_number_of_objectives(); i++)
    {
      sprintf(string, "  Best_obj[%zu]", i);
      fputs(string, file); if (print_output) std::cout << string;
    }
    //if (enable_hyper_volume)
    //{
      sprintf(string, "  Hypervol. ");
      fputs(string, file); if (print_output) std::cout << string;
    //}
    
    if(fitness_function->igd_available)
    {
      sprintf(string, "  IGD");
      fputs(string, file); if (print_output) std::cout << string;
    }
    
    if(fitness_function->igdx_available)
    {
      sprintf(string, "         IGDX");
      fputs(string, file); if (print_output) std::cout << string;
    }
    
    if(fitness_function->sr_available)
    {
      sprintf(string, "        SR");
      fputs(string, file); if (print_output) std::cout << string;
    }
    
    sprintf(string, "    Appr.set");
    fputs(string, file); if (print_output) std::cout << string;
    
    //sprintf(string, "    El.A size");
    //fputs(string, file); if (print_output) std::cout << string;
    
    sprintf(string, "  [ Pop.index  Subgen.  Pop.size  Subpops ]\n");
    fputs(string, file); if (print_output) std::cout << string;
    
  }
  else {
    sprintf(string, "%sstatistics%s.txt", write_directory.c_str(), file_appendix.c_str());
    file = fopen(string, "a");
  }

  sprintf(string, "  %10d %11d %11.3f", total_number_of_generations, number_of_evaluations, getTimer());
  fputs(string, file); if (print_output) std::cout << string;

  for (size_t i = 0; i < fitness_function->get_number_of_objectives(); i++)
  {
    sprintf(string, " %11.3e", approximation_set->get_best_objective_values_in_elitist_archive(i));
    fputs(string, file); if (print_output) std::cout << string;
  }

  if (enable_hyper_volume)
  {
    sprintf(string, " %11.3e", approximation_set->compute2DHyperVolume());
    fputs(string, file); if (print_output) std::cout << string;
  }
  else
  {
    sprintf(string, " %11.3e", 0.0);
    fputs(string, file); if (print_output) std::cout << string;
  }
  
  if (fitness_function->igd_available)
  {
    sprintf(string, " %11.3e", approximation_set->computeIGD(fitness_function->pareto_set));
    fputs(string, file); if (print_output) std::cout << string;
  }

  if (fitness_function->igdx_available)
  {
    sprintf(string, " %11.3e", approximation_set->computeIGDX(fitness_function->pareto_set));
    fputs(string, file); if (print_output) std::cout << string;
  }
  
  if (fitness_function->sr_available)
  {
    vec_t ones(fitness_function->pareto_sets_max_igdx.size(), 1.0);
    
    double threshold;
    if(fitness_function->number_of_objectives == 2) {
      threshold = 5e-2;
    }
    else {
      threshold = 1e-1;
    }
    
    sprintf(string, " %7.4f", approximation_set->computeSR(fitness_function->pareto_sets, threshold, ones));
    fputs(string, file); if (print_output) std::cout << string;
  }

  sprintf(string, " %11zu", approximation_set->size());
  fputs(string, file); if (print_output) std::cout << string;
  
  sprintf(string, "  [ %4zu %6d %10zu %10zu ]\n", population_number, population.number_of_generations, population.population_size, population.subpopulations.size());
  fputs(string, file); if (print_output) std::cout << string;

  fclose(file);
}

void hillvallea::recursion_scheme_t::writeGenerationalSolutions(const elitist_archive_t & approximation_set, bool final, unsigned int & number_of_solution_sets_written) const
{
  char  string[1000];

  // Approximation set
  if (final) {
    sprintf(string, "%sapproximation_set_generation%s_final.txt", write_directory.c_str(), file_appendix.c_str());
  }
  else {
    sprintf(string, "%sapproximation_set_generation%s_%05d.txt", write_directory.c_str(), file_appendix.c_str(), number_of_solution_sets_written);
  }

  approximation_set.writeToFile(string);

  number_of_solution_sets_written++;
}

void hillvallea::recursion_scheme_t::writePopulation(size_t population_index, const optimizer_t & optimizer) const
{
  char  string[1000];
  
  sprintf(string, "%spopulation_%05zu_generation_%05d%s.txt", write_directory.c_str(), population_index, optimizer.number_of_generations-1, file_appendix.c_str());
  
  population_t population;
  population.collectSolutions(optimizer.subpopulations);
  
  population.writeToFile(string);
  
}


