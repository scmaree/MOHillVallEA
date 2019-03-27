

/*
 
 MO-HillVallEA
 
 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "cluster.h"
#include "solution.h"
#include "amalgam.h"
#include "iamalgam.h"

hillvallea::cluster_pt hillvallea::init_cluster(const int cluster_index, fitness_pt fitness_function, rng_pt rng)
{

  // parse settings
  switch (cluster_index)
  {
    case 0: return std::make_shared<amalgam_t>(fitness_function, rng, false); break; // use_univariate = false
    case 1: return std::make_shared<amalgam_t>(fitness_function, rng, true); break; // use_univariate = true
    case 10: return std::make_shared<iamalgam_t>(fitness_function, rng, false); break; // use_univariate = false
    case 11: return std::make_shared<iamalgam_t>(fitness_function, rng, true); break; // use_univariate = true
    default: return std::make_shared<amalgam_t>(fitness_function, rng, false); break;
  }

}

// initialization of the general parameters of the EDAs
hillvallea::cluster_t::cluster_t(fitness_pt fitness_function, rng_pt rng) : population_t()
{
  this->fitness_function = fitness_function;
  this->rng = rng;

  this->previous = nullptr;
  this->init_bandwidth = 0.01;

}

hillvallea::cluster_t::cluster_t(const cluster_t & other) : population_t(other)
{
  this->fitness_function = other.fitness_function;
  this->rng = other.rng;
  this->previous = other.previous;
  this->objective_mean = other.objective_mean;
  this->parameter_mean = other.parameter_mean;
  this->init_bandwidth = other.init_bandwidth;
}

hillvallea::cluster_t::~cluster_t() {}

bool hillvallea::cluster_t::checkTerminationCondition()
{
  std::cout << "cluster error, 'checkTerminiationCriteria' not implemented" << std::endl;
  assert(false);
  return true;
}

void hillvallea::cluster_t::estimateParameters()
{
  std::cout << "cluster_t error 'estimateParameters' not implemented" << std::endl;
  assert(false);

}


bool hillvallea::cluster_t::updateStrategyParameters(const elitist_archive_t & elitist_archive, size_t & no_improvement_stretch, const size_t maximum_no_improvement_stretch)
{
  std::cout << "cluster_t error 'updateStrategyParameters' not implemented" << std::endl;
  assert(false);

  return false;
}


void hillvallea::cluster_t::computeParametersForSampling()
{
  std::cout << "cluster_t error 'computeParametersForSampling' not implemented" << std::endl;
  assert(false);
}



void hillvallea::cluster_t::generateNewSolutions(std::vector<solution_pt> & solutions, size_t number_of_solutions, size_t number_of_ams_solutions, rng_pt & rng)
{
  std::cout << "cluster_t error 'generateNewSolutions' not implemented" << std::endl;
  assert(false);
}



hillvallea::vec_t hillvallea::cluster_t::getVec(std::string variable_name)
{
  std::cout << "cluster_t error 'getVec' not implemented" << std::endl;
  assert(false);

  return vec_t();
}

double hillvallea::cluster_t::getDouble(std::string variable_name)
{
  std::cout << "cluster_t error 'getDouble' not implemented" << std::endl;
  assert(false);

  return 0.0;
}

std::string hillvallea::cluster_t::name() const
{
  std::cout << "cluster_t error 'name' not implemented" << std::endl;
  assert(false);

  return "<no name specified>";
}


hillvallea::matrix_t hillvallea::cluster_t::getMatrix(std::string variable_name)
{
  std::cout << "cluster_t error 'variable_name' not implemented" << std::endl;
  assert(false);
  
  return matrix_t();
}

size_t hillvallea::cluster_t::recommended_popsize(size_t number_of_parameters) const
{
  std::cout << "cluster_t error 'recommended_popsize' not implemented" << std::endl;
  assert(false);
  
  return 0;
}
