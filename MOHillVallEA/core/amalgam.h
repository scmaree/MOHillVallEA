#pragma once

/*
 
 AMaLGaM - implementation as part of MO-HillVallEA
 
 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "hillvallea_internal.h"
#include "cluster.h"

namespace hillvallea
{

  class amalgam_t : public cluster_t
  {

  public: 

    // C++ Rule of Three
    //-------------------------------------------
    amalgam_t(fitness_pt fitness_function, rng_pt rng, bool use_univariate);
    amalgam_t(const amalgam_t & other);
    ~amalgam_t();

    // Run-time control
    //---------------------------------------------------------------------------------
    bool checkTerminationCondition();
    void estimateParameters();
    std::string name() const;

    bool updateStrategyParameters(const elitist_archive_t & elitist_archive, size_t & no_improvement_stretch, const size_t maximum_no_improvement_stretch);
    void computeParametersForSampling();
    void generateNewSolutions(std::vector<solution_pt> & solutions, size_t number_of_solutions, size_t number_of_ams_solutions, rng_pt & rng);
    size_t recommended_popsize(size_t number_of_parameters) const;

    vec_t getVec(std::string variable_name);
    double getDouble(std::string variable_name);


  private:


    // Essential data members
    //-------------------------------------------
    vec_t mean;                       // sample mean
    double multiplier;                // distribution multiplier
    
    // full matrices
    matrix_t covariance;              // sample covariance matrix C
    matrix_t cholesky;                // decomposed covariance matrix C = LL^T
    matrix_t inverse_cholesky;        // inverse of the cholesky decomposition

    // univariate 'matrices'
    vec_t univariate_covariance;
    vec_t univariate_cholesky;
    vec_t univariate_inverse_cholesky;

    bool use_boundary_repair;
    bool use_univariate;

    // Stopping criteria
    //-------------------------------------------
    bool terminated;
    bool improvement;
    double out_of_bounds_sample_ratio;

    //// AMS & SDR
    ////-------------------------------------------
    void apply_ams(std::vector<solution_pt> & solutions, const size_t number_of_elites, const size_t number_of_ams_solutions, const double ams_factor, const vec_t & ams_direction, const vec_t & lower_param_bounds, const vec_t & upper_param_bounds) const;
    double getSDR(const vec_t & params) const;
    
    
    bool adaptDistributionMultiplier(const elitist_archive_t & elitist_archive, size_t & no_improvement_stretch, const size_t maximum_no_improvement_stretch);
    bool generationalImprovementForOneCluster(double & st_dev_ratio, const elitist_archive_t & elitist_archive) const;
    };

}
