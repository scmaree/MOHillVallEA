#pragma once

/*
 
 MO-HillVallEA

 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "../core/fitness.h"

namespace hillvallea
{

  class OmniTest_t : public fitness_t
  {

  public:

    OmniTest_t()
    {
      number_of_objectives = 2; // fixed
      number_of_parameters = 5; // default, can be adapted
    }
    ~OmniTest_t() {}

    // number of objectives 
    void set_number_of_objectives(size_t & number_of_objectives)
    {
      this->number_of_objectives = 2;
      number_of_objectives = this->number_of_objectives;
    }

    // any positive value
    void set_number_of_parameters(size_t & number_of_parameters)
    {
      this->number_of_parameters = number_of_parameters;
    }


    void get_param_bounds(vec_t & lower, vec_t & upper) const
    {

      lower.clear();
      lower.resize(number_of_parameters, 0.0);
      
      upper.clear();
      upper.resize(number_of_parameters, 6.0);
    }

    void define_problem_evaluation(solution_t & sol)
    {
    
      // f1
      sol.obj[0] = 0.0;
      
      for (size_t i = 0; i < number_of_parameters; ++i) {
        sol.obj[0] += sin(PI*sol.param[i]);
      }
      
      // f2
      sol.obj[1] = 0.0;
      for (size_t i = 0; i < number_of_parameters; ++i) {
        sol.obj[1] += cos(PI*sol.param[i]);
      }
      
      
      // penalty
      sol.penalty = 0.0;
    }

    std::string name() const
    {
      return "OmniTest";
    }

    // compute VTR in terms of the D_{\mathcal{P}_F}\rightarrow\mathcal{S}
    bool get_pareto_set()
    {
      
      size_t pareto_set_size = 10000;
      
      // generate default front
      if (pareto_set.size() != pareto_set_size)
      {
        
        pareto_set.sols.clear();
        pareto_set.sols.reserve(pareto_set_size);
        
        rng_pt rng = std::make_shared<rng_t>(100);
        std::uniform_real_distribution<double> unif(0, 1);
        
        int tilePosition;
        double max_var;
        double min_var;
        double randomVar;
        // the front
        for (size_t i = 0; i < pareto_set_size; ++i)
        {
          solution_pt sol = std::make_shared<solution_t>(number_of_parameters, number_of_objectives);
          sol->param.fill(0.0);
          
          randomVar = unif(*rng);
          
          for(size_t j = 0; j < sol->param.size(); ++j)
          {
            tilePosition = (int) (3*unif(*rng));
            max_var = 2.0*tilePosition + 1.5;
            min_var = 2.0*tilePosition + 1.0;

            sol->param[j] = min_var + (max_var - min_var) * randomVar;
          }
          
          define_problem_evaluation(*sol); // runs a feval without registering it.
          
          pareto_set.sols.push_back(sol);
        }
        
        igdx_available = true;
        igd_available = true;
        
      }
      
      return true;
    }

  };
}
