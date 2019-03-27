#pragma once

/*
 
 MO-HillVallEA

 By S.C. Maree, 2019
 s.c.maree[at]amc.uva.nl
 github.com/scmaree
 
 */

#include "../core/fitness.h"
#include "../core/elitist_archive.h"

namespace hillvallea
{

  class ZDT3_t : public fitness_t
  {

  public:

    ZDT3_t()
    {
      number_of_objectives = 2; // fixed
      number_of_parameters = 2; // default, can be adapted
    }
    ~ZDT3_t() {}

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
      upper.resize(number_of_parameters, 1.0);

    }

    void define_problem_evaluation(solution_t & sol)
    {
    
      // f1
      sol.obj[0] = sol.param[0];
    
      // f2
      double g = 0.0;
      for (size_t i = 1; i < number_of_parameters; i++) {
        g += sol.param[i] / (number_of_parameters - 1.0);
      }
      g = 1.0 + 9.0*g;

      double h = 1.0 - sqrt(sol.obj[0] / g) - (sol.obj[0] / g)*sin(10.0*3.1415926536*sol.obj[0]);

      sol.obj[1] = g*h;
      sol.penalty = 0;
    }

    std::string name() const
    {
      return "ZDT3";
    }

    // compute VTR in terms of the D_{\mathcal{P}_F}\rightarrow\mathcal{S}
    bool get_pareto_set()
    {
      
      size_t pareto_set_size = 5000;
      
      // generate default front
      // this one is annyoing.
      // the front is only part of the described curve.
      // so we create an archive out of it and discard the dominated solutions.
      if (pareto_set.size() < 10)
      {
        rng_pt rng = std::make_shared<rng_t>(100); // not used anyways as the archive size is never adapted here
        elitist_archive_t temp_archive(5000, rng);
        
        size_t temp_pareto_set_size = 18975;
        
        pareto_set.sols.clear();
        pareto_set.sols.reserve(pareto_set_size);
        
        // the front
        for (size_t i = 0; i < temp_pareto_set_size; ++i)
        {
          solution_pt sol = std::make_shared<solution_t>(number_of_parameters, number_of_objectives);
          
          sol->param.fill(0.0);
          sol->param[0] = (i / ((double)temp_pareto_set_size - 1.0));
          
          define_problem_evaluation(*sol); // runs a feval without registering it.
          
          temp_archive.updateArchive(sol);
          
        }
        
        pareto_set.sols = temp_archive.sols;
        
        // pareto_set.writeToFile("./ZDT3.txt");
        igdx_available = true;
        igd_available = true;
      }
      
      return true;
    }


  };
}
