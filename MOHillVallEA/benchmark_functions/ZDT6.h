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

  class ZDT6_t : public fitness_t
  {

  public:

    ZDT6_t()
    {
      number_of_objectives = 2; // fixed
      number_of_parameters = 2; // default, can be adapted
    }
    ~ZDT6_t() {}

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
      sol.obj[0] = 1.0 - exp(-4.0*sol.param[0])*pow(sin(6.0*3.1415926536*sol.param[0]), 6.0);
    
      // f2
      double g = 0.0;
      for (size_t i = 1; i < number_of_parameters; i++) {
        g += sol.param[i] / (number_of_parameters - 1.0);
      }
      g = 1.0 + (number_of_parameters - 1.0)*pow(g, 0.25);

      double h = 1.0 - (sol.obj[0] / g)*(sol.obj[0] / g);

      sol.obj[1] = g*h;
      sol.penalty = 0;
    }

    std::string name() const
    {
      return "ZDT6";
    }

    // compute VTR in terms of the D_{\mathcal{P}_F}\rightarrow\mathcal{S}
    bool get_pareto_set()
    {
      size_t pareto_set_size = 5000;
      
      // generate default front
      if (pareto_set.size() != pareto_set_size)
      {
        
        pareto_set.sols.clear();
        pareto_set.sols.reserve(pareto_set_size);
        
        // the front
        for (size_t i = 0; i < pareto_set_size; ++i)
        {
          solution_pt sol = std::make_shared<solution_t>(number_of_parameters, number_of_objectives);
          
          sol->param.fill(0.0);
          sol->param[0] = (i / ((double)pareto_set_size - 1.0));
          
          define_problem_evaluation(*sol); // runs a feval without registering it.
          
          pareto_set.sols.push_back(sol);
        }
        
        // pareto_set.writeToFile("./ZDT6");
        igdx_available = true;
        igd_available = true;
      }
      
      return true;
    }


  };
}
