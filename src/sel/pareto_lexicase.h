/* FEWTWO
copyright 2017 William La Cava
license: GNU/GPL v3
*/

#ifndef PARETOLEXICASE_H
#define PARETOLEXICASE_H

#include "selection_operator.h"


namespace FT{
namespace Sel{
    enum class eDominance : int { Left = -1,
                                  None  = 0,
                                  Right = 1};
                                    
    ////////////////////////////////////////////////////////////// Declarations
    /*!
     * @class ParetoLexicase
     * @brief ParetoLexicase selection operator.
     */
    struct ParetoLexicase : SelectionOperator
    {
        ParetoLexicase(bool surv);
        
        ~ParetoLexicase();

        // function returns a set of selected indices from pop 
        vector<size_t> select(Population& pop,  
                const Parameters& params, const Data& d); 
        
        // lexicase survival
        vector<size_t> survive(Population& pop,  
                const Parameters& params, const Data& d); 
          
        // number of test cases used to select each of the selected individuals
        vector<size_t> n_cases_used;

        private:

            //< dominance comparison with epison relaxation
            auto epsi_dominated(vector<float> const &lhs, 
                vector<float> const &rhs, vector<float> const &eps)
                const -> eDominance;

            // pareto front of rank 0 using Fast non-dominated sorting, based
            // on the epsilon for test case t (eps_t) and epsilon for
            // population complexity (eps_c)
            void fast_eNDS(
                vector<Individual>& individuals, vector<size_t>& pool,
                size_t case_id, vector<float> eps);  

            // the Pareto front (just individuals in the rank 1). Used by fast_eNDS
            vector<size_t> front;   
    };
}
}

#endif