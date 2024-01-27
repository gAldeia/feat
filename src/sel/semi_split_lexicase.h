/* FEWTWO
copyright 2017 William La Cava
license: GNU/GPL v3
*/

#ifndef SemiDynamicSPLITLEXICASE_H
#define SPLITLEXICASE_H

#include "selection_operator.h"


namespace FT{
namespace Sel{
    ////////////////////////////////////////////////////////////// Declarations
    /*!
     * @class SemiDynamicSplitLexicase
     * @brief SemiDynamicSplitLexicase selection operator.
     */
    struct SemiDynamicSplitLexicase : SelectionOperator
    {
        //static: get epsilon first with population, and then create bool array to apply normal lexicase
        //semi: get all epsilon first with population, but criteria is based on eps+error(pool), not population
        //dynamic: epsilon uses the pool
        
        SemiDynamicSplitLexicase(bool surv);
        ~SemiDynamicSplitLexicase();

        // function returns a set of selected indices from pop 
        vector<size_t> select(Population& pop,  
                const Parameters& params, const Data& d); 
        
        // lexicase survival
        vector<size_t> survive(Population& pop,  
                const Parameters& params, const Data& d); 
          
        // number of test cases used to select each of the selected individuals
        vector<size_t> n_cases_used;

        // split or epsilon threshold to remain in the pool before picking the
        // final individual
        vector<float> thresholds;

        private:
            /// Uses a heuristic to set a splitting threshold.
            float find_threshold(const ArrayXf& x);

            /// returns the gain of a split 
            float gain(const VectorXf& lsplit, const VectorXf& rsplit);
    };
}
}

#endif