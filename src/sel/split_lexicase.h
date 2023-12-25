/* FEWTWO
copyright 2017 William La Cava
license: GNU/GPL v3
*/

#ifndef SPLITLEXICASE_H
#define SPLITLEXICASE_H

#include "selection_operator.h"


namespace FT{
namespace Sel{
    ////////////////////////////////////////////////////////////// Declarations
    /*!
     * @class SplitLexicase
     * @brief SplitLexicase selection operator.
     */
    struct SplitLexicase : SelectionOperator
    {
        SplitLexicase(bool surv);
        ~SplitLexicase();

        // function returns a set of selected indices from pop 
        vector<size_t> select(Population& pop,  
                const Parameters& params, const Data& d); 
        
        // lexicase survival
        vector<size_t> survive(Population& pop,  
                const Parameters& params, const Data& d); 
          
        // number of test cases used to select each of the selected individuals
        vector<size_t> n_cases_used;

        private:
            /// Uses a heuristic to set a splitting threshold.
            float find_threshold(const ArrayXf& x);

            /// returns the gain of a split 
            float gain(const VectorXf& lsplit, const VectorXf& rsplit);
    };
}
}

#endif