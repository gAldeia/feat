/* FEWTWO
copyright 2017 William La Cava
license: GNU/GPL v3
*/

#ifndef PARETOLEXICASE_H
#define PARETOLEXICASE_H

#include "selection_operator.h"

namespace FT{
namespace Sel{

    ////////////////////////////////////////////////////////////// Declarations
    /*!
     * @class ParetoLexicase
     * @brief ParetoLexicase selection operator.
     */
    struct ParetoLexicase : SelectionOperator
    {
        ParetoLexicase(bool surv);
        
        ~ParetoLexicase();

        /// function returns a set of selected indices from pop 
        vector<size_t> select(Population& pop,  
                const Parameters& params, const Data& d); 
        
        /// lexicase survival
        vector<size_t> survive(Population& pop,  
                const Parameters& params, const Data& d); 

    };
}
}

#endif