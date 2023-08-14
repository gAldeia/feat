/* FEWTWO
copyright 2017 William La Cava
license: GNU/GPL v3
*/

#include "pareto_lexicase.h"

namespace FT{
namespace Sel{

ParetoLexicase::ParetoLexicase(bool surv)
{ 
    name = "pareto_lexicase"; 
    survival = surv; 
}

ParetoLexicase::~ParetoLexicase(){}

vector<size_t> ParetoLexicase::select(Population& pop,  
        const Parameters& params, const Data& d)
{
    /*! Selection according to lexicase selection for 
     * binary outcomes and epsilon-lexicase selection for continuous. 
     * @param pop: population
     * @param params: parameters.
     *
     * @return selected: vector of indices corresponding to pop that 
     * are selected.
     *
     */
    
    vector<size_t> pool;
    pool.resize(1);
    pool[0] = 0;
    
    return pool;
}

vector<size_t> ParetoLexicase::survive(Population& pop, 
        const Parameters& params, const Data& d)
{
    /* ParetoLexicase survival */
    THROW_RUNTIME_ERROR("Lexicase survival not implemented");
    return vector<size_t>();
}

}
}