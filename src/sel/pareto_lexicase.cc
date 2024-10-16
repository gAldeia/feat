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
    // TODO: write proper docstring
    /*! Selection according to lexicase selection for 
     * binary outcomes and epsilon-lexicase selection for continuous. 
     * @param pop: population
     * @param params: parameters.
     *
     * @return selected: vector of indices corresponding to pop that 
     * are selected.
     */
    
    // set objectives
    #pragma omp parallel for
    for (unsigned int i=0; i<pop.size(); ++i)
        pop.individuals.at(i).set_obj(params.objectives);

    //< number of samples
    unsigned int N = pop.individuals.at(0).error.size(); 
    //< number of individuals
    unsigned int P = pop.individuals.size();
    
    // selection pool
    vector<size_t> starting_pool;
    for (int i = 0; i < pop.individuals.size(); ++i)
    {
        starting_pool.push_back(i);
    }
    assert(starting_pool.size() == P);   

    // selected individuals
    vector<size_t> selected(P,0);
    
    // tracking how many test cases each individual took before being selected
    n_cases_used.resize(P);
    std::iota(n_cases_used.begin(), n_cases_used.end(), 0);
    
    // threshold used to pick each individual
    thresholds.resize(P);
    std::iota(thresholds.begin(), thresholds.end(), 0);

    // selection loop
    #pragma omp parallel for 
    for (unsigned int i = 0; i<P; ++i)
    {
        vector<size_t> cases; // cases (samples)
        if (params.classification && !params.class_weights.empty()) 
        {
            // for classification problems, weight case selection 
            // by class weights
            vector<size_t> choices(N);
            std::iota(choices.begin(), choices.end(),0);
            vector<float> sample_weights = params.sample_weights;
            for (unsigned i = 0; i<N; ++i)
            {
                vector<size_t> choice_idxs(N-i);
                std::iota(choice_idxs.begin(),choice_idxs.end(),0);
                size_t idx = r.random_choice(choice_idxs,
                        sample_weights);
                cases.push_back(choices.at(idx));
                choices.erase(choices.begin() + idx);
                sample_weights.erase(sample_weights.begin() + idx);
            }
        }
        else
        {   // otherwise, choose cases randomly
            cases.resize(N); 
            std::iota(cases.begin(),cases.end(),0);
            r.shuffle(cases.begin(),cases.end());   // shuffle cases
        }

        vector<size_t> pool = starting_pool;    // initial pool   
        vector<size_t> winner;                  // winners
    
        bool pass = true;     // checks pool size and number of cases
        unsigned int h = 0;   // case count, index used in cases[h]

        float error_epsilon;

        while(pass){    // main loop
            // Calculating epsilons on demand
            error_epsilon = 0;

            // if output is continuous, use epsilon lexicase            
            if (!params.classification || params.scorer_.compare("log")==0 
            ||  params.scorer_.compare("multi_log")==0)
            {
                VectorXf pool_error(pool.size());
                for (int j = 0; j<pool.size(); ++j)
                {
                    pool_error(j) = pop.individuals.at(pool[j]).error(cases[h]);
                }
                error_epsilon = mad(pool_error);
            }
            
            float complexity_epsilon;
            VectorXf pool_complexity(pool.size());
            for (int j = 0; j<pool.size(); ++j)
            {
                pool_complexity(j) = 
                    (float) pop.individuals.at(pool[j]).get_complexity();
            }
            complexity_epsilon = mad(pool_complexity);
            
            winner.resize(0);   // winners     

            // fast non-dominated epsilon-sort
            vector<float> eps{error_epsilon, complexity_epsilon};
            fast_eNDS(pop.individuals, pool, cases[h], eps);

            // get winners based on relative index used in pool
            winner.reserve(winner.size() + distance(front.begin(),front.end()));
            winner.insert(winner.end(), front.begin(), front.end());
            
            ++h; // next case 
            // only keep going if needed
            pass = (winner.size()>1 && h<cases.size()); 
            
            if(winner.size() == 0)
            {
                if(h >= cases.size())
                    winner.push_back(r.random_choice(pool));
                else
                    pass = true;
            }
            else
                pool = winner; // reduce pool to remaining individuals
        }       

        assert(winner.size()>0);

        n_cases_used[i] = h;
        thresholds[i]   = error_epsilon;

        //if more than one winner, pick randomly
        selected.at(i) = r.random_choice(winner);   
    }

    if (selected.size() != pop.individuals.size())
    {
        std::cout << "selected: " ;
        for (auto s: selected) std::cout << s << " "; std::cout << "\n";
        THROW_LENGTH_ERROR("Pareto lexicase did not select correct \
                number of parents");
    }

    return selected;
}

vector<size_t> ParetoLexicase::survive(Population& pop, 
        const Parameters& params, const Data& d)
{
    /* ParetoLexicase survival */
    THROW_RUNTIME_ERROR("Lexicase survival not implemented");
    return vector<size_t>();
}

auto ParetoLexicase::epsi_dominated(vector<float> const &lhs, 
    vector<float> const &rhs, vector<float> const &eps) const -> eDominance
{
    bool better { false };
    bool worse { false };

    auto lhs_it = lhs.begin();
    auto rhs_it = rhs.begin();
    auto eps_it = eps.begin();
        
    auto epsi_less = []<bool CheckNan = false, typename T>(T a, T b, T eps)
    {
        if (CheckNan) {
            if (std::isnan(a)) return false;
            if (std::isnan(b)) return true;
        }
        // needs to be smaller, but not equal. 
        return a < b && abs(a-b) > eps;
    };

    for (; lhs_it!=lhs.end() && rhs_it!=rhs.end() && eps_it!=eps.end(); 
        ++lhs_it, ++rhs_it, ++eps_it)
    {
        better |= epsi_less(*lhs_it, *rhs_it, *eps_it);
        worse  |= epsi_less(*rhs_it, *lhs_it, *eps_it);
    }

    if (better && !worse) return eDominance::Left;
    if (worse && !better) return eDominance::Right;

    // not better and not worse -> equal (difference is within eps range)
    return eDominance::None;

    // example to verify the result in a logical statement
    // bool accept = d != Dominance::None;
}

void ParetoLexicase::fast_eNDS(
    vector<Individual>& individuals, vector<size_t>& pool, 
    size_t case_id, vector<float> eps)
{
    front.clear();

    // create vectors of individuals objectives here
    #pragma omp parallel for
    for (int i = 0; i < pool.size(); ++i) {
    
        std::vector<unsigned int> dom;
        int dcount = 0;
    
        // TODO: work with any secundary objective, not just complexity
        Individual& p = individuals.at(pool[i]);
        vector<float> p_obj{p.error(case_id), (float)p.complexity};
        // p.dcounter  = 0;
        // p.dominated.clear();
    
        for (int j = 0; j < pool.size(); ++j)
        {
            Individual& q = individuals.at(pool[j]);
            vector<float> q_obj{q.error(case_id), (float)q.complexity};
        
            // TODO: unify e, epsi, eps (maybe epsilon)
            auto compare = epsi_dominated(p_obj, q_obj, eps);

            if (compare == eDominance::Right) { // p dominates q
                //p.dominated.push_back(j);
                dom.push_back(pool[j]);
            } else if (compare == eDominance::Left) { // q dominates p
                //p.dcounter += 1;
                dcount += 1;
            }
        }
    
        #pragma omp critical
        {
            p.dcounter  = dcount;
            p.dominated.clear();
            p.dominated = dom;
        
            if (p.dcounter == 0) {
                p.set_rank(1);
                front.push_back(pool[i]);
            }
        }
    }
}

} // namespace FT
} // namespace Sel
