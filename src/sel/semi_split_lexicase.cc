/* FEWTWO
copyright 2017 William La Cava
license: GNU/GPL v3
*/

#include "semi_split_lexicase.h"

namespace FT{
namespace Sel{

SemiDynamicSplitLexicase::SemiDynamicSplitLexicase(bool surv)
{ 
    name = "semi_split_lexicase"; 
    survival = surv; 
}

SemiDynamicSplitLexicase::~SemiDynamicSplitLexicase(){}

vector<size_t> SemiDynamicSplitLexicase::select(Population& pop,  
        const Parameters& params, const Data& d)
{
    // Choose threshold t such that the sum of the variance of the errors in
    // the pool on either side of the threshold is minimized. 

    // TODO: write proper docstring
    
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
    std::fill(n_cases_used.begin(), n_cases_used.end(), 0);

    // threshold used to pick each individual (this is used for logging, but its
    // also static and can be calculated once)
    thresholds.resize(P);
    std::fill(thresholds.begin(), thresholds.end(), 0);

    // used to retrieve the split threshold used in each case (for logging)
    vector<float> split_thresholds(N);
    std::fill(split_thresholds.begin(), split_thresholds.end(), 0);
         
    // define epsilon
    ArrayXf epsilon = ArrayXf::Zero(N);

    // if output is continuous, use epsilon lexicase            
    if (!params.classification || params.scorer_.compare("log")==0 
    ||  params.scorer_.compare("multi_log")==0)
    {
        // for each sample, calculate epsilon
        for (int i = 0; i<epsilon.size(); ++i)
        {
            VectorXf case_errors(pop.individuals.size());
            for (int j = 0; j<pop.individuals.size(); ++j)
            {
                case_errors(j) = pop.individuals.at(j).error(i);
            }

            float split_threshold = find_threshold(case_errors);

            split_thresholds[i] = split_threshold;   
            epsilon(i) = split_threshold;
        }
    }

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
        
        while(pass){    // main loop
            winner.resize(0);   // winners     

            // criteria to stay in pool
            float epsilon_threshold = epsilon[cases[h]];

            // select best
            for (size_t j = 0; j<pool.size(); ++j)
                if (pop.individuals.at(pool[j]).error(cases[h]) 
                        <= epsilon_threshold)
                winner.push_back(pool[j]);    

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
        
        // last threshold value before ending the selection of individual i 
        n_cases_used[i] = h;
        thresholds[i]   = epsilon[cases[--h]];

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

vector<size_t> SemiDynamicSplitLexicase::survive(Population& pop, 
        const Parameters& params, const Data& d)
{
    /* Lexicase survival */
    THROW_RUNTIME_ERROR("Lexicase survival not implemented");
    return vector<size_t>();
}

float SemiDynamicSplitLexicase::find_threshold(const ArrayXf& x)
{
    /* cout << "setting threshold\n"; */
    // for each unique value in x, calculate the reduction in the 
    // heuristic brought about by
    // splitting between that value and the next. 
    // set threshold according to the biggest reduction. 
    vector<float> s = unique(x);
    
    if (s.size() <= 1)
        return 0.0;

    sort(s.begin(), s.end()); 

    vector<int> idx(x.size());
    std::iota(idx.begin(),idx.end(), 0);
    
    Map<ArrayXi> midx(idx.data(),idx.size());
    
    float score      = 0;
    float best_score = 0;
    float threshold  = 0;

    /* cout << "s: " ; */ 
    /* for (auto ss : s) cout << ss << " " ; cout << "\n"; */
    /* cout << "x: " << x << "\n"; */
    /* cout << "y: " << y << "\n"; */
    /* cout << "threshold,score\n"; */

    for (unsigned i =0; i<s.size()-1; ++i)
    {
        // we always handle this as equivalent to arity['f']==0. TODO: confirm that
        float val = (s.at(i) + s.at(i+1)) / 2;
        ArrayXi split_idx = (x < val).select(midx,-midx-1);
        
        /* cout << "split val: " << val << "\n"; */

        // split data
        vector<float> d1, d2; 
        for (unsigned j=0; j< split_idx.size(); ++j)
        {
            if (split_idx(j) <0) // TODO: confirm if we just change y with x
                d2.push_back(x(-1-split_idx(j)));
            else
                d1.push_back(x(split_idx(j)));
        }
        if (d1.empty() || d2.empty())
            continue;

        Map<VectorXf> map_d1(d1.data(), d1.size());  
        Map<VectorXf> map_d2(d2.data(), d2.size());  
        
        /* cout << "d1: " << map_d1.transpose() << "\n"; */
        /* cout << "d2: " << map_d2.transpose() << "\n"; */
        
        // If regression and classification must be handled differently, change here
        score = gain(map_d1, map_d2);

        /* cout << "score: " << score << "\n"; */
        if (score < best_score || i == 0)
        {
            best_score = score;
            threshold  = val;
        }
        /* cout << val << "," << score << "\n"; */
    }

    threshold = std::isinf(threshold)? 
        0 : std::isnan(threshold)? 
        0 : threshold;

        /* cout << "final threshold set to " << threshold */ 
        /*      << " with score " */
        /*      << best_score << "\n"; */
    
    return threshold;
}

float SemiDynamicSplitLexicase::gain(const VectorXf& lsplit, const VectorXf& rsplit)
{
    float lscore, rscore, score;
    
    lscore = variance(lsplit.array())/float(lsplit.size());
    rscore = variance(rsplit.array())/float(rsplit.size());

    score = lscore + rscore; 

    return score;
}

} // namespace FT
} // namespace Sel
