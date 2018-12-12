/* FEAT
copyright 2017 William La Cava
license: GNU/GPL v3
*/

#include "variation.h"

namespace FT{
    
    namespace Vary{

        /// constructor
        Variation::Variation(float cr): cross_rate(cr) {}
                   
        /// update cross rate
        void Variation::set_cross_rate(float cr)
        {
        	cross_rate = cr;
        }
        
        /// return current cross rate
        float Variation::get_cross_rate()
        {
        	return cross_rate;
        }
        
         /// destructor
        Variation::~Variation(){}

        std::unique_ptr<Node> random_node(const NodeVector & v)
        {
           /*!
            * return a random node from a list of nodes.
            */          
            assert(v.size()>0 && " attemping to return random choice from empty vector");
            std::vector<size_t> vi(v.size());
            std::iota(vi.begin(), vi.end(), 0);
            size_t idx = r.random_choice(vi);
            return v.at(idx)->clone();
        }
     
        void Variation::vary(Population& pop, const vector<size_t>& parents, const Parameters& params,
                             const Data& d)
        {
            /*!
             * performs variation on the current population. 
             *
             * @param   pop: current population
             * @param  	parents: indices of population to use for variation
             * @param  	params: feat parameters
             *
             * @return  appends params.pop_size offspring derived from parent variation
             */
            unsigned start= pop.size();
            pop.resize(2*params.pop_size);
            #pragma omp parallel for
            for (unsigned i = start; i<pop.size(); ++i)
            {
                bool pass=false;                      // pass check for children undergoing variation     
       
                while (!pass)
                {
                    Individual child; // new individual
                    child.set_id(params.current_gen*params.pop_size+i-start);           

                    if ( r() < cross_rate)      // crossover
                    {
                        // get random mom and dad 
                        Individual& mom = pop.individuals.at(r.random_choice(parents));
                        Individual& dad = pop.individuals.at(r.random_choice(parents));
                        /* int dad = r.random_choice(parents); */
                        // create child
                        params.msg("crossing " + mom.get_eqn() + " with " + 
                                   dad.get_eqn(), 3);
                       
                        if (params.semantic_xo)
                            pass = semantic_cross(mom, dad, child, params, d);
                        else
                            pass = cross(mom, dad, child, params);
                        
                        params.msg("crossing " + mom.get_eqn() + " with " + 
                               dad.get_eqn() + " produced " + child.get_eqn() + 
                               ", pass: " + std::to_string(pass),3);    

                        child.set_parents({mom, dad});
                    }
                    else                        // mutation
                    {
                        // get random mom
                        Individual& mom = pop.individuals.at(r.random_choice(parents));
                        /* int mom = r.random_choice(parents); */                
                        params.msg("mutating " + mom.get_eqn() + "(" + 
                                mom.program_str() + ")", 3);
                        // create child
                        pass = mutate(mom,child,params);
                        
                        params.msg("mutating " + mom.get_eqn() + " produced " + 
                                child.get_eqn() + ", pass: " + std::to_string(pass),3);
                        child.set_parents({mom});
                    }
                    if (pass)
                    {
                        assert(child.size()>0);
                        assert(pop.open_loc.size()>i-start);
                        params.msg("assigning " + child.program_str() + " to pop.individuals[" + 
                            std::to_string(i) + "] with pop.open_loc[" + std::to_string(i-start) + 
                            "]=" + std::to_string(pop.open_loc[i-start]),3);

                        pop.individuals[i] = child;
                        pop.individuals[i].loc = pop.open_loc[i-start];                   
                    }
                }    
           }
          
           pop.update_open_loc();
        }

        bool Variation::mutate(Individual& mom, Individual& child, const Parameters& params)
        {
            /*!
             * chooses uniformly between point mutation, insert mutation and delete mutation 
             * 
             * @param   mom: parent
             * @param   child: offspring produced by mutating mom 
             * @param   params: parameters
             * 
             * @return  true if valid child, false if not 
             */    

            // make child a copy of mom
            mom.clone(child, false);  
            
            float rf = r();
            if (rf < 1.0/3.0 && child.get_dim() > 1){
                delete_mutate(child,params); 
                assert(child.program.is_valid_program(params.num_features, params.longitudinalMap));
            }
            else if (rf < 2.0/3.0 && child.size() < params.max_size)
            {
                insert_mutate(child,params);
                assert(child.program.is_valid_program(params.num_features, params.longitudinalMap));
            }
            else
            {        
                point_mutate(child,params);
                assert(child.program.is_valid_program(params.num_features, params.longitudinalMap));
            }
     
            // check child depth and dimensionality
            return child.size()>0 && child.size() <= params.max_size 
                    && child.get_dim() <= params.max_dim;
        }

        void Variation::point_mutate(Individual& child, const Parameters& params)
        {
            /*! 1/n point mutation. 
             * @param child: individual to be mutated
             * @param params: parameters 
             * @return modified child
             * */
            params.msg("\tpoint mutation",3);
            float n = child.size(); 
            unsigned i = 0;
            // loop thru child's program
            for (auto& p : child.program)
            {
                if (r() < child.get_p(i)/n)  // mutate p. 
                {
                    params.msg("\t\tmutating node " + p->name, 3);
                    NodeVector replacements;  // potential replacements for p

                    if (p->total_arity() > 0) // then it is an instruction
                    {
                        // find instructions with matching in/out types and arities
                        for (const auto& f: params.functions)
                        {
                            if (f->otype == p->otype &&
                                f->arity['f']==p->arity['f'] && 
                                f->arity['b']==p->arity['b'] &&
                                f->arity['c']==p->arity['c'] &&
                                f->arity['z']==p->arity['z'])
                                replacements.push_back(f->rnd_clone());
                        }
                    }
                    else                    // otherwise it is a terminal
                    {
                        // TODO: add terminal weights here
                        // find terminals with matching output types
                                          
                        for (const auto& t : params.terminals)
                        {
                            if (t->otype == p->otype)
                                replacements.push_back(t->clone());
                        }                                       
                    }
                    // replace p with a random one
                    p = random_node(replacements);  
                }
                ++i; 
            }

        }

        void Variation::insert_mutate(Individual& child, const Parameters& params)
        {        
            /*! insertion mutation. 
             * @param child: indiviudal to be mutated
             * @param params: parameters
             * @return modified child
             * */
            
            params.msg("\tinsert mutation",3);
            float n = child.size(); 
            
            if (r()<0.5 || child.get_dim() == params.max_dim)
            {
                // loop thru child's program
                for (unsigned i = 0; i< child.program.size(); ++i)
                {
                    
                    if (r() < child.get_p(i)/n)  // mutate with weighted probability
                    {
                        params.msg("\t\tinsert mutating node " + child.program[i]->name +
                                   " with probability " + std::to_string(child.get_p(i)) + 
                                   "/" + std::to_string(n), 3);
                        NodeVector insertion;  // inserted segment
                        NodeVector fns;  // potential fns 
                        
                        // find instructions with matching output types and a matching arity to i
                        for (const auto& f: params.functions)
                        { 
                            // find fns with matching output types that take this node type as arg
                            if (f->arity[child.program[i]->otype] > 0 && 
                                    f->otype==child.program[i]->otype )
                            { 
                                // make sure there are satisfactory types in terminals to fill fns' 
                                // args
                                if (child.program[i]->otype=='b')
                                    if (in(params.dtypes,'b') || f->arity['b']==1)
                                        fns.push_back(f->rnd_clone());
                                else if (child.program[i]->otype=='f')
                                    if (f->arity['b']==0 || in(params.dtypes,'b') )
                                        fns.push_back(f->rnd_clone());              
                            }
                        }

                        if (fns.size()==0)  // if no insertion functions match, skip
                            continue;

                        // choose a function to insert                    
                        insertion.push_back(random_node(fns));
                        
                        unsigned fa = insertion.back()->arity['f']; // float arity
                        unsigned ca = insertion.back()->arity['c']; //categorical arity
                        unsigned ba = insertion.back()->arity['b']; // bool arity
                        // decrement function arity by one for child node 
                        if (child.program[i]->otype=='f') --fa;
                        else if (child.program[i]->otype=='c') --ca;
                        else --ba; 
                        
                        // push back new arguments for the rest of the function
                        for (unsigned j = 0; j< fa; ++j)
                            insertion.make_tree(params.functions,params.terminals,0,
                                      params.term_weights,'f',params.ttypes);
                        if (child.program[i]->otype=='f')    // add the child now if float
                            insertion.push_back(child.program[i]->clone());
                        for (unsigned j = 0; j< ca; ++j)
                            insertion.make_tree(params.functions,params.terminals,0,
                                      params.term_weights,'c',params.ttypes);
                        if (child.program[i]->otype=='c')    // add the child now if categorical
                            insertion.push_back(child.program[i]->clone());
                        for (unsigned j = 0; j< ba; ++j)
                            insertion.make_tree(params.functions,params.terminals,0,
                                      params.term_weights,'b',params.ttypes);
                        if (child.program[i]->otype=='b') // add the child now if bool
                            insertion.push_back(child.program[i]->clone());
                        // post-fix notation
                        std::reverse(insertion.begin(),insertion.end());
                       
                        string s; 
                        for (const auto& ins : insertion) s += ins->name + " "; 
                        params.msg("\t\tinsertion: " + s + "\n", 3);
                        NodeVector new_program; 
                        splice_programs(new_program, 
                                        child.program, i, i, 
                                        insertion, size_t(0), insertion.size()-1);
                        child.program=new_program;
                        /* // put the new code in the program */ 
                        /* child.program.erase(child.program.begin()+i); */
                        /* for (unsigned j =i ; i<insertion.size(); ++i) */
                        /*     child.program.push_back(insertion[i]->clone()); */
                        /* /1* child.program.insert(child.program.begin()+i, insertion.begin(), *1/ */ 
                        /* /1*                      insertion.end()); *1/ */
                        i += insertion.size()-1;
                    }
                    /* std::cout << "i: " << i << "\n"; */ 
                }
            }
            else    // add a dimension
            {            
                NodeVector insertion; // new dimension
                insertion.make_program(params.functions, params.terminals, 1,  
                             params.term_weights,1,r.random_choice(params.otypes),
                             params.longitudinalMap, params.ttypes);
                /* child.program.insert(child.program.end(),insertion.begin(),insertion.end()); */
                for (const auto& ip : insertion) 
                    child.program.push_back(ip->clone());
            }
        }

        void Variation::delete_mutate(Individual& child, const Parameters& params)
        {

            /*! deletion mutation. works by pruning a dimension. 
             * @param child: individual to be mutated
             * @param params: parameters  
             * @return mutated child
             * */
            params.msg("\tdeletion mutation",3);
            params.msg("\t\tprogram: " + child.program_str(),3);
            vector<size_t> roots = child.program.roots();
            
            size_t end = r.random_choice(roots,child.p); 
            size_t start = child.program.subtree(end);  
            if (params.verbosity >=3)
            { 
                std::string s="";
                for (unsigned i = start; i<=end; ++i) s+= child.program[i]->name + " ";
                params.msg("\t\tdeleting " + std::to_string(start) + " to " + std::to_string(end) 
                           + ": " + s, 3);
            }    
            child.program.erase(child.program.begin()+start,child.program.begin()+end+1);
            // delete program from start to end by doing a crossover with an empty program at those
            // locations
            /* vector<std::unique_ptr<Node>>blanks; */
            /* splice_programs(child.program, child.program, start, end, */ 
            /*                 blanks,size_t(0),size_t(-1)); */
            params.msg("\t\tresult of delete mutation: " + child.program_str(), 3);
        }
        
        

        bool Variation::cross(Individual& mom, Individual& dad, Individual& child, 
                              const Parameters& params)
        {
            /*!
             * crossover by either subtree crossover or swapping of dimensions. 
             *
             * @param   mom: root parent
             * @param   dad: parent from which subtree is chosen
             * @param   child: result of cross
             * @param   params: parameters
             *
             * @return  child: mom with dad subtree graft
             */
                        
            bool subtree = r() <0.5;     // half the time, do subtree xo. 
                                         // half the time, swap dimensions.
            vector<size_t> mlocs, dlocs; // mom and dad locations for consideration
            size_t i1, j1, i2, j2;       // i1-j1: mom portion, i2-j2: dad portion
            
            if (subtree) 
            {
                params.msg("\tsubtree xo",3);
                // limit xo choices to matching output types in the programs. 
                vector<char> d_otypes;
                for (const auto& p : dad.program)
                    if(!in(d_otypes,p->otype))
                        d_otypes.push_back(p->otype);
                
                // get valid subtree locations
                for (size_t i =0; i<mom.size(); ++i) 
                    if (in(d_otypes,mom.program[i]->otype)) 
                        mlocs.push_back(i);       
                if (mlocs.size()==0)        // mom and dad have no overlapping types, can't cross
                {
                    params.msg("WARNING: no overlapping types between " + mom.program_str() + "," 
                                 + dad.program_str() + "\n", 3);
                    return 0;               
                }
                j1 = r.random_choice(mlocs,mom.get_p(mlocs));    

                // get locations in dad's program that match the subtree type picked from mom
                for (size_t i =0; i<dad.size(); ++i) 
                    if (dad.program[i]->otype == mom.program[j1]->otype) dlocs.push_back(i);
            } 
            else             // half the time, pick a root node
            {
                params.msg("\troot xo",3);
                mlocs = mom.program.roots();
                dlocs = dad.program.roots();
                params.msg("\t\trandom choice mlocs (size "+
                           std::to_string(mlocs.size())+"), p size: "+std::to_string(mom.p.size()),3);
                j1 = r.random_choice(mlocs,mom.get_p(mlocs));   // weighted probability choice    
            }
            /* cout << "mom subtree\t" << mom.program_str() << "\n"; */
            // get subtree              
            i1 = mom.program.subtree(j1);
                                
            /* cout << "dad subtree\n" << dad.program_str() << "\n"; */
            // get dad subtree
            j2 = r.random_choice(dlocs);
            i2 = dad.program.subtree(j2); 
                   
            // make child program by splicing mom and dad
            splice_programs(child.program, mom.program, i1, j1, dad.program, i2, j2 );
                         
            if (params.verbosity >= 3)
                print_cross(mom,i1,j1,dad,i2,j2,child);     

            assert(child.program.is_valid_program(params.num_features, params.longitudinalMap));
            // check child depth and dimensionality
            return child.size()>0 && child.size() <= params.max_size 
                        && child.get_dim() <= params.max_dim;
        }

        /// semantic crossover
        bool Variation::semantic_cross(Individual& mom, Individual& dad, Individual& child,
                            const Parameters& params, const Data& d)
        {
            /*!
             * crossover by either swapping in a dimension most correlated with the residual of mom. 
             *
             * @param   mom: root parent
             * @param   dad: parent from which subtree is chosen
             * @param   child: result of cross
             * @param   params: parameters
             *
             * @return  child: mom with dad subtree graft
             */
                       
            vector<size_t> mlocs, dlocs; // mom and dad locations for consideration
            size_t i1, j1, j1_idx, i2, j2;       // i1-j1: mom portion, i2-j2: dad portion
            
            params.msg("\troot xo",3);
            mlocs = mom.program.roots();
            vector<int> mlocs_indices(mlocs.size());
            std::iota(mlocs_indices.begin(),mlocs_indices.end(),0);

            dlocs = dad.program.roots();
            params.msg("\t\trandom choice mlocs (size "+
                       std::to_string(mlocs.size())+"), p size: "+std::to_string(mom.p.size()),3);
            j1_idx = r.random_choice(mlocs_indices,mom.get_p(mlocs));   // weighted probability choice    
            j1 = mlocs.at(j1_idx); 
            // get subtree              
            i1 = mom.program.subtree(j1);
                                
            // get dad subtree
            // choose root in dad that is most correlated with the residual of
            // j2 = index in dad.Phi that maximizes R2(d.y, w*mom.Phi - w_i1*Phi_i1)
            /* cout << "mom: " << mom.get_eqn() << "\n"; */
            /* cout << "mom yhat: " << mom.yhat.transpose() << "\n"; */
            VectorXd tree = mom.ml->get_weights().at(j1_idx)*mom.Phi.row(j1_idx).array();
            /* cout << "tree (idx=" << j1_idx << "): " << tree.transpose() << "\n"; */
            VectorXd mom_pred_minus_tree = mom.yhat - tree; 
    /* #pragma omp critical */
            /* { */
            /* VectorXd mom_pred_minus_tree = mom.predict_drop(d,params,j1_idx); */ 
            /* } */
            /* cout << "mom_pred_minus_tree: " << mom_pred_minus_tree.transpose() << "\n"; */
            VectorXd mom_residual = d.y - mom_pred_minus_tree;
            /* cout << "mom_residual: " << mom_residual.transpose() << "\n"; */
           
            // get correlations of dad's features with the residual from mom, less the swap choice
            vector<double> corrs(dad.Phi.rows());
            int best_corr_idx = 0;
            double best_corr = 0;
            double corr; 

            for (int i = 0; i < dad.Phi.rows(); ++i)
            {
                corr = pearson_correlation(mom_residual.array(), dad.Phi.row(i).array());
                corrs.at(i) = corr; // this can be removed
                if (corr > best_corr || i == 0 )
                {
                    best_corr_idx = i; 
                    best_corr = corr;
                }
            }
            /* cout << "corrs: "; */
            /* for (auto c : corrs) cout << c << ", "; cout << "\n"; */
            /* cout << "chose corr at " << best_corr_idx << "\n"; */
            j2 = dlocs.at(best_corr_idx);
            i2 = dad.program.subtree(j2); 
            
            if (params.verbosity >= 3)
                print_cross(mom,i1,j1,dad,i2,j2,child,false);     
           
            // make child program by splicing mom and dad
            splice_programs(child.program, mom.program, i1, j1, dad.program, i2, j2 );
            
            if (params.verbosity >= 3)
                print_cross(mom,i1,j1,dad,i2,j2,child,true);     
                 
            
            assert(child.program.is_valid_program(params.num_features, params.longitudinalMap));
            // check child depth and dimensionality
            return child.size()>0 && child.size() <= params.max_size 
                        && child.get_dim() <= params.max_dim;
            
        }
        
        // swap vector subsets with different sizes. 
        void Variation::splice_programs( NodeVector& vnew,
                                         const NodeVector& v1, size_t i1, size_t j1, 
                                         const NodeVector& v2, size_t i2, size_t j2)
        {
            /*!
             * swap vector subsets with different sizes. 
             * constructs a vector made of v1[0:i1], v2[i2:j2], v1[i1:end].
             *
             * @param   v1: root parent 
             * @param       i1: start of splicing segment 
             * @param       j1: end of splicing segment
             * @param   v2: donating parent
             * @param       i2: start of donation
             * @param       j2: end of donation
             *
             * @return  vnew: new vector 
             */
            /* std::cout << "in splice_programs\n"; */
            /* std::cout << "i1: " << i1 << ", j1: " << j1  << ", i2: " << i2 << ", j2:" << j2 << "\n"; */
            // size difference between subtrees  
            for (unsigned i = 0; i < i1 ; ++i)                  // beginning of v1
                vnew.push_back(v1.at(i)->clone());
            for (unsigned i = i2; i <= j2 ; ++i)                 // spliced in v2 portion
                vnew.push_back(v2.at(i)->clone());
            for (unsigned i = j1+1; i < v1.size() ; ++i)        // end of v1
                vnew.push_back(v1.at(i)->clone());
            /* for (const auto& p : vnew) std::cout << p->name << " "; */
            /* std::cout << "\n"; */
            /* vnew.insert(vnew.end(),v1.begin(),v1.begin()+i1); */          
            /* vnew.insert(vnew.end(),v2.begin()+i2,v2.begin()+j2+1);  // spliced in v2 portion */
            /* vnew.insert(vnew.end(),v1.begin()+j1+1,v1.end()); */       
        }
        
        void Variation::print_cross(Individual& mom, size_t i1, size_t j1, Individual& dad, size_t i2, 
                                    size_t j2, Individual& child, bool after)
        {
            std::cout << "\t\tattempting the following crossover:\n\t\t";
            for (int i =0; i<mom.program.size(); ++i){
               if (i== i1) 
                   std::cout << "[";
               std::cout << mom.program[i]->name << " ";
               if (i==j1)
                   std::cout <<"]";
            }
            std::cout << "\n\t\t";
           
            for (int i =0; i<dad.program.size(); ++i){
                if (i== i2) 
                    std::cout << "[";
                std::cout << dad.program[i]->name << " ";
                if (i==j2)
                    std::cout <<"]";
            }
            std::cout << "\n\t\t";
            if (after)
            {
                std::cout << "child after cross: ";
                for (unsigned i = 0; i< child.program.size(); ++i){
                    if (i==i1) std::cout << "[";
                    std::cout << child.program[i]->name << " ";
                    if (i==i1+j2-i2) std::cout << "]";
                }
                std::cout << "\n";
            }
        }
    }
}

