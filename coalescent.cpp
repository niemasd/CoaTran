#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <tuple>
#include "coalescent.h"
#include "common.h"
#include <unordered_set> // TODO REMOVE

// global variables related to just implementing coalescent
vector<int> coalescent_root;

// helper iterative post-order traversal
vector<int> postorder(int const & seed) {
    // prep for iterative post-order traversal
    unsigned int const & MAX_DEPTH = 2*num2name.size(); // max possible depth
    vector<int> s1; s1.reserve(MAX_DEPTH);              // first stack of iterative post-order
    vector<int> s2; s2.reserve(MAX_DEPTH);              // second stack of iterative post-order
    vector<int> out; out.reserve(MAX_DEPTH);            // output ordering of individuals
    unordered_set<int> visited;

    // perform iterative post-order traversal and return output
    s1.push_back(seed);
    visited.insert(seed);
    while(!s1.empty()) {
        int const & curr = s1.back(); s1.pop_back(); s2.push_back(curr);
        for(int const & child : infected[curr]) {
            s1.push_back(child);
            visited.insert(child);
        }
    }
    while(!s2.empty()) {
        out.push_back(s2.back()); s2.pop_back();
    }
    return out;
}

// run the actual logic of the coalescent
void coalescent_logic(int const & seed, vector<tuple<int,int,double,int>> & phylo) {
    // store things that are used multiple times
    double const & SEED_INF_TIME = infection_time[seed];

    // add node(s) for sample time(s) of the seed
    vector<int> leaves; // vector of phylo indices of leaves of this segment
    for(double const & t : sample_times[seed]) {
        leaves.push_back(phylo.size()); phylo.push_back(make_tuple(-1,-1,t,seed));
    }

    // first check that this has already been called on children
    for(int const & child : infected[seed]) {
        if(coalescent_root[child] == -1) {
            if(!sample_times[child].empty()) {
                cerr << "Coalescent not run in post-order" << endl;
                cerr << "parent: " << seed << " (" << num2name[seed] << ")" << endl;
                cerr << "child: " << child << " (" << num2name[child] << ")" << endl;
                exit(1);
            }
        } else {
            leaves.push_back(coalescent_root[child]);
        }
    }

    // if no leaves, nothing to do
    if(leaves.empty()) {
        return;
    }

    // sort leaves in decreasing order of time
    sort(leaves.begin(), leaves.end(), [&phylo](int const & lhs, int const & rhs){return get<2>(phylo[lhs]) > get<2>(phylo[rhs]);});

    // precompute values that will be repeatedly used
    #if defined EXPGROWTH   // exponential effective population size growth
        // can precompute -2*r*S0, but not worth the complicated function (doesn't save much))
    #elif defined TRANSTREE // latest possible coalescence (time of transmission)
        // no precomputed values needed
    #elif defined INFTIME   // earliest possible coalescence (time of infection)
        // no precomputed values needed
    #else                   // constant effective population size
        const double TWO_TIMES_C = 2 * eff_pop_size;
    #endif

    // coalesce leaves
    vector<int> lineages = {leaves[0]}; double curr_time = -1;
    for(unsigned int i = 1; i < leaves.size(); ++i) {
        // prepare for coalescing
        lineages.push_back(leaves[i]); // add the next leaf
        curr_time = get<2>(phylo[leaves[i]]); // move time to next leaf

        // if we've added the last lineage, just break and do truncated coalescence
        if(i == leaves.size()-1) {
            break;
        }

        // coalesce as much as possible before time of next next leaf
        while(lineages.size() != 1) {
            // sample the time of the next coalescent event
            double const & coal_time = curr_time
            #if defined EXPGROWTH   // exponential effective population size growth
                - sample_coal_time_expgrowth(curr_time, lineages.size(), SEED_INF_TIME, init_eff_pop_size, eff_pop_growth)
            #elif defined TRANSTREE // latest possible coalescence (time of transmission)
                // do nothing
            #elif defined INFTIME   // earliest possible coalescence (time of infection)
                - DOUBLE_INFINITY
            #else                   // constant effective population size
                - sample_expon(lineages.size()*(lineages.size()-1)/TWO_TIMES_C)
            #endif
            ;

            // if next coalescent event is earlier than next leaf, failed to coalesce
            double const & cutoff_time = get<2>(phylo[leaves[i+1]]);
            if(coal_time < cutoff_time) {
                curr_time = cutoff_time; break;
            }

            // coalesce 2 random lineages
            const int & parent = phylo.size();
            const int & lin1 = vector_pop(lineages);
            const int & lin2 = vector_pop(lineages);
            phylo.push_back(make_tuple(lin1,lin2,coal_time,-1));
            lineages.push_back(parent); curr_time = coal_time;
        }
    }

    // coalesce remaining lineages, constrained to coalesce between curr_time and infection_time[seed]
    while(lineages.size() != 1) {
        // check for validity
        if(curr_time < 0) {
            cerr << "Negative curr_time" << endl; exit(1);
        }
        double coal_time;

        // check if we've hit the seed infection time
        if(abs(curr_time-SEED_INF_TIME) < ZERO_TOLERANCE_TIME) {
            coal_time = SEED_INF_TIME;
        }

        // if not, sample delta under truncated distribution
        else {
            coal_time =
            #if defined EXPGROWTH   // exponential effective population size growth
                curr_time - sample_coal_time_expgrowth_trunc(curr_time, lineages.size(), SEED_INF_TIME, init_eff_pop_size, eff_pop_growth)
            #elif defined TRANSTREE // latest possible coalescence (time of transmission)
                curr_time
            #elif defined INFTIME   // earliest possible coalescence (time of infection)
                SEED_INF_TIME
            #else                   // constant effective population size
                curr_time - sample_trunc_expon(lineages.size()*(lineages.size()-1)/TWO_TIMES_C, curr_time-SEED_INF_TIME)
            #endif
            ;
        }

        // coalesce 2 random lineages
        const int & parent = phylo.size();
        const int & lin1 = vector_pop(lineages);
        const int & lin2 = vector_pop(lineages);
        phylo.push_back(make_tuple(lin1,lin2,coal_time,-1));
        lineages.push_back(parent); curr_time = coal_time;
    }

    // add dummy root node at time of transmission
    const int & parent = phylo.size();
    const int & child = lineages[0];
    phylo.push_back(make_tuple(child,child,SEED_INF_TIME,-1));
    coalescent_root[seed] = parent;
}

// organize how coalescent is run (to avoid recursion)
int coalescent(int const & seed, vector<tuple<int,int,double,int>> & phylo) {
    // prepare for iterative post-order traversal
    /*
    if(coalescent_root.empty()) {
        coalescent_root = vector<int>(num2name.size(), -1);
    }

    // run the coalescent logic on each individual in post-order
    for(int const & curr : postorder(seed)) {
        coalescent_logic(curr, phylo);
    }
    */
    if(coalescent_root.empty()) {
        coalescent_root = vector<int>(num2name.size(), -1);
        for(int curr = num2name.size()-1; curr >= 0; --curr) {
            coalescent_logic(curr, phylo);
        }
    }

    // finished running coalescent on all subtrees, so return the overall root
    return coalescent_root[seed];
}
