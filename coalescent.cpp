#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <tuple>
#include "coalescent.h"
#include "common.h"

int coalescent(int const & seed, vector<tuple<int,int,double,int>> & phylo) {
    // store things that are used multiple times
    double const & SEED_INF_TIME = infection_time[seed];

    // add node(s) for sample time(s) of the seed
    vector<int> leaves; // vector of phylo indices of leaves of this segment
    for(double const & t : sample_times[seed]) {
        leaves.push_back(phylo.size()); phylo.push_back(make_tuple(-1,-1,t,seed));
    }

    // first call this function recursively on children
    for(int const & child : infected[seed]) {
        const int & tmp = coalescent(child, phylo);
        if(tmp != -1) {
            leaves.push_back(tmp);
        }
    }

    // if no leaves, nothing to do
    if(leaves.empty()) {
        return -1;
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
    return parent;
}
