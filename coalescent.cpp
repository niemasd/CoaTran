#include <algorithm>
#include <iostream>
#include <random>
#include <tuple>
#include "coalescent.h"
#include "common.h"

int coalescent(
#ifdef EXPGROWTH // exponential effective population size
    double const & eff_pop_growth
#else // constant effective population size
    double const & eff_pop_size
#endif
, int const & seed, vector<double> const & infection_time, vector<vector<int>> const & infected, vector<vector<double>> const & sample_times, vector<tuple<int,int,double,int>> & phylo) {

    // add node(s) for sample time(s) of the seed
    vector<int> leaves; // vector of phylo indices of leaves of this segment
    for(double const & t : sample_times[seed]) {
        leaves.push_back(phylo.size()); phylo.push_back(make_tuple(-1,-1,t,seed));
    }

    // first call this function recursively on children
    for(int const & child : infected[seed]) {
        const int & tmp = coalescent(
        #ifdef EXPGROWTH // exponential effective population size
            eff_pop_growth
        #else // constant effective population size
            eff_pop_size
        #endif
        , child, infection_time, infected, sample_times, phylo);
        if(tmp != -1) {
            leaves.push_back(tmp);
        }
    }

    // if no leaves, nothing to do
    if(leaves.empty()) {
        return -1;
    }

    // if 1 leaf, just return that leaf
    if(leaves.size() == 1) {
        return leaves[0];
    }

    // sort leaves in decreasing order of time
    sort(leaves.begin(), leaves.end(), [&phylo](int const & lhs, int const & rhs){return get<2>(phylo[lhs]) > get<2>(phylo[rhs]);});

    // precompute values that will be repeatedly used
    #ifdef EXPGROWTH // exponential effective population size
        // TODO ADD EXPGROWTH
    #else // constant effective population size
        const double TWO_TIMES_C = 2 * eff_pop_size;
    #endif

    // coalesce leaves
    vector<int> lineages = {leaves[0]}; double curr_time;
    for(unsigned int i = 1; i < leaves.size(); ++i) {
        lineages.push_back(leaves[i]); // add the next leaf
        curr_time = get<2>(phylo[leaves[i]]); // move time to next leaf
        // coalesce as much as possible before time of next next leaf
        while(lineages.size() != 1) {
            // sample the time of the next coalescent event
            const int & N = lineages.size();
            #ifdef EXPGROWTH // exponential effective population size
                exponential_distribution<double> rv(N); // TODO REPLACE WITH CORRECT ONE FOR EXP GROWTH
            #else // constant effective population size
                exponential_distribution<double> rv(N*(N-1)/TWO_TIMES_C);
            #endif
            double const & coal_time = curr_time - rv(RNG);

            // if next coalescent event is earlier than next leaf, failed to coalesce
            double cutoff_time;
            if(i == leaves.size()-1) {
                cutoff_time = infection_time[seed];
            } else {
                cutoff_time = get<2>(phylo[leaves[i+1]]);
            }
            if(coal_time < cutoff_time) {
                break;
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
        // TODO sample delta under truncated distribution
        // TODO coalesce 2 random lineages
        cerr << "NEED TO HANDLE TRUNCATED DISTRIBUTION COALESCENCE" << endl; exit(1);
    }
    return lineages[0];
}