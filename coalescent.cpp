#include <algorithm>
#include <iostream>
#include <tuple>
#include "coalescent.h"

int coalescent(
#ifdef EXPGROWTH // exponential effective population size
    double const & eff_pop_growth
#else // constant effective population size
    double const & eff_pop_size
#endif
, int const & seed, vector<double> const & infection_time, vector<vector<int>> const & infected, vector<vector<double>> const & sample_times, vector<tuple<int,int,int,double>> & phylo) {

    // add node(s) for sample time(s) of the seed
    cout << "HANDLING NODE " << seed << endl;
    vector<int> leaves; // vector of phylo indices of leaves of this segment
    for(double const & t : sample_times[seed]) {
        leaves.push_back(phylo.size()); phylo.push_back(make_tuple(-1,-1,-1,t));
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
    sort(leaves.begin(), leaves.end(), [&phylo](int const & lhs, int const & rhs){return get<3>(phylo[lhs]) > get<3>(phylo[rhs]);});

    // coalesce leaves
    //double curr_time = get<3>(phylo[leaves[1]]); // start at time of 2nd-latest leaf
    cout << "NEED TO IMPLEMENT THE ACTUAL COALESCENCE!" << endl; exit(1);
    return leaves[0];
}