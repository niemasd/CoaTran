#ifndef COALESCENT_H
#define COALESCENT_H
#include <vector>
using namespace std;

/**
 * Sample a coalescent tree under constant effective population size
 * @param eff_pop_size Effective population size (if constant)
 * @param eff_pop_growth The effective population size growth rate (if exponential growth)
 * @param seed The seed individual whose transmission chain phylogeny we want to build
 * @param infection_time Each person's infection time
 * @param infected Keep track of the people person x infected
 * @param sample_times Keep track of each person's sample time(s)
 * @param phylo The (initially empty) vector to fill with the Newick tree as <left,right,time,person> tuples
 * @return The node (as an index of phylo) corresponding to the root of the (sub)tree
 */
int coalescent(
#if defined EXPGROWTH   // exponential effective population size growth
    double const & init_eff_pop_size, double const & eff_pop_growth,
#elif defined TRANSTREE // latest possible coalescence (time of transmission)
    // no parameters needed
#elif defined INFTIME   // earliest possible coalescence (time of infection)
    // no parameters needed
#else                   // constant effective population size
    double const & eff_pop_size,
#endif
int const & seed, vector<double> const & infection_time, vector<vector<int>> const & infected, vector<vector<double>> const & sample_times, vector<tuple<int,int,double,int>> & phylo);
#endif