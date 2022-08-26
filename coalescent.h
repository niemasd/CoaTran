#ifndef COALESCENT_H
#define COALESCENT_H
#include <vector>
using namespace std;

// global variables related to coalescent
#if defined EXPGROWTH   // exponential effective population size growth
extern double init_eff_pop_size; // Initial effective population size
extern double eff_pop_growth;    // Effective population size growth rate
#elif defined TRANSTREE // latest possible coalescence (time of transmission)
// no parameters needed
#elif defined INFTIME   // earliest possible coalescence (time of infection)
// no parameters needed
#else                   // constant effective population size
extern double eff_pop_size;      // Effective population size
#endif

/**
 * Sample a coalescent tree under constant effective population size
 * @param seed The seed individual whose transmission chain phylogeny we want to build
 * @param phylo The (initially empty) vector to fill with the Newick tree as <left,right,time,person> tuples
 * @return The node (as an index of phylo) corresponding to the root of the (sub)tree
 */
int coalescent(int const seed, vector<tuple<int,int,double,int>> & phylo);
#endif
