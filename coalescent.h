#ifndef COALESCENT_H
#define COALESCENT_H
#include <string>
#include <vector>
using namespace std;

/**
 * Sample a coalescent tree under constant effective population size
 * @param eff_pop_size Effective population size (constant)
 * @param seed The seed individual whose transmission chain phylogeny we want to build
 * @param infection_time Each person's infection time
 * @param infected Keep track of the people person x infected
 * @param phylo The (initially empty) string object to fill with the Newick tree
 */
void coalescent_constant(double const & eff_pop_size, int const & seed, vector<double> const & infection_time, vector<vector<int>> const & infected, string & phylo);

/**
 * Sample a coalescent tree under exponential effective population growth
 * @param eff_pop_growth The effective population size growth rate
 * @param seed The seed individual whose transmission chain phylogeny we want to build
 * @param infection_time Each person's infection time
 * @param infected Keep track of the people person x infected
 * @param phylo The (initially empty) string object to fill with the Newick tree
 */
void coalescent_expgrowth(double const & eff_pop_growth, int const & seed, vector<double> const & infection_time, vector<vector<int>> const & infected, string & phylo);
#endif