#ifndef COALESCENT_H
#define COALESCENT_H
#include <vector>
using namespace std;

/**
 * Sample a coalescent tree under constant effective population size
 * @param eff_pop_size Effective population size (constant)
 * @param tree_pointers The tree to build; tree_pointers[u] = <parent, left_child, right_child> of node u
 */
void coalescent_constant(double const & eff_pop_size, vector<tuple<int,int,int>> & tree_pointers);

/**
 * Sample a coalescent tree under exponential effective population growth
 * @param eff_pop_growth The effective population size growth rate
 * @param tree_pointers The tree to build; tree_pointers[u] = <parent, left_child, right_child> of node u
 */
void coalescent_expgrowth(double const & eff_pop_growth, vector<tuple<int,int,int>> & tree_pointers);
#endif