#ifndef COMMON_H
#define COMMON_H
#include <chrono>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

// random number generation
extern const int RNG_SEED;
extern default_random_engine RNG;

/**
 * Check if a file exists
 * @param fn The filename to check
 * @return `true` if `fn` exists, otherwise `false`
 */
bool file_exists(char* const & fn);

/**
 * Load the transmission network from file
 * @param fn The filename of the transmission network (TSV)
 * @param num2name A mapping from integer label to transmission network name; need to populate
 * @param name2num A mapping from transmission network names to integer labels; need to populate
 * @param seeds The seed nodes; need to populate
 * @param infection_time Each person's infection time; need to populate
 * @param infected Keep track of the people person x infected; need to populate
 */
void parse_transmissions(char* const & fn, vector<string> & num2name, unordered_map<string,int> & name2num, vector<int> & seeds, vector<double> & infection_time, vector<vector<int>> & infected);

/**
 * Load the sample times from file
 * @param fn The filename of the sample times (TSV)
 * @param name2num An already-filled mapping of names to numbers
 * @param sample_times Keep track of the sample times of person x; need to populate
 * @param num_sample_times The total number of sample times (sum of sample_times values)
 */
void parse_sample_times(char* const & fn, unordered_map<string,int> const & name2num, vector<vector<double>> & sample_times);

/**
 * Build a Newick string from a phylo vector of <left,right,time,person> tuples
 * @param root The root of the subtree
 * @param phylo The phylo vector
 * @param num2name A mapping from integer label to transmission network name
 * @return A Newick string of the tree represented by phylo
 */
void newick(int const & root, vector<tuple<int,int,double,int>> const & phylo, vector<string> const & num2name, string & s);

/**
 * Pop a random element from an unsorted vector
 * @param vec The unsorted vector from which to pop
 * @return A random element, after removing it from the vector (order will change)
 */
template<class T>
T vector_pop(vector<T> & vec) {
    const int & last_ind = vec.size() - 1;
    uniform_int_distribution<int> uniform_rv(0, last_ind);
    const int & ind_to_remove = uniform_rv(RNG);
    T tmp = vec[ind_to_remove];
    vec[ind_to_remove] = vec[last_ind];
    vec.pop_back();
    return tmp;
}
#endif