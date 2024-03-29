#ifndef COMMON_H
#define COMMON_H
#include <chrono>
#include <limits>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

// define 0 tolerance for Poisson rates
#ifndef ZERO_TOLERANCE_RATE
#define ZERO_TOLERANCE_RATE 0.00000000001
#endif

// define 0 tolerance for times
#ifndef ZERO_TOLERANCE_TIME
#define ZERO_TOLERANCE_TIME 0.0000001
#endif

// define 0 tolerance for effective population sizes (S0)
#ifndef ZERO_TOLERANCE_S0
#define ZERO_TOLERANCE_S0 0.00000000001
#endif

// global variables related to coalescent
extern vector<double> infection_time;           // Each person's infection time
extern unordered_map<string,int> name2num;      // Map names to integers
extern vector<string> num2name;                 // Map integers to names
extern vector<tuple<int,int,double,int>> phylo; // Newick tree as <left,right,time,person> tuples
extern vector<int> seeds;                       // Seed individuals (as integers)
extern vector<vector<int>> infected;            // The individuals infected by a given individual
extern vector<vector<double>> sample_times;     // Keep track of each person's sample time(s)

// random number generation
extern int RNG_SEED;
extern default_random_engine RNG;
extern uniform_real_distribution<double> UNIFORM_0_1;

// infinity
extern const double DOUBLE_INFINITY;

/**
 * Check if a file exists
 * @param fn The filename to check
 * @return `true` if `fn` exists, otherwise `false`
 */
bool file_exists(char* const & fn);

/**
 * Sample from an exponential distribution
 * @param rate The rate parameter (lambda) of the exponential distribution
 * @return A random sample from the user-defined exponential distribution
 */
double sample_expon(double const rate);

/**
 * Sample from a truncated exponential distribution
 * @param rate The rate parameter (lambda) of the truncated exponential distribution
 * @param T The point at which to truncate the exponential distribution, i.e., the maximum sample value
 * @return A random sample from the user-defined truncated exponential distribution
 */
double sample_trunc_expon(double const rate, double const T);

/**
 * Sample from the probability distribution of coalescent time with exponential population growth
 * @param tau Current time
 * @param N Number of lineages at current time tau
 * @param tauI Time of infection
 * @param S0 Initial effective population size
 * @param r Growth rate
 * @return A random sample of a coalescent time under exponential effective population growth
 */
double sample_coal_time_expgrowth(double const tau, int const N, double const tauI, double const S0, double const r);

/**
 * Sample from the probability distribution of truncated coalescent time with exponential population growth
 * Note that T (the truncation time) = tau - tauI
 * @param tau Current time
 * @param N Number of lineages at current time tau
 * @param tauI Time of infection
 * @param S0 Initial effective population size
 * @param r Growth rate
 * @return A random sample of a truncated coalescent time under exponential effective population growth
 */
double sample_coal_time_expgrowth_trunc(double const tau, int const N, double const tauI, double const S0, double const r);

/**
 * Load the transmission network from file
 * @param fn The filename of the transmission network (TSV)
 */
void parse_transmissions(char* const & fn);

/**
 * Load the sample times from file
 * @param fn The filename of the sample times (TSV)
 */
void parse_sample_times(char* const & fn);

/**
 * Build a Newick string from a phylo vector of <left,right,time,person> tuples
 * @param root The root of the subtree
 * @param phylo The phylo vector
 * @param s The string to build
 * @return A Newick string of the tree represented by phylo
 */
void newick(int const root, vector<tuple<int,int,double,int>> const & phylo, string & s);

/**
 * Pop a random element from an unsorted vector
 * @param vec The unsorted vector from which to pop
 * @return A random element, after removing it from the vector (order will change)
 */
template<class T>
T vector_pop(vector<T> & vec) {
    const int last_ind = vec.size() - 1;
    uniform_int_distribution<int> uniform_rv(0, last_ind);
    const int ind_to_remove = uniform_rv(RNG);
    T tmp = vec[ind_to_remove];
    vec[ind_to_remove] = vec[last_ind];
    vec.pop_back();
    return tmp;
}
#endif
