#include <cstdlib>
#include <iostream>
#include <string.h>
#include "coalescent.h"
#include "common.h"
using namespace std;

// CoaTran version
#ifndef COATRAN_VERSION
#define COATRAN_VERSION "0.0.1"
#endif

// RNG seed environment variable
#ifndef RNG_SEED_ENV_VAR
#define RNG_SEED_ENV_VAR "COATRAN_RNG_SEED"
#endif

// description
#ifndef DESCRIPTION
#define DESCRIPTION string("CoaTran v") + string(COATRAN_VERSION)
#endif

// opening message
#ifndef OPEN_MESSAGE
const string OPEN_MESSAGE = DESCRIPTION + string(
#if defined EXPGROWTH   // exponential effective population size growth
" (exponential effective population size growth)"
#elif defined TRANSTREE // latest possible coalescence (time of transmission)
" (time of transmission)"
#elif defined INFTIME   // earliest possible coalescence (time of infection)
" (time of infection)"
#else                   // constant effective population size
" (constant effective population size)"
#endif
);
#endif

// number of user args
#ifndef NUM_USER_ARGS
#if defined EXPGROWTH   // exponential effective population size growth
#define NUM_USER_ARGS 5
#elif defined TRANSTREE // latest possible coalescence (time of transmission)
#define NUM_USER_ARGS 3
#elif defined INFTIME   // earliest possible coalescence (time of infection)
#define NUM_USER_ARGS 3
#else                   // constant effective population size
#define NUM_USER_ARGS 4
#endif
#endif

// declare extern global vars from common.h
vector<double> infection_time;
unordered_map<string,int> name2num;
vector<string> num2name;
vector<tuple<int,int,double,int>> phylo;
vector<int> seeds;
vector<vector<int>> infected;
vector<vector<double>> sample_times;

// declare extern global vars from coalescent.h
#if defined EXPGROWTH   // exponential effective population size growth
double init_eff_pop_size;
double eff_pop_growth;
#elif defined TRANSTREE // latest possible coalescence (time of transmission)
// no parameters needed
#elif defined INFTIME   // earliest possible coalescence (time of infection)
// no parameters needed
#else                   // constant effective population size
double eff_pop_size;
#endif

// main driver
int main(int argc, char** argv) {
    // check usage
    if(argc != NUM_USER_ARGS || strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0) {
        cerr << OPEN_MESSAGE << endl << "USAGE: " << argv[0] << " <trans_network> <sample_times>"
        #if defined EXPGROWTH   // exponential effective population size growth
            << " <init_eff_pop_size> <eff_pop_growth>"
        #elif defined TRANSTREE // latest possible coalescence (time of transmission)
            // no parameters needed
        #elif defined INFTIME   // earliest possible coalescence (time of infection)
            // no parameters needed
        #else                   // constant effective population size
            << " <eff_pop_size>"
        #endif
        << endl; exit(1);
    }

    // check if user provided a seed
    const char* const rng_seed_env = getenv(RNG_SEED_ENV_VAR);
    if(rng_seed_env != nullptr) {
        int tmp = atoi(rng_seed_env);
        if(tmp != 0) {
            RNG_SEED = tmp; RNG = default_random_engine(RNG_SEED);
        }
    }

    // check if files exist
    if(!file_exists(argv[1])) {
        cerr << "File not found: " << argv[1] << endl; exit(1);
    }
    if(!file_exists(argv[2])) {
        cerr << "File not found: " << argv[2] << endl; exit(1);
    }

    // parse parameter(s)
    #if defined EXPGROWTH   // exponential effective population size growth
        init_eff_pop_size = atof(argv[3]);
        eff_pop_growth = atof(argv[4]);
    #elif defined TRANSTREE // latest possible coalescence (time of transmission)
        // no parameters needed
    #elif defined INFTIME   // earliest possible coalescence (time of infection)
        // no parameters needed
    #else                   // constant effective population size
        eff_pop_size = atof(argv[3]);
    #endif

    // parse transmission network
    parse_transmissions(argv[1]);
    const unsigned int NUM_PEOPLE = num2name.size();
    const unsigned int NUM_SEEDS = seeds.size();

    // parse sample times
    sample_times = vector<vector<double>>(NUM_PEOPLE, vector<double>());
    parse_sample_times(argv[2]);

    // sample coalescent phylogenies; phylos[i] is a vector of <left,right,time,person> nodes for seed i
    vector<vector<tuple<int,int,double,int>>> phylos(NUM_SEEDS, vector<tuple<int,int,double,int>>());
    vector<int> roots(NUM_SEEDS, -1); // roots[i] is the root index of phylos[i]
    for(unsigned int i = 0; i < NUM_SEEDS; ++i) {
        roots[i] = coalescent(seeds[i], phylos[i]);
    }

    // output Newick strings for each phylogeny
    for(unsigned int i = 0; i < NUM_SEEDS; ++i) {
        int const & root = roots[i]; vector<tuple<int,int,double,int>> const & phylo = phylos[i];
        if(!phylo.empty()) {
            string s; newick(root, phylo, s); s += ':'; s += to_string(get<2>(phylo[root])); s += ';';
            cout << s << endl;
        }
    }
    return 0;
}
