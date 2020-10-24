#include <cstdlib>
#include <iostream>
#include <string.h>
#include "coalescent.h"
#include "common.h"
using namespace std;

// RNG seed environment variable
#ifndef RNG_SEED_ENV_VAR
#define RNG_SEED_ENV_VAR "COATRAN_RNG_SEED"
#endif

// description
#ifndef DESCRIPTION
#define DESCRIPTION string("CoaTran: Coalescent tree simulation along a transmission network")
#endif

// opening message
#ifndef OPEN_MESSAGE
#ifdef EXPGROWTH // exponential effective population size
const string OPEN_MESSAGE = DESCRIPTION + string(" (exponential effective population size)");
#else // constant effective population size
const string OPEN_MESSAGE = DESCRIPTION + string(" (constant effective population size)");
#endif
#endif

// number of user args
#ifndef NUM_USER_ARGS
#ifdef EXPGROWTH // exponential effective population size
#define NUM_USER_ARGS 5
#else // constant effective population size
#define NUM_USER_ARGS 4
#endif
#endif

// main driver
int main(int argc, char** argv) {
    // check usage
    if(argc != NUM_USER_ARGS || strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0) {
        cerr << OPEN_MESSAGE << endl << "USAGE: " << argv[0] << " <trans_network> <sample_times>" <<
        #ifdef EXPGROWTH // exponential effective population size
            " <init_eff_pop_size> <eff_pop_growth>"
        #else // constant effective population size
            " <eff_pop_size>"
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
    #ifdef EXPGROWTH // exponential effective population size
        double const INIT_EFF_POP_SIZE = atof(argv[3]);
        double const EFF_POP_GROWTH = atof(argv[4]);
    #else // constant effective population size
        double const EFF_POP_SIZE = atof(argv[3]);
    #endif

    // parse transmission network
    vector<string> num2name;
    unordered_map<string,int> name2num;
    vector<int> seeds;
    vector<double> infection_time;
    vector<vector<int>> infected;
    parse_transmissions(argv[1], num2name, name2num, seeds, infection_time, infected);
    const unsigned int NUM_PEOPLE = num2name.size();
    const unsigned int NUM_SEEDS = seeds.size();

    // parse sample times
    vector<vector<double>> sample_times(NUM_PEOPLE, vector<double>());
    parse_sample_times(argv[2], name2num, sample_times);

    // sample coalescent phylogenies; phylos[i] is a vector of <left,right,time,person> nodes for seed i
    vector<vector<tuple<int,int,double,int>>> phylos(NUM_SEEDS, vector<tuple<int,int,double,int>>());
    vector<int> roots(NUM_SEEDS, -1); // roots[i] is the root index of phylos[i]
    for(unsigned int i = 0; i < NUM_SEEDS; ++i) {
        roots[i] = coalescent(
        #ifdef EXPGROWTH // exponential effective population size
            INIT_EFF_POP_SIZE, EFF_POP_GROWTH
        #else // constant effective population size
            EFF_POP_SIZE
        #endif
        , seeds[i], infection_time, infected, sample_times, phylos[i]);
    }

    // output Newick strings for each phylogeny
    for(unsigned int i = 0; i < NUM_SEEDS; ++i) {
        int const & root = roots[i]; vector<tuple<int,int,double,int>> const & phylo = phylos[i];
        if(!phylo.empty()) {
            string s; newick(root, phylo, num2name, s); s += ':'; s += to_string(get<2>(phylo[root])); s += ';';
            cout << s << endl;
        }
    }
    return 0;
}