#include <iostream>
#include <string.h>
#include "coalescent.h"
#include "common.h"
using namespace std;

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
#define NUM_USER_ARGS 4
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
            " <eff_pop_growth>"
        #else // constant effective population size
            " <eff_pop_size>"
        #endif
        << endl; exit(1);
    }

    // check if files exist
    if(!file_exists(argv[1])) {
        cerr << "File not found: " << argv[1] << endl; exit(1);
    }
    if(!file_exists(argv[2])) {
        cerr << "File not found: " << argv[2] << endl; exit(1);
    }

    // parse transmission network
    vector<string> num2name;
    unordered_map<string,int> name2num;
    vector<int> seeds;
    vector<pair<int,double>> infected_by;
    vector<vector<int>> infected;
    parse_transmissions(argv[1], num2name, name2num, seeds, infected_by, infected);
    const unsigned int NUM_PEOPLE = num2name.size();

    // parse sample times
    vector<vector<double>> sample_times(NUM_PEOPLE, vector<double>());
    unsigned int tmp_uint;
    parse_sample_times(argv[2], name2num, sample_times, tmp_uint);
    const unsigned int NUM_LEAVES = tmp_uint; // number of leaves in phylogeny
    const unsigned int NUM_NODES = 2 * NUM_LEAVES - 1;

    // sample coalescent phylogeny; tree_pointers[u] = <parent, left_child, right_child> of node u
    vector<tuple<int,int,int>> tree_pointers(NUM_NODES, make_tuple(-1,-1,-1));
    #ifdef EXPGROWTH // exponential effective population size
        coalescent_expgrowth(atof(argv[3]), tree_pointers);
    #else // constant effective population size
        coalescent_constant(atof(argv[3]), tree_pointers);
    #endif
    return 0;
}