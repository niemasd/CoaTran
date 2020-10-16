#include <iostream>
#include <string.h>
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
        cout << OPEN_MESSAGE << endl << "USAGE: " << argv[0] << " <trans_network> <sample_times>" <<
        #ifdef EXPGROWTH // exponential effective population size
            " <eff_pop_growth>"
        #else // constant effective population size
            " <eff_pop_size>"
        #endif
        << endl; return 1;
    }

    // check if files exist
    if(!file_exists(argv[1])) {
        cout << "File not found: " << argv[1] << endl; return 1;
    }
    if(!file_exists(argv[2])) {
        cout << "File not found: " << argv[2] << endl; return 1;
    }

    // parse transmission network
    vector<string> num2name;
    unordered_map<string,int> name2num;
    vector<int> seeds;
    vector<pair<int,double>> infected_by;
    vector<vector<int>> infected;
    parse_transmissions(argv[1], num2name, name2num, seeds, infected_by, infected);

    //atof(argv[2]) // this parses the 2nd argument as a double
    return 0;
}