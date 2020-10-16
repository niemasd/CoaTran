// include statements
#include <iostream>
#include <string.h>

// description
#ifndef DESCRIPTION
#define DESCRIPTION std::string("CoaTran: Coalescent tree simulation along a transmission network")
#endif

// opening message
#ifndef OPEN_MESSAGE
#ifdef EXPGROWTH // exponential effective population size
const std::string OPEN_MESSAGE = DESCRIPTION + std::string(" (exponential effective population size)");
#else // constant effective population size
const std::string OPEN_MESSAGE = DESCRIPTION + std::string(" (constant effective population size)");
#endif
#endif

// number of user args
#ifndef NUM_USER_ARGS
#ifdef EXPGROWTH // exponential effective population size
#define NUM_USER_ARGS 2
#else // constant effective population size
#define NUM_USER_ARGS 2
#endif
#endif

// main driver
int main(int argc, char* argv[]) {
    if(argc != NUM_USER_ARGS || strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0) {
        std::cout << OPEN_MESSAGE << std::endl;
        #ifdef EXPGROWTH // exponential effective population size
            std::cout << "USAGE: " << argv[0] << " <eff_pop_growth>" << std::endl; return 1;
        #else // constant effective population size
            std::cout << "USAGE: " << argv[0] << " <eff_pop_size>" << std::endl; return 1;
        #endif
    }
    //std::cout << 2*atof(argv[1]) << std::endl;
    return 0;
}