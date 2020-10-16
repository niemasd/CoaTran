#ifndef COMMON_H
#define COMMON_H
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

/**
 * Check if a file exists
 * @param fn The filename to check
 * @return `true` if `fn` exists, otherwise `false`
 */
bool file_exists(char* const & fn);

/**
 * Load the transmission network from file
 * @param name2num A mapping from transmission network names to integer labels; need to populate
 * @param infected_by Keep track of who person x was infected by (-1 if seed) and the time <infector,time>; need to populate
 * @param infected Keep track of the people person x infected; need to populate
 */
void parse_transmissions(char* const & fn, vector<string> & num2name, unordered_map<string,int> & name2num, vector<int> & seeds, vector<pair<int,double>> & infected_by, vector<vector<int>> & infected);

/**
 * Load the sample times from file
 * @param name2num An already-filled mapping of names to numbers
 * @param sample_times Keep track of the sample times of person x; need to populate
 */
void parse_sample_times(char* const & fn, unordered_map<string,int> const & name2num, vector<vector<double>> & sample_times);
#endif