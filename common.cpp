#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "common.h"
#include <iostream>

bool file_exists(char* const & fn) {
    struct stat tmp;
    return (stat(fn, &tmp) == 0);
}

void parse_transmissions(char* const & fn, vector<string> & num2name, unordered_map<string,int> & name2num, vector<int> & seeds, vector<pair<int,double>> & infected_by, vector<vector<int>> & infected) {
    ifstream file(fn); string line; string tmp;
    while(getline(file,line)) {
        // check for empty line and set up stringstream
        if(line.size() == 0 || line[0] == '#' || line[0] == '\n') {
            continue;
        }
        istringstream is(line);

        // parse u
        int u; getline(is, tmp, '\t');
        if(tmp == "None") {
            u = -1;
        } else {
            auto itr = name2num.find(tmp);
            if(itr == name2num.end()) {
                cerr << "Infection from person not previously infected: " << tmp << endl; exit(1);
            } else {
                u = itr->second;
            }
        }

        // parse v
        int v; getline(is, tmp, '\t');
        if(tmp == "None") {
            cerr << "\"None\" cannot get infected" << endl; exit(1);
        } else if(name2num.find(tmp) == name2num.end()) {
            v = num2name.size(); num2name.push_back(tmp); name2num[tmp] = v; infected.push_back({});
        } else {
            cout << "Reinfection event: " << tmp << endl; exit(1);
        }

        // parse t
        getline(is, tmp, '\n'); double t = stof(tmp);

        // add transmission
        infected_by.push_back(make_pair(u,t));
        if(u == -1) {
            seeds.push_back(v);
        } else {
            infected[u].push_back(v);
        }
    }
}

void parse_sample_times(char* const & fn, unordered_map<string,int> const & name2num, vector<vector<double>> & sample_times) {
    ifstream file(fn); string line; string tmp;
    while(getline(file,line)) {
        // check for empty line and set up stringstream
        if(line.size() == 0 || line[0] == '#' || line[0] == '\n') {
            continue;
        }
        istringstream is(line);

        // parse u
        int u; getline(is, tmp, '\t');
        if(tmp == "None") {
            cerr << "\"None\" cannot be sampled" << endl; exit(1);
        } else {
            auto itr = name2num.find(tmp);
            if(itr == name2num.end()) {
                cerr << "Sample time of person not in transmission network: " << tmp << endl; exit(1);
            } else {
                u = itr->second;
            }
        }

        // parse t and add to sample_times
        getline(is, tmp, '\n'); sample_times[u].push_back(stof(tmp));
    }
}