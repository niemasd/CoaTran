#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <sys/stat.h>
#include "common.h"

// initialize extern variables from common.h
int RNG_SEED = chrono::system_clock::now().time_since_epoch().count();
default_random_engine RNG(RNG_SEED);
uniform_real_distribution<double> UNIFORM_0_1(0., 1.);
const double DOUBLE_INFINITY = numeric_limits<double>::infinity();

bool file_exists(char* const & fn) {
    struct stat tmp;
    return (stat(fn, &tmp) == 0);
}

double sample_expon(double const & rate) {
    // if rate is 0, return infinity
    if(rate < ZERO_TOLERANCE_RATE) {
        return DOUBLE_INFINITY;
    }

    // otherwise, sample from exponential r.v.
    const double & P = UNIFORM_0_1(RNG);
    return (-log(1.-P))/rate;
}

double sample_trunc_expon(double const & rate, double const & T) {
    // if rate is 0, return truncation point
    if(rate < ZERO_TOLERANCE_RATE) {
        return T;
    }

    // otherwise, sample from truncated exponential r.v.
    const double & P = UNIFORM_0_1(RNG);
    return (-log(1.-(P*(1.-exp((-rate)*T)))))/rate;
}

double sample_coal_time_expgrowth(double const & tau, int const & N, double const & tauI, double const & S0, double const & r) {
    // if initial effective population size is 0, return 0
    if(S0 < ZERO_TOLERANCE_S0) {
        return 0;
    }

    // if growth rate is 0, return what I do with constant effective population size
    if(r < ZERO_TOLERANCE_RATE) {
        return sample_trunc_expon(N*(N-1)/(2*S0), tau-tauI);
    }

    // otherwise, sample from exponential population growth distribution
    const double & P = UNIFORM_0_1(RNG);
    cout << ((log((-2.*r*S0*log(1.-P))/(N*(N-1)))+1.)/r)+tau-tauI << endl;
    return ((log((-2.*r*S0*log(1.-P))/(N*(N-1)))+1.)/r)+tau-tauI;
}

double sample_coal_time_expgrowth_trunc(double const & tau, int const & N, double const & tauI, double const & S0, double const & r) {
    // if initial effective population size is 0, return 0
    if(S0 < ZERO_TOLERANCE_S0) {
        return 0;
    }

    // if growth rate is 0, return what I do with constant effective population size
    if(r < ZERO_TOLERANCE_RATE) {
        return sample_trunc_expon(N*(N-1)/(2*S0), tau-tauI);
    }

    // otherwise, sample from truncated exponential population growth distribution
    double const & T = tau - tauI; // truncation time
    const double & P = UNIFORM_0_1(RNG);
    return ((log((-2.*r*S0*log(1.-(P*(1.-exp((N*(N-1)*exp(r*(T+tauI-tau)-1.))/(-2.*r*S0))))))/(N*(N-1)))+1)/r)+tau-tauI;
}

void parse_transmissions(char* const & fn, vector<string> & num2name, unordered_map<string,int> & name2num, vector<int> & seeds, vector<double> & infection_time, vector<vector<int>> & infected) {
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
        } else {
            auto itr = name2num.find(tmp);
            if(itr == name2num.end()) {
                v = num2name.size(); num2name.push_back(tmp); name2num[tmp] = v; infected.push_back({});
            } else if(u != itr->second) { // ignore recovery events
                cout << "Reinfection event: " << tmp << endl; exit(1);
            }
        }

        // parse t
        getline(is, tmp, '\n'); double t = stof(tmp);

        // add transmission
        infection_time.push_back(t);
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

void newick(int const & root, vector<tuple<int,int,double,int>> const & phylo, vector<string> const & num2name, string & s) {
    // store node values for convenience
    tuple<int,int,double,int> const & node = phylo[root];
    int const & left = get<0>(node);
    int const & right = get<1>(node);
    double const & time = get<2>(node);
    int const & person = get<3>(node);
    if(time < 0) {
        cerr << "Encountered negative time" << endl; exit(1);
    }

    // if leaf, output NODE|PERSON|TIME
    if(left == -1 && right == -1) {
        if(person == -1) {
            cerr << "Encountered a leaf not associated with a person" << endl; exit(1);
        }
        s += to_string(root); s += "|"; s += num2name[person]; s += "|"; s += to_string(time);
    }

    // if dummy transmission event node, output unifurcation
    else if(left == right) {
        s += "(";
        newick(left, phylo, num2name, s); // child subtree
        s += ":"; s += to_string(get<2>(phylo[left]) - time); // child branch length
        s += ")";
    }

    // if internal node, don't output any label
    else {
        s += "(";
        newick(left, phylo, num2name, s); // left subtree
        s += ":"; s += to_string(get<2>(phylo[left]) - time); // left branch length
        s += ",";
        newick(right, phylo, num2name, s); // right subtree
        s += ":"; s += to_string(get<2>(phylo[right]) - time); // right branch length
        s += ")";
    }
}
