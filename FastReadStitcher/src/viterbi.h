#ifndef viterbi_H
#define viterbi_H
#include "read.h"
#include <map>
#include <vector>
#include <string>
class state{
public:
	string ID;
	int t;
	int k;
	int start, stop;
	string chrom;
	double max;
	double prob;
	state * next;
	state * prev;
	state * ptr;
	state();
	state(int , double, int, int, int,string);
};

vector<int> runViterbi(vector<vector<double>> , vector<double> , vector<vector<double>>);
vector<vector<double>> learn_transition_parameters(vector<double> ,vector<vector<double>> , vector<int> );
map<string, map<string, vector<segment *> >> run_viteribi_across(map<string, segment * > , vector<double> , vector<vector<double>> );
#endif
