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

map<string, state *> runViterbi(map<string,contig *> , vector<double> , vector<vector<double>>, int );
#endif
