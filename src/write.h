#ifndef write_H
#define write_H
#include <string>
#include <map>
#include "BaumWelch.h"
#include "viterbi.h"
void writeTrainingFile(string ,vector<double>  , vector<vector<double>>  ,
	double  , double  , double  );
void writeViterbiPaths(string , map<string, map<string, vector<segment *> >> );
#endif