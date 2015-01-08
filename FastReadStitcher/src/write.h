#ifndef write_H
#define write_H
#include <string>
#include <map>
#include "BaumWelch.h"
#include "viterbi.h"
void writeTrainingFile(string, BW_OUT,double, double, double);
void writeViterbiPaths(string, map<string, state*>, string, string);
#endif