#ifndef BW_H
#define BW_H
#include "read.h"
#include <string>
#include <vector>
#include <cmath>
#include <map>
class BW_OUT{
public:
	bool converged 	= 0;
	double LL 		= 0;
	vector<double> W;
	vector<vector<double>>A;
	BW_OUT(bool , double , vector<double> , double ** );
	BW_OUT();
};

BW_OUT runBW(map<string,contig *>, vector<double>, double, double, double, bool, int, int, bool );
void GAMMA(double ** , double ** , int ,double **);
void backward(double ** ,  double ** , int ,double ** );
void forward(double ** , double ** , int , double **  );
void emissions(contig *, vector<double> , int , double **, bool);

#endif