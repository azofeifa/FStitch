/*
 * simuateGMM.h
 *
 *  Created on: Feb 7, 2014
 *      Author: joeyazo
 */

#ifndef SIMUATEGMM_H_
#define SIMUATEGMM_H_

#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <random>
#include "histogram.h"
#include <math.h>
#include <cmath>
using namespace std;


results simulate(int nrolls=1000, int G = 2){
	vector<double> X;
	vector<double> Y;
	vector<double> MUS;
	MUS.push_back(2);
	MUS.push_back(100);

	double rand;
	uniform_real_distribution<double> distribution(0.0,1.0);
	default_random_engine generator;

	for (int j = 0; j < MUS.size(); j++){
		for (int i=0; i < nrolls; i++){
			rand = sqrt(-2*log(distribution(generator))) * cos(2*M_PI*distribution(generator)) ;
			X.push_back((rand * 2.0) + MUS[j]);
			Y.push_back(1);
			}
	}



	return results(X,Y);
}


#endif /* SIMUATEGMM_H_ */
