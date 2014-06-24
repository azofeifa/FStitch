/*
 * functions.h
 *
 *  Created on: Jan 31, 2014
 *      Author: joeyazo
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_
#include <vector>
#include <string>
#include <numeric>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

class logistic_equation{
public:
	vector<double> 	W;
	double weight	= 1.0;
	int label 		= 0;
	logistic_equation(){
		weight =1;
	};
	logistic_equation(double val,vector<double> Ts, int lbl){
		W 		= Ts;
		weight 	= val;
		label 	= lbl;

	}
	double evaluate(vector<double> xs){
		double sum=0;
		if (xs.size() != W.size()){
			cout<<"Vector Dimensions on Logistic Regression do not match..."<<endl;
			return -1.0;
		}
		for (int i = 0 ; i <xs.size(); i++){
			sum=sum+(xs[i]*W[i]);
		}
		double prob 	= 1.0 / (1.0 + exp(-sum));
		if (label){
			return prob;
		}
		return 1-prob;
	}
};





#endif /* FUNCTIONS_H_ */
