/*
 * histogram.h
 *
 *  Created on: Feb 6, 2014
 *      Author: joeyazo
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
class results{
public:
	vector<double> edges;
	vector<double> counts;
	results(vector<double> a, vector<double> b){
		edges = a;
		counts = b;
	}
	results(){};
};




results bin(vector<double> X, vector<double> Y, int bin_number, bool expand = 0){
	vector< vector<double> > bins;
	sort(X.begin(), X.end());
	double bin_size = (double)(X[X.size()-1] - X[0]) / (double)bin_number;
	vector<double> empty;
	for (int i = 0; i < bin_number; i++){
		bins.push_back(empty);
	}
	double start=0;
	double stop=bin_size;
	int j = 0;
	for (int i = 0; i < bin_number; i++){
		while (X[j] < stop and j < X.size()){
			bins[i].push_back(X[j]);
			j+=1;
		}
		start+=bin_size;
		stop +=bin_size;
	}

	start=0;
	stop=bin_size;
	vector <double> edges;
	vector <double> counts;
	double x;
	for (int i =0; i < bin_number; i++){
		if (expand){
			x =(stop + start) / 2.0;
			for (int j = 0; j < bins[i].size(); j++){
				edges.push_back(x);
				counts.push_back(1);
			}
		}
		else{
			edges.push_back((stop + start) / 2.0);
			counts.push_back(bins[i].size());
		}
		start+=bin_size;
		stop +=bin_size;
	}




	return results(edges, counts);
}




#endif /* HISTOGRAM_H_ */
