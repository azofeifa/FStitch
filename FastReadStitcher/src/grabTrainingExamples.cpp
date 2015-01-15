#include <string>
#include <vector>
#include "grabTrainingExamples.h"
#include "read.h"
#include <map>
#include <iostream>
#include "interval_tree.h"
#include <random>
using namespace std;
map<string,T> makeIntervalTree(map<string, interval *> intervals){
	typedef map<string, interval *>::iterator it_type;
	map<string,T> R;	
	for (it_type it = intervals.begin(); it!=intervals.end(); it++){
		R[it->first] 	  = T(it->second);
	}
	return R;
}
run_out::run_out(vector< vector<double> > x, vector<int> y){
	Y=y, X=x;
	EXIT=false;
}
run_out::run_out(){
	EXIT=true;
}
run_out run_grabTrainingExamples(map<string,T> intervals, map<string,contig *> Data, bool ChIP){
	typedef map<string,contig *>::iterator it;
	vector<interval> F;
	int i = 0;
	int k;
	vector<int> Y; 
	vector<vector<double>>X;
						
	for (it data_it = Data.begin(); data_it!=Data.end(); data_it++){
		if (intervals.find(data_it->first) != intervals.end() ){
			contig * C 	= data_it->second;
			T  * tree 	= &intervals[data_it->first];
			while (C!=NULL){
				F 	= tree->search_interval(C->start, C->stop);
				
				if (not F.empty()){
					i++;
					if (F.size()>1){
						for (int i =0; i < F.size(); i++){
							cout<<F[i].start<<"-"<<F[i].stop<<endl;
						}
						cout<<"One or more of your training intervals overlap...\nexiting"<<endl;
						run_out RO;
						return RO;
					}
					else{
						k 	= stoi(F[0].info);
						Y.push_back(k);
						X.push_back(C->getVect(ChIP));
					}
				}
				C 	= C->next;
			}
		}
	}
	//subsample to get save number of ON/OFF examples?
	double one 	= 0;
	double zero = 0;
	for (int i =0 ; i < Y.size(); i++){
		if (Y[i]==1){
			one++;
		}else{
			zero++;
		}
	}
	double N = one+zero;
	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0,1.0);
	double U;
	double ratio;
	if (one > zero){
		ratio 	= one / N;
	}else{
		ratio 	= zero / N;
	}
	vector<vector<double>> newX;
	vector<int> newY;
	double one_one 		= 0;
	double zero_zero  	= 0;
	for (int i =0 ; i < Y.size(); i++){
		U = distribution(generator);
		if (U > ratio && ( (one > zero && Y[i]==1 ) || (one < zero && Y[i]==0 )  )){
			newY.push_back(Y[i]);
			newX.push_back(X[i]);	
			if (Y[i]){
				one_one++;
			}else{
				zero_zero++;
			}
		}else if (( (one > zero && Y[i]==0 ) || (one < zero && Y[i]==1 )  ) ) {
			newY.push_back(Y[i]);
			newX.push_back(X[i]);
			if (Y[i]){
				one_one++;
			}else{
				zero_zero++;
			}
		}
	}
	return run_out(newX,newY);


}