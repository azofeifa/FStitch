// #include <string>
// #include <vector>
// #include "grabTrainingExamples.h"
// #include "read.h"
// #include <map>
// #include <iostream>
// #include <random>
// using namespace std;
// map<string,T> makeIntervalTree(map<string, interval *> intervals){
// 	typedef map<string, interval *>::iterator it_type;
// 	map<string,T> R;	
// 	for (it_type it = intervals.begin(); it!=intervals.end(); it++){
// 		R[it->first] 	  = T(it->second);
// 	}
// 	return R;
// }
// run_out::run_out(vector< vector<double> > x, vector<int> y){
// 	Y=y, X=x;
// 	EXIT=false;
// }
// run_out::run_out(){
// 	EXIT=true;
// }
// run_out run_grabTrainingExamples(map<string,T> intervals, map<string,contig *> Data, bool ChIP){
// 	typedef map<string,contig *>::iterator it;
// 	vector<interval> F;
// 	int i = 0;
// 	int k;
// 	vector<int> Y; 
// 	vector<vector<double>>X;
						
// 	for (it data_it = Data.begin(); data_it!=Data.end(); data_it++){
// 		if (intervals.find(data_it->first) != intervals.end() ){
// 			contig * C 	= data_it->second;
// 			T  * tree 	= &intervals[data_it->first];
// 			while (C!=NULL){
// 				F 	= tree->search_interval(C->start, C->stop);
				
// 				if (not F.empty()){
// 					i++;
// 					if (F.size()>1){
// 						for (int i =0; i < F.size(); i++){
// 							cout<<F[i].start<<"-"<<F[i].stop<<endl;
// 						}
// 						cout<<"One or more of your training intervals overlap...\nexiting"<<endl;
// 						run_out RO;
// 						return RO;
// 					}
// 					else{
// 						k 	= stoi(F[0].info);
// 						Y.push_back(k);
// 						X.push_back(C->getVect(ChIP));
// 					}
// 				}
// 				C 	= C->next;
// 			}
// 		}
// 	}
// 	//subsample to get save number of ON/OFF examples?
// 	double one 	= 0;
// 	double zero = 0;
// 	for (int i =0 ; i < Y.size(); i++){
// 		if (Y[i]==1){
// 			one++;
// 		}else{
// 			zero++;
// 		}
// 	}
// 	printf("%f,%f\n",one,zero );
// 	// for (int i = 0 ; i < Y.size(); i++){
// 	// 	printf("%f,%f,%f,%f\n",Y[i],X[i][0],X[i][1],X[i][2] );
// 	// }
// 	return run_out(X,Y);


// }