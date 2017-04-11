#include "read.h"
#include "viterbi.h"
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "NewtonsMethod.h"
#include "BaumWelch.h"
#include <omp.h>
using namespace std;
state::state(){
	next=NULL, ptr=NULL, prev=NULL;	
}
state::state(int T, double MAX, int i, int st, int sp, string chr){
	start = st,stop=sp;
	chrom 	= chr;
	t 		= T;
	max 	= MAX;
	if (not i){
		ID 	= "OFF";
	}else{
		ID 	= "ON";
	}
	k=i;
	next=NULL, ptr=NULL, prev=NULL;	

}

double getEmit(vector<double> x, vector<double> W, int k){
	if (not k){
		return 1-g(x , W);
	}else{
		return g(x , W);
	}
}





vector<vector<double>> runViterbi(vector<vector<double>> X, 
				vector<double> W, vector<vector<double>> a ){
	typedef map<string,contig *>::iterator c_it;
	static double pi 	= 0.5;
	static double *  A[2];
	A[0] 	= new double[2], A[1] 	= new double[2];
	A[0][0] = a[0][0],A[1][0] = a[1][0],A[0][1] = a[0][1],A[1][1] = a[1][1];
	vector<vector<double>> results;
	if (not X.size()){
	  return results;
	}
	
	int T  				= X.size();
	//====================================================================
	// Allocate necessary arrays
	
	state * trellis[2];
	vector<double> x;
	double  sum, max, emit;
	sum=0, max=0, emit 	= 0;
	state * argmax 		= NULL;
	string chrom 				= "";
	for (int t = 0 ; t < X.size(); t++){
		x 		= X[t];
		for (int i =0; i <2;i++){
			emit= getEmit(x,W,i);
			if (not  t){
				trellis[i] 				= new state(t, emit*0.5,i,0, 0, chrom);
			}else{
				sum=-1;
				max=-1;
				argmax 	= trellis[0];
				for (int j =0; j <2; j++){
					if ((A[j][i]*emit*trellis[j]->max) > max){
						max 	= (A[j][i]*emit*trellis[j]->max);
					}
					if ((A[j][i]*trellis[j]->max) > sum){
						argmax 	= trellis[j];
						sum 		= (A[j][i]*trellis[j]->max);
					}

				}	
				trellis[i]->next 			= new state(t, max,i,0, 0, chrom);
				trellis[i]->next->ptr 	= argmax;
			}
		}
		if (t){
			for (int i=0;i <2;i++){
				trellis[i] 	= trellis[i]->next;
			}
		}

		sum=0, max=0;
		for (int i =0; i <2; i++){
			max+=(trellis[i]->max);
		}
		for (int i =0 ; i <2; i++){
			trellis[i]->max   /= max;
		}

	}
	int 	t= X.size()-1;
	state * R;
	if (trellis[0]->max > trellis[1]->max){
		R 	= trellis[0];
	}else{
		R 	= trellis[1];			
	}
	state * path;
	state * rootPath 	= path;
	state * prev 		= NULL;
	t--;
	vector <double> current 	= {(double)R->k, R->max};
	results.push_back(current);
	while (R!=NULL and R->ptr!=NULL){
		path 			= R;
		current 		= {(double)R->k, R->max};
		path->next 	= prev;
		R 				= R->ptr;

		results.push_back(current);
		prev 		= path;
		path 	 	= path->prev;
		t--;
	}
	vector<vector<double>> switched;
	for (int t = results.size()-1; t >-1; t--){
		switched.push_back(results[t]);
	}
	//====================================================================
	// Deallocate necessary arrays
	return switched;
}
double get_f1_second(vector<int> Y, vector<vector<double>> pred){
	double tp=0, fp=0,tn=0,fn=0;
	double prob 	=0;
	double P = 0, N = 0;
	for (int i =0 ; i  < Y.size(); i++ ){
		if (pred[i][0]>0 and Y[i] > 0){
			P++;
			tp++;
		}else if(pred[i][0]>0 and Y[i] < 1){
			N++;
			fp++;
		}else if (pred[i][0]<1 and Y[i] < 1){
			N++;
			tn++;
		}else{
			fn++;
			P++;
		}
	}
	double precision 	= tp / P , recall 	= tn / N ;
	double f1 	= 2*precision*recall/(precision+ recall);
	return f1;
}


map<string, map<string, vector<segment *> >> run_viteribi_across(map<string, segment * > G, vector<double> W, vector<vector<double>> A){
	map<string, map<string, vector<segment *> >> S ;
	typedef map<string, segment * >::iterator it_type;
	for (it_type c = G.begin(); c!=G.end(); c++){
		vector<vector<double>> st_sp;
		vector<vector<double>> X 		= (c->second)->get_contig_info(1,st_sp) ;
		vector<vector<double>> results 			= runViterbi(X,W,A);

		double start=0, stop=0,running_sum=0,running_count=0;
		int prev 		= -10;
		string strand 	= ".";
		for (int s = 0 ; s < results.size();s++){
			if (prev!=(int)results[s][0]){
				if (prev!= -10){
					segment * T 	= new segment(c->first, start, stop, prev,strand , running_sum/running_count);
					S[strand][c->first].push_back(T);
				}
				prev 	= (int)results[s][0];
				start = st_sp[s][0];
				running_sum=0,running_count=0;
			}
			stop 		= st_sp[s][1];
			running_count+=1, running_sum+=results[s][1];
		}
	}
	return S;
}




vector<vector<double>> learn_transition_parameters(vector<double> W ,
				vector<vector<double>> X, vector<int> Y){
	int res 	= 100;
	double delta 	= 1.0 / res;
	double a 		= 0;
	vector<vector<double>> final_A;
	double bestf1 			= 0, best_a = 0;
	for (int s = 0 ; s < res; s++)	{
		cout.flush();
		a+=delta;
		vector<vector<double>> A;
		vector<double> ON 	= {1.0 - a, a};
		vector<double> OFF 	= {a, 1.0-a};
		A.push_back(ON);
		A.push_back(OFF);
		vector<vector<double>> results = runViterbi( X,   W, A );
		cout.flush();
		double f1 			= get_f1_second(Y, results);
		cout.flush();
		if (f1 > bestf1){
			best_a 			= a, bestf1 	= f1;
		}
	}
	vector<vector<double>> A;
	vector<double> ONN 	= {1.0 - best_a, best_a};
	vector<double> OFFF 	= {best_a, 1.0-best_a};
	A.push_back(ONN);
	A.push_back(OFFF);
	return A;

}






























