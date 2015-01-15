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

map<string, state *> runViterbi(map<string,contig *> ContigData, vector<double> W, vector<vector<double>> a, int np, bool ChIP){
	typedef map<string,contig *>::iterator c_it;
	static double pi 	= 0.5;
	static double *  A[2];
	A[0] 	= new double[2], A[1] 	= new double[2];
	A[0][0] = a[0][0],A[1][0] = a[1][0],A[0][1] = a[0][1],A[1][1] = a[1][1];
	map<string, state *> results;
	vector<string> chromosomes;

	for (c_it chrom = ContigData.begin(); chrom!=ContigData.end(); chrom++){
		chromosomes.push_back(chrom->first);
	}	
	int kn 	= chromosomes.size();
	contig * data[kn];
	for (int k=0; k < kn; k++){
		data[k] 	= ContigData[chromosomes[k]];
	}

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(np); // Use 4 threads for all consecutive parallel regions
	#pragma omp parallel for
	for (int k =0; k< kn; k++){
		int T  				= 0;
		contig * C 			= data[k];
		contig * root 		= C;
		string chrom 		= chromosomes[k];
		while(C){
			C=C->next;
			T++;
		}
		C=root;
		//====================================================================
		// Allocate necessary arrays
		double * bj[2];
		bj[0] 	= new double[T], bj[1] 	= new double[T];
		double * beta[2];
		beta[0] 	= new double[T], beta[1] 	= new double[T];
		double * alpha[2];
		alpha[0] 	= new double[T], alpha[1] 	= new double[T];
		double * G[2];
		G[0] 	= new double[T], G[1] = new double[T];
		emissions(root, W, T,bj, ChIP);
		forward(A,bj, T,alpha);
		backward(A, bj, T,beta);
		GAMMA(alpha, beta, T,G);
		int t 				= 0;
		state * trellis[2];
		vector<double> x;
		double  sum, max, emit;
		sum=0, max=0, emit 	= 0;
		state * argmax 		= NULL;
		while(C){
			x 		= C->getVect(ChIP);
			for (int i =0; i <2;i++){
				emit= getEmit(x,W,i);
				if (not  t){
					trellis[i] 				= new state(t, emit*0.5,i,C->start, C->stop, chrom);
				}else{
					sum=-1;
					max=-1;
					for (int j =0; j <2; j++){
						if ((A[j][i]*emit*trellis[j]->max) > max){
							max 	= (A[j][i]*emit*trellis[j]->max);
						}
						if ((A[j][i]*trellis[j]->max) > sum){
							argmax 	= trellis[j];
							sum 	= (A[j][i]*trellis[j]->max);
						}

					}	
					trellis[i]->next 		= new state(t, max,i,C->start, C->stop, chrom);
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
			
			t++;	
			C 	= C->next;
		}

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
		while (R){
			path 		= R;
			path->next 	= prev;
			path->prob 	= G[path->k][t];
			R 			= R->ptr;
			prev 		= path;
			path 	 	= path->prev;
			t--;
		}
		results[chromosomes[k]] 	= prev;

		//====================================================================
		// Deallocate necessary arrays
		for (int i = 0; i < 2; i++){
			delete bj[i], beta[i],alpha[i], G[i];
		}
	}
	return results;
}
































