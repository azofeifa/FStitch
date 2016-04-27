#include "BaumWelch.h"
#include "read.h"
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <iostream>
#include <random>
#include "NewtonsMethod.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h> 
#include <omp.h>
using namespace std;
class b{
public:
	int t;
	float p;
	b * next;
	b(int T,float P){
		t=T;
		p=P;
	} 
};
void emissions(contig * dd, vector<double>W, int T, double ** bj, bool ChIP){
	contig * d = dd;
	for (int t =0; t< T; t++){
		for (int i =0;i <2; i++){
			if (i==1){
				bj[i][t] 	= g(d->getVect(ChIP), W);
			}else{
				bj[i][t] 	= 1.0 - g(d->getVect(ChIP), W);
			}
		}
		d=d->next;
	}
	
}
void forward(double ** A, double ** bj, int T, double **  alpha){
	int t 		= 0;
	double S;
	for (int t =0; t< T; t++){
		if (t==0){
			alpha[0][t] 	= bj[0][t]*0.5;
			alpha[1][t] 	= bj[1][t]*0.5;
		}else{
			for (int i =0; i <2; i++){
				S=0;
				for (int j = 0; j< 2; j++){
					S+=alpha[j][t-1] * A[j][i];
				}
				alpha[i][t]=(S*bj[i][t]);
			}
			//normalize, help with underflow problems
			S=(alpha[0][t]+alpha[1][t]);
			alpha[0][t]=(alpha[0][t]/S),alpha[1][t]=(alpha[1][t]/S);
			
		}
	}
}
void backward(double ** A,  double ** bj, int T,double ** beta){
	int t 		= T;
	double S;
	for (int t =T-1; t>= 0; t--){
		if (t==T-1){
			beta[0][t] 	= 1.0;
			beta[1][t] 	= 1.0;
		}else{
			for (int i =0; i <2; i++){
				S=0;
				for (int j = 0; j< 2; j++){
					S+=(beta[j][t+1] * A[i][j]*bj[j][t+1]);
				}
				beta[i][t]=S;
			}
			//normalize, help with underflow problems
			S=(beta[0][t]+beta[1][t]);
			beta[0][t]=(beta[0][t]/S),beta[1][t]=(beta[1][t]/S);
			
		}
	}
}
void GAMMA(double **alpha, double **beta, int T,double **G){
	double S;
	for (int t=0; t<T; t++){
		for (int i=0; i<2;i++){
			S 	= 0;
			for (int j=0;j<2;j++){
				S+=(alpha[j][t]*beta[j][t]);
			}
			G[i][t] 	= (alpha[i][t]*beta[i][t])/S;
		}
	}	
}
void EPSILON(double **alpha, double **beta, double ** bj, double ** A, int T, double *** E){
	double S;
	for (int t=0; t < T-1; t++){
		S 	= 0;
		for (int i = 0; i < 2;i++){
			for (int j=0; j<2;j++){
				S+=(alpha[i][t]*A[i][j]*bj[j][t+1]*beta[j][t+1]);
			}
		}
		for (int i = 0; i < 2;i++){
			for (int j=0; j<2;j++){
				E[i][j][t] 	= (alpha[i][t]*A[i][j]*bj[j][t+1]*beta[j][t+1]) / S;
			}
		}
	}
}
double sumArray(double * array, int T){
	double S = 0;
	for (int t = 0; t<T;t++){
		S+=array[t];
	}
	return S;
}
vector<int> decode(double ** G, int T){
	vector<int> newY;
	for (int t=0; t< T; t++){
		if (G[0][t] > G[1][t] and G[0][t]>0.6){
			newY.push_back(0);
		}else if(G[1][t]>0.6){
			newY.push_back(1);
		}else{
			newY.push_back(2);
		}
	}
	return newY;
}
class subSampleOut{
public:
	vector<int> YY;
	vector<vector<double>> XX;
	subSampleOut(vector<vector<double>> XXX, vector<int> YYY){
		YY=YYY;
		XX=XXX;
	}
	subSampleOut(){};
};
subSampleOut subsample(vector<int> YY, vector<vector<double>> XX, int MAX){
	int zeroCT=0;
	int oneCT=0;
	srand (time(NULL));
	default_random_engine generator(time(NULL));
	uniform_real_distribution<double> distribution(0.0,1.0);
	vector<int> YYY;
	vector<vector<double>> XXX;
	int t=0;
	double U;
	U = distribution(generator);
	while ((zeroCT < MAX or oneCT < MAX) and t < YY.size()){
		U = distribution(generator);
	
		if (zeroCT < MAX and  U < 0.1 and YY[t]==0){
			YYY.push_back(YY[t]);
			XXX.push_back(XX[t]);
			zeroCT++;
		}
		if (oneCT < MAX and  U < 0.1 and YY[t]==1){
			YYY.push_back(YY[t]);
			XXX.push_back(XX[t]);
			oneCT++;
		}
		t++;
	}
	return subSampleOut(XXX,YYY);
}
double getLL(double ** G, int T){
	double LL=0;
	for (int t=0; t<T;t++){
		if (G[0][t] > G[1][t]){
			LL+=log(G[0][t]);
		}else{
			LL+=log(G[1][t]);	
		}
	}
	return LL;
}
BW_OUT::BW_OUT(bool c, double ll, vector<double> w, double ** a){
	converged 	= c;
	LL 			= ll;
	W 			= w;
	vector<vector<double>>aa;
	for (int i =0; i < 2; i ++){
		vector<double> aaa;
		for (int j =0; j < 2; j ++){
			aaa.push_back(a[i][j]);
		}	
		delete a[i];
		aa.push_back(aaa);
	}
	A 			= aa;
}
BW_OUT::BW_OUT(){}
BW_OUT runBW(map<string,contig *> D, vector<double> W, double cm, double ct, double learning_rate, bool verbose, int np, int maxSeed, bool ChIP){
	vector<BW_OUT> ALL(maxSeed);
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(np); // Use 4 threads for all consecutive parallel regions
	#pragma omp parallel for
	for (int i = 0; i < maxSeed; i++){

		string chrom 		= D.begin()->first;
		contig * X 			= D[chrom];
		contig * root 		= X;
		int t 				= 0;
		vector<vector<double>> dataX;
		//====================================================================
		// Get length of linked list and put into vector<vector<double>>  
		while (X){
			dataX.push_back(X->getVect(ChIP));
			X 				= X->next;
			t++;
		}
		int T = t;
		//====================================================================
		// Allocate necessary arrays
		double * A[2];
		A[0] 	= new double[2], A[1] 	= new double[2];
		A[0][0] = 0.8,A[1][0] = 0.2,A[0][1] = 0.2,A[1][1] = 0.8;
		double * bj[2];
		bj[0] 	= new double[T], bj[1] 	= new double[T];
		double * beta[2];
		beta[0] 	= new double[T], beta[1] 	= new double[T];
		double * alpha[2];
		alpha[0] 	= new double[T], alpha[1] 	= new double[T];
		double * G[2];
		G[0] 	= new double[T], G[1] = new double[T];
		double ** E[2];
		for (int i =0; i <2;i++){
			E[i] 	= new double*[2];
			for (int j = 0;j<2;j++){
				E[i][j] 	= new double[T-1];
			}
		}
		//====================================================================

		int k = 0;
		vector<double> newY;
		vector<vector<double>> newX;
		bool begin=0;
		double prevLL 	= 0.;
		double LL 		= 1000;
		bool converged 	= 0;
		
		while (not converged and k < cm){
			emissions(root, W, t,bj, ChIP);
			forward(A,bj, t,alpha);
			backward(A, bj, t,beta);
			GAMMA(alpha, beta, t,G);
			EPSILON(alpha, beta, bj, A, t,E);
			// if (begin){
			// 	subSampleOut SS 	= subsample(decode( G, t), dataX, 2000);
			// 	//this is to get a broader more diverse range of values to learn on
			// 	W 					= learn(SS.XX, SS.YY, 0, learning_rate); //learned logistic regression parameters
			// 	begin 				= 0;
			// }
			double NN,DD;
			for (int i =0; i <2; i++){
				NN=0, DD=sumArray(G[i],t-1);
				for (int j = 0; j < 2; j++){
					NN 	= sumArray(E[i][j],t-1);
					A[i][j] 	= NN/DD;
				}
			}
			LL 	= getLL(G, t);
			if (abs(prevLL - LL) < ct){
				converged 	= 1;
			}
			k++;
			prevLL 	= LL;
		}
		//====================================================================
		// Deallocate necessary arrays
		for (int i = 0; i < 2; i++){
			delete bj[i], beta[i],alpha[i], G[i];
			for (int j = 0; j<2;j++){
				delete E[i][j];
			}
			delete E[i];	
		}
		
		//====================================================================
		BW_OUT BWO 		= BW_OUT(converged, LL, W,A);
		ALL[i] 				= BWO;
	}
	BW_OUT maxBWO 			= ALL[0];
	for (int i=1; i < ALL.size(); i ++){
		if (ALL[i].LL > maxBWO.LL){
			maxBWO 			= ALL[i];
		}
	}

	return maxBWO;
	
}





















