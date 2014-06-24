/*
 * NewtonsMethod.h
 *
 *  Created on: Feb 18, 2014
 *      Author: Joey Azofeifa
 */

#ifndef NEWTONSMETHOD_H_
#define NEWTONSMETHOD_H_
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <numeric>
#include "split.h"
#include <unistd.h>
#include <string>
#include <cmath>
#include <math.h>
using namespace std;


vector< vector<double> > transpose(vector<vector<double>> H){
	//==============================================================
	//	!!!!!!!!!!!!!!!This function is incomplete!!!!!!!!!!!!!!!
	//	Only handles square matrices and vectors...
	//==============================================================
	int dim = H.size();
	int k,l, temp;
	if (dim ==1 ){ //this is just a row vector...make a column vector
		vector<vector<double>> newH;
		for (int i = 0; i < H[0].size(); i ++){
			newH.push_back(vector<double>(1));
			newH[i][0] 	= H[0][i];
		}
		return newH;
	}
	if (H[0].size() == 1 ){//this is just a column vector...make a row vector
		vector<vector<double>> newH;
		newH.push_back(vector<double>(H.size()));
		for (int i =0; i < H.size(); i++){
			newH[0][i]=H[i][0];
		}
		return newH;
	}

	//This should only happend below if its a square matrix....
	for (int r = 0; r < dim; r++){
		for (int i =r; i < dim; i++){
			for (int j = i; j < dim; j++){
				if (i !=j){
					k 		= j;
					l		= i;
					temp 	= H[i][j];
					H[i][j] = H[k][l];
					H[k][l] = temp;
				}
			}
			break;
		}
	}
	return H;
}
double g(vector<double> x, vector<double> w){
	if (x.size() != w.size()){
		cout<<"!!weight and x vectors are not same dimension!!"<<endl;
		return 0;
	}
	double sum = 0;
	for (int i = 0; i < x.size(); i++){
		sum+=x[i]*w[i];
	}
	return exp(sum) / (1.0 + exp(sum));
}
vector<vector<double> > gradient(vector< vector<double> > X, vector<int> Y,vector<double> W){
	vector<vector<double>> gradient;
	int N	= X.size(); //Number of Training Examples
	int K	= W.size(); //Feature Dimensions
	int i, j;
	gradient.push_back(vector<double>(K));
	for (i = 0; i <K; i++){
		gradient[0][i]=0;
	}


	for (i = 0; i < N; i++ ){
		for (j = 0; j < K; j++ ){
			gradient[0][j]+=X[i][j]*(Y[i] - g(X[i], W));
		}
	}
	return gradient;
}
vector< vector<double> > hessian(vector< vector<double> > X, vector<int> Y, vector<double> W){
	vector< vector<double> > H;
	int N	= X.size(); //Number of Training Examples
	int K	= W.size(); //Feature Dimensions
	int i, j, n;
	//Initialize H
	for (j = 0; j < K; j++){
		H.push_back(vector<double>(K));
		for (n = 0; n < K; n++){
			H[j][n] = 0;
		}
	}

	for (i = 0; i<N;i++){
		for (j = 0; j < K; j++){
			for (n = 0; n < K; n++){
				H[j][n]+=X[i][j]*X[i][n]*g(X[i], W)*(1-g(X[i], W));
			}
		}
	}
	for (j = 0; j < K; j++){
		for (n = 0; n < K; n++){
			H[j][n]=(-1*H[j][n]);
		}
	}
//	cout<<"Not Inverted Hessian"<<endl;
//	for (int h =0; h < K; h++ ){
//		for (int o =0; o < K; o++ ){
//			cout<<H[h][o]<<"\t";
//		}
//		cout<<endl;
//	}


	return H;
}
vector<double> mult(vector<double> A, double X){
	vector<double> a(A.size());
	for(int i = 0; i < A.size(); i++){
		a[i] = A[i] *X;
	}
	return a;
}
vector<double> sub(vector<double> A, vector<double> B){
	vector<double> a(A.size());
	for(int i = 0; i < A.size(); i++){
		a[i] = A[i] - B[i];
	}
	return a;
}
vector<double> add(vector<double> A, vector<double> B){
	vector<double> a(A.size());
	for(int i = 0; i < A.size(); i++){
		a[i] = A[i] + B[i];
	}
	return a;
}
vector< vector<double> > identity(int dim){
	vector< vector<double> > I(dim);
	for (int i =0; i< dim; i++){
		I[i]	= vector<double>(dim);
		I[i][i] = 1;
	}
	return I;
}
vector<vector<double> > add_matrices(vector<vector<double>> A, vector<vector<double>> B ){
	vector<vector<double> > C 	= identity(A.size());
	for (int i =0; i <A.size(); i++){
		for (int j =0; j <A[i].size();j++){
			C[i][j] 	= A[i][j] + B[i][j];
		}
	}
	return C;


}


double dot(vector<double> A, vector<double> B){
	int dim		= A.size();
	double sum	= 0;
	if (A.size() != B.size()){
		cout<<"Vectors need to be the same dimension..."<<endl;
		return 0;
	}
	for (int i = 0; i < dim; i++){
		sum+=(A[i]*B[i]);
	}
	return sum;
}

vector< vector<double> > concatenate(vector< vector<double> > A, vector< vector<double> > B){
	int dim1	= A.size();
	int dim2	= B.size();

	for (int i=0;i<dim1;i++){
		for (int j =0;j<dim2;j++){
			A[i].push_back(B[i][j]);
		}
	}
	return A;

}
vector<vector<double>> matrix_mult(vector<vector<double>> A, vector<vector<double>> B){
	int A_i 	= A.size();
	int A_j 	= A[0].size();
	//dimensions of matrix A is A_i x A_j
	int B_i 	= B.size();
	int B_j 	= B[0].size();
	//dimensions of matrix B is B_i x B_j
	//check to make sure that A_j == B_i
	if (A_j != B_i){
		cout<<"Matrices need to be in the right dimensions..."<<endl;
	}
	//The new matrix becomes dimensions A_i j B_j
	vector< vector<double> > H(A_i);
	for (int i =0; i< A_i; i ++){
		H[i] = vector<double>(B_j);
	}
	//take the tranpose of B...it will be easier to the dot product of the column vectors
	B 	= transpose(B);
	//Now fil out the new matrix by doing a series of dot products
	for (int i = 0; i < A_i; i++){
		H[i]=vector<double>(B_j);
		for (int j = 0; j< B_j; j++){
			H[i][j] = dot(A[i], B[j]);
		}
	}
	return H;
}
vector< vector<double> > inv(vector<vector<double>> H){ //uses gauss-jordon elimination

	//=======================================================================================================
	//	This was kind of a bitch to write, but the matrix data structure is a vector of vectors
	//	Need to put in some checks like make sure its a square matrix and has a non-zero determinat
	//  But yeah it works like every other algorithm out there, place next to the matrix of interest its
	// 	identity matrix and do row eliminations on the matrix of interest until you get it down to identity
	//	then just return the now modified identity matrix.
	//=======================================================================================================

	int dim = H.size();
	int i	= 0;
	int j	= 0;
	int k	= 0;
	int l	= 0;

	vector< vector<double>> I = identity(dim);
	vector< vector<double>> H_T = transpose(H);
	H = concatenate(H,I);
	double multiple;
	vector<double> row1,row2, erow;


	while (i < dim){
		if (H[i][i] == 0){
			k=j;
			while (k<dim){ 			//find a non zero entry in the ith column and swap
				if (H[k][i] != 0){  //if all non zero than just increase i,j by one
					break;
				}
				k++;
			}
			if (k == dim){ //just need to increase i,j by one...everything below was zeros
				H[i][i]++;
			}
			else{ //here is the swapping...
				row1 	= H[i];
				row2 	= H[k];
				H[i] 	= row2;
				H[k]	= row1;
			}
		}
		//now we need to do the row reduction to get it in reduced row echelon form
		//first make i,j equal to one by deviding the row by i,j
		H[i] 	= mult(H[i], (1.0 / H[i][i]));
		//now eliminate all other entries in the jth column by subtracking suitable
		//multiples of the ith row from other rows
		j=i;
		j+=1;
		while (j < dim){
			if (H[j][i] != 0){
				multiple	= H[i][i] * H[j][i];
				row1 		= mult(H[i], (multiple / H[i][i]));
				row2 		= mult(H[j], (multiple / H[j][i]));
				if ((row1[i] + row2[i]) == 0){
					erow 	= add(row1,row2);
				}
				else if ((row1[i] - row2[i]) == 0){
					erow 	= sub(row2,row1);
				}
				else{
					cout<<"ERROR Forward Pass"<<endl;
//					for (int h =0; h < dim; h++ ){
//						for (int o =0; o < dim; o++ ){
//							cout<<H[h][o]<<"\t";
//						}
//						cout<<endl;
//					}

				}
				H[j] 	= erow;
			}
			j+=1;
		}
		i+=1;
		j=i;
	}
	//So now that its in reduced echelon form, we have to go back and make the other side of triangle all zeros
	i-=1;
	while (i>-1){
		j=i-1;

		while (j > -1){
			if (H[j][i] != 0){
				multiple 	= H[i][i]*H[j][i];
				row1 		= mult(H[i], (multiple / H[i][i]));
				row2 		= mult(H[j], (multiple / H[j][i]));

				if ((row1[i] + row2[i]) == 0){
					erow 	= add(row1,row2);
				}else if((row1[i] - row2[i]) == 0){
					erow 	= sub(row2,row1);
				}
				else{
					cout<<"WHAT(Backward Pass)"<<endl;
				}
				H[j] 	= erow;
			}
			j-=1;
		}
		i-=1;
	}

	//Now Just return the inverse

	vector< vector<double>> INV;
	for (i = 0; i < H.size();i++){
		INV.push_back(vector<double>(dim));
		for (j = dim; j < H[i].size();j++){
			INV[i][j-dim] 	= H[i][j];
		}
	}

	return INV;
}
double loglikelihood(vector<vector<double>>X,  vector<int> Y, vector<double> W){
	double ll = 0.0;

	for (int i =0; i < X.size(); i++){
		ll+=(	(log(1-g(X[i], W) )) + (Y[i]*(dot(X[i], W)))	);
	}
	return ll;
}
double cost(vector<vector<double>>X,  vector<int> Y, vector<double> W){
	double COST		= 0.0;
	int prediction 	= 0;
	for (int i =0; i < X.size(); i++){
		if ((g(X[i], W)) < 0.5){
			prediction = 0;
		}
		else{
			prediction = 1;
		}

		COST+=(abs(prediction-Y[i]));
	}
	COST 	= (COST / X.size());
	return COST;
}

double difference(vector<double> W, vector<double> P_W){
	double sum=0;
	for (int i =0; i < W.size(); i++){
		sum+=abs(W[i]-P_W[i]);
	}
	return sum;
}


vector<double> NewtonsMethod(vector< vector<double> > X, vector<int> Y, bool verbose, double alpha){
	vector<double> W;
	for (int i = 0; i < X[0].size(); i++){
		W.push_back(0.0); 	//Initialize Starting Weights to Zero
	}
	vector< vector<double>> reg = identity(X[0].size());



	//========================================================================
	//	Initialize some important data structures...basically a bunch of
	//	vectors of vectors to replicate matrices and vectors
	//	The below loop runs for a fixed number of iterations or stops when the
	//	loglikelihood converges
	//========================================================================

	vector<vector<double>> H;
	vector<vector<double>> INV_H;
	vector<vector<double>> G;
	vector<vector<double>> current_gradient;

	int T							= 200;
	int t							= 0;
	double LL 						= 0;
	double previous					= 0;
	double *ptr 					= &previous;
	double convergence_threshold 	= 1.0 / 1000000;
	double prev_T					= 0;
	vector<double> P_W(W.size());
	vector<double> gamma(W.size());
	int ct=0;
	while (t < T){

		G 					= gradient(X,Y,W);

		H 					= hessian(X,Y,W);
		for (int h =0; h < H.size(); h++ ){
			for (int o =0; o < H.size(); o++ ){
				if (H[h][o]!=H[h][o]){
					cout<<"Logistic Regression Blew Up...Did not converge..."<<endl;
					cout<<"Final Parameter Set...."<<endl;
					for (int i=0; i < W.size(); i++){
						cout<<to_string(W[i]).substr(0,10)<<"\t";
					}
					cout<<endl;
					return W;
				}
			}
		}
		INV_H 				= inv(H);
		current_gradient	= transpose(matrix_mult(INV_H, transpose(G)));
		gamma				= sub(W,current_gradient[0]);
		W 					= add(mult(W, 1-alpha), mult(gamma, alpha) );

		LL 					= loglikelihood(X,Y,W);
		if (verbose){
			cout<<"...Learning..."<<flush;
			if (ct > 3){
				cout<<endl;
				ct=0;
			}

//			for (int i=0; i < W.size(); i++){
//				cout<<to_string(W[i]).substr(0,10)<<"\t";
//			}
//			cout<<endl;
		}
		ct++;
		t++;
		if(ptr != NULL && abs(difference(W,P_W)) < convergence_threshold){
			if (verbose){



				cout<<"\n...Newton's Method converged in "<<t<<" iterations"<<endl;
				cout<<"...Learned Logistic Regression Parameters: ";
				for (int i=0; i < W.size(); i++){
					cout<<to_string(W[i]).substr(0,10)<<"\t";
				}

				cout<<endl;
			}
			return W;
		}
		*ptr = LL;
		P_W 	= W;
		if (( t / (double)T ) > (prev_T+0.25) and verbose){
			prev_T 	= t / (double)T;
		}
	}
	cout<<"Warning: Learning Algorithm did not converge..."<<endl;
	cout<<"Consider increasing number of iterations or lowering convergence threshold..."<<endl;
	cout<<"Learned Logistic Regression Parameters: ";
	for (int i=0; i < W.size(); i++){
		cout<<W[i]<<"\t";
	}
	cout<<endl;

	return W;


}








vector<double> learn(vector< vector<double> > X, vector<int> Y, bool verbose, double alpha){

	vector<double> W 	= NewtonsMethod(X,Y, verbose,alpha);

	return W;
}




int readInTestNewton(){
	string fileName = "/Users/joeyazo/Desktop/MachineLearning/ML_assign_3_azofeifa_joey/src/ex2data1.txt";
	string line;
	vector<string> lineInfo;
	ifstream myfile(fileName);
	vector< vector<double> > X;
	vector<int> Y;
	vector<double> W;



	int t = 0;
	if (myfile.is_open()){
		while ( getline (myfile,line) ){
			lineInfo = split(line, ",");
			X.push_back(vector<double>(lineInfo.size()) );
			for (int i = 0; i < lineInfo.size(); i++){
				if (i == 0){
					X[t][i] = 1;
				}
				else{
					X[t][i] = atof(lineInfo[i-1].c_str());
				}
			}
			Y.push_back(atoi(lineInfo[lineInfo.size()-1].c_str()));
			t+=1;
		}
	}
	else{
		cout<<"Couldn't Open: "<< fileName <<endl;
		return 0;
	}
	NewtonsMethod(X,Y,1,1);

	return 1;
}



























#endif /* NEWTONSMETHOD_H_ */
