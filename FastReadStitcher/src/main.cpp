//============================================================================
// Name        : main.cpp
// Author      : Joey Azofeifa
// Version     : 1.1
// Description : Main file for running FastReadStitcher Package
//============================================================================
#include <iostream>
#include <omp.h>
#include "read_in_parameters.h"
#include "main_train.h"
#include "main_segment.h"
using namespace std;

int main(int argc, char* argv[]) {
	paramWrapper * P = new paramWrapper();
	readInParameters(argv, P);
	if (P==NULL){
		cout<<"exiting..."<<endl;
		delete P;
		return 0;
	}
	if (P->train){
		paramsTrain PT = P->PT;
		if (not PT.params["-v"].empty()){
			P->display();
		}
	
		run_main_train(PT);
	}else if (P->segment){
		paramsSegment PT = P->PS;
		if (not PT.params["-v"].empty()){
			P->display();
		}
		run_main_segment(PT);
	}else if (P->eRNA){
		cout<<"need to finish segment code"<<endl;
	}


	delete P;
	return 1;
}
