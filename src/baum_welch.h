/*
 * baum_welch.h
 *
 *  Created on: May 2, 2014
 *      Author: joeyazo
 */

#ifndef BAUM_WELCH_H_
#define BAUM_WELCH_H_
#include <unistd.h>
#include <string>
#include <vector>
#include "MEMM_class.h"
#include "readIn.h"
#include "viterbi.h"
#include "functions.h"
#include <stdio.h>
#include "validate_file.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "extract_file_from_path.h"
#include "search_for_files.h"
#include "readInParameters.h"
#include <numeric>
#include "writeOutFile.h"
#include "link.h"
#include <time.h>
#include "example.h"
#include "NewtonsMethod.h"
#include "grabTrainingData.h"
#include "SeperateByStrand.h"
#include "split.h"
#include "MAP_Decoding.h"


class coordinates{
public:
	int start, stop;
	coordinates(int st, int sp){
		start=st;
		stop =sp;
	};
};


vector<MEMM> getBWStruct_array(string Training_FILE, string out_dir){
	string line;
	vector<MEMM> training_struct_array;
	vector<string> lineInfo;
	int start=0;
	int stop =0;
	bool header 	= 1;
	string file="";
	string chrom="";
	string prev_chrom="";
	ifstream TrainingFH(Training_FILE);
	if (TrainingFH.is_open()){
		while ( getline(TrainingFH,line) ){
			if (not header){

				lineInfo 	= split(line, "\t");
				if (lineInfo.size() > 2 and not start and not stop){
					chrom	=	lineInfo[1].c_str();
					start 	= atoi(lineInfo[2].c_str());
					stop  	= atoi(lineInfo[3].c_str());
				}else if( lineInfo.size() > 2 and stop < atoi(lineInfo[3].c_str())){
					stop=atoi(lineInfo[3].c_str());
				}
			}
			else{
				file	= split(line,"\t")[1];
				header=0;
			}
			if (not chrom.empty() and prev_chrom.empty()){
				prev_chrom=chrom;
			}else if(not prev_chrom.empty() and chrom!=prev_chrom){
				break;
			}
			prev_chrom=chrom;


		}
	}else{
		cout<<"Could not open"<<Training_FILE<<endl;
	}

	vector<string> paths ;
	paths.push_back("");
	paths.push_back(out_dir);
	string pathToSrc = "/"+join(split(__FILE__, "/"), split(__FILE__, "/").size()-2,"/");
	paths.push_back(pathToSrc+"segmentation_outputs/");

	string linkedFileName	= "linked_" + split(file, "/")[split(file, "/").size()-1] + ".dat";
	if (not  search(paths, linkedFileName).empty()){
		file	= search(paths, linkedFileName);
		ifstream binaryFile(file, ios::binary|ios::in);
		binaryFile.seekg(0, ios::beg);
		MEMM memmObj;
		bool FOUND	= 0;
		if (binaryFile.is_open()){
			while (1){
				binaryFile.read((char*)&memmObj, sizeof(MEMM));
				if (binaryFile.eof()){
					break;
				}
				if (memmObj.chrom ==chrom and (start < memmObj.begin) and (memmObj.end < stop)){
					training_struct_array.push_back(memmObj);
					FOUND=1;
				}else if(FOUND){
					break;
				}


			}
		}
		else{
			cout<<"Unable to Open Binary File: "<<file<<endl;
		}
	}else{
		cout<<"Couldn't Find: "<<linkedFileName<<endl;
	}
	//now open that bin file and grab coordinates




	return training_struct_array;
}


vector<MEMM> getData(vector<MEMM> struct_array, int start, int stop){
	vector<MEMM> data;
	for (vector<MEMM>::iterator M=struct_array.begin(); M!=struct_array.end(); ++M){
		if (M->begin > stop){
			break;
		}else if (start < M->begin <stop or start< M->end< stop){
			data.push_back(*M);
		}
	}
	return data;

}


vector< vector<double> > runBaumWelch(vector<MEMM> struct_array, vector<logistic_equation> emissions, vector<int> feature_dimension, bool SHOW,int T){
	vector<MEMM> data=struct_array;
	if (SHOW){
		cout<<"...Running Baum Welch to estimate transition parameters..."<<endl;
		cout<<"...Number of training examples for Baum Welch: "<<data.size()<<endl;
	}
	double A	= 0.6;
	int num_of_col = 2;
	int num_of_row = 2;
	vector< vector< double> > gamma(data.size());
	//Initialize the gamma trellis
	for (int i =0; i < data.size(); i++){
		gamma[i]=vector<double>(2);
	}
	int HEIGHT	= data.size();
	int WIDTH	= 2;
	int DEPTH	= 2;
	double ***eta;

	// Allocate memory
	eta = new double**[HEIGHT];
	for (int i = 0; i < HEIGHT; ++i) {
		eta[i] = new double*[WIDTH];
		for (int j = 0; j < WIDTH; ++j)
			eta[i][j] = new double[DEPTH];
	}

	vector< vector<double> > A_matrix(2);
	for (int i =0; i <2; i++){
		A_matrix[i]	= vector<double>(2);
	}


	A_matrix[0][0] = A;
	A_matrix[0][1] = 1.0-A;
	A_matrix[1][0] = 1.0-A;
	A_matrix[1][1] = A;

	int t 	= 0;
	bool converged	= 0;
	double SUM		= 0;
	if (data.empty()){
		cout<<"Error: Did not gather data to estimate baum welch"<<endl;

		return A_matrix;
	}
	else{
		vector< vector<double> > alpha;
		vector< vector<double> > beta;
		vector<double> SUMS;
		int ct=0;
		while (t < T and not converged){

			alpha 	= runForward(data, A_matrix, emissions,feature_dimension, SHOW);
			beta 	= runBackward(data, A_matrix, emissions, feature_dimension, SHOW);
			//Compute normalized gramma trellis
			for (int i = 0; i < alpha.size(); i++){
				SUM	=0.;
				for (int j = 0; j < alpha[0].size(); j++){
					SUM+=alpha[i][j]*beta[i][j];
				}
				SUMS.push_back(SUM);
				for (int j = 0; j < alpha[0].size(); j++){
					gamma[i][j]	= alpha[i][j]*beta[i][j] / SUM;
				}

			}


			//Compute eta
			for (int i = 0; i < gamma.size()-1; i++){
				for (int j = 0; j < gamma[0].size(); j++){
					for (int k = 0; k < gamma[0].size(); k++){
						eta[i][j][k]=(alpha[i][j]*A_matrix[j][k] * beta[i+1][k]* emissions[k].evaluate(data[i].getFeatures(feature_dimension)) )/ SUMS[i];
					}
				}
			}

			//Set new estimates for A_matrix
			double SUM1, SUM2;
			for (int i = 0; i < A_matrix.size(); i++){
				for (int j = 0; j < A_matrix[i].size(); j++){
					SUM1	= 0;
					SUM2	= 0;
					for (int n = 0; n < gamma.size()-1;n++){
						SUM2+=eta[n][i][j];
						SUM1+=gamma[n][i];

					}
					A_matrix[i][j]=SUM2/SUM1;
				}
			}

			//Normalize A_matrix
			for (int i = 0; i < A_matrix.size(); i++){
				SUM1=0.;
				for (int j = 0; j < A_matrix[i].size(); j++){
					SUM1+=A_matrix[i][j];
				}
				for (int j = 0; j < A_matrix[i].size(); j++){
					A_matrix[i][j]/=SUM1;
				}

			}
			if (SHOW){
				cout<<"...Learning..."<<flush;
				if (ct > 3){
					cout<<endl;
					ct=0;
				}
			}
			ct++;

			SUMS.clear();

			t++;

		}
	}
	// De-Allocate memory to prevent memory leak
	for (int i = 0; i < HEIGHT; ++i) {
		for (int j = 0; j < WIDTH; ++j){
			delete eta[i][j];
		}
		delete eta[i];
	}
	delete eta;
	if (SHOW){
		cout<<"\n...Final Transition Matrix Parameters..."<<endl;
		for (int j = 0; j < 2;j ++){
			for (int k =0; k <2;k++){
				cout<<A_matrix[j][k]<<",";
			}
			cout<<endl;
		}

	}

	return A_matrix;

}





#endif /* BAUM_WELCH_H_ */
