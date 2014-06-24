/*
 * MEMM_class.h
 *
 *  Created on: Jan 22, 2014
 *      Author: joeyazo
 */
#include <vector>
#include <stdlib.h>
#include <string>
#include <cstring>
#ifndef MEMM_CLASS_H_
#define MEMM_CLASS_H_
using namespace std;
class MEMM{

public:
	int begin;
	int end;
	int length;
	int SegCall;
	char strand[10];
	bool GENE;
	char annotation[150];
	char chrom[10];
	double median;
	double mean;
	double mode;
	double max;
	double min;
	double variance;
	double cov_number 	=	0;
	double density  	= 	0.0;
	double score 		= 	0.;


	MEMM(){};
	vector<double> getFeatures(vector<int> code){
		vector<double> X;
		X.push_back(1);
		for(vector<int>::iterator c = code.begin(); c!=code.end(); ++c){
			if (*c==1){
				X.push_back(length);
			}else if(*c==2){
				X.push_back(mean);

			}else if(*c==3){
				X.push_back(mode);

			}else if(*c==4){
				X.push_back(max);

			}else if(*c==5){
				X.push_back(min);

			}else if(*c==6){
				X.push_back(cov_number);

			}else if(*c==7){
				X.push_back(density);

			}else if(*c==8){
				X.push_back(variance);

			}else if(*c==9){
				X.push_back(GENE);

			}
		}
		return X;

	}
	MEMM(string input_strand , string begin_chrom, int val1, int val2,  string gene_ID, int length_val, double coverages_int , double max_cov_dbl, double variance_coverage){
		begin 		= val1;
		end 		= val2;
		length 		= length_val;
		max 		= max_cov_dbl;
		mean		= coverages_int;
		SegCall		= 0;
		variance	= variance_coverage;
		strcpy(strand,input_strand.c_str());
		strcpy(chrom,begin_chrom.c_str());

		strcpy(annotation,gene_ID.c_str());
		if (gene_ID=="N_A"){
			GENE = 0;
		}
		else{
			GENE = 1;
		}
	}

	void setDensity(double v1){
		density 	= v1;
	}

	void show(){
		cout<<chrom<<"\t"<<begin<<"\t"<<end<<"\t"<<"\t"<<length<<"\t"<<annotation<<"\t"<<mean<<"\t"<<density<<"\t"<<endl;
	}
	vector<double> getPosCov(){
		vector<double> pos_cov(2);
		pos_cov[1]	= cov_number;
		pos_cov[0]	= begin;
		return pos_cov;
	}
	void setStats(double v1, double v2, double v3, double v4, double v5, double v6,double v7, double v8){
		mean		= v1;
		median 		= v2;
		mode		= v3;
		max			= v4;
		variance	= v5;
		cov_number	= v6;
		density 	= v7;
		min			= v8;

	}


	void set_values(int, int, int,double,bool, string, string, string);
	void turnOn();
	vector<double> get_features(int);
};

void MEMM::set_values(int val1, int val2, int val3, double val4, bool val5, string val6, string val7, string val8){
	begin = val1;
	end   = val2;
	length = val3;
	mean = val4;
	GENE = val5;
	strcpy(annotation, val6.c_str());
	strcpy(strand,val7.c_str());
	strcpy(chrom,val8.c_str());

}




#endif /* MEMM_CLASS_H_ */
