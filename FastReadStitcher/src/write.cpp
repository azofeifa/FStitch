#include <iostream>
#include <fstream>
#include "BaumWelch.h"
#include <time.h>
#include "viterbi.h"
#include "read.h"
#include "interval_tree.h"
#define DTTMFMT "%Y-%m-%d %H:%M:%S "
#define DTTMSZ 21
using namespace std;
static char *getDtTm (char *buff) {
    time_t t = time (0);
    strftime (buff, DTTMSZ, DTTMFMT, localtime (&t));
    return buff;
}
void writeTrainingFile(string OUT,BW_OUT BWO, double alpha, double cm, double ct){
	ofstream FHW;
	FHW.open(OUT);
	if (FHW){	
		char buff[DTTMSZ];
		FHW<<"#####################################################"<<endl;
		FHW<<"#                  Fast Read Stitcher"<<endl;
		FHW<<"#Parameter Estimation Output"<<endl;
		FHW<<"#Date/Time                       :"<<getDtTm(buff)<<endl;
		FHW<<"#Learning Rate                   :"<<to_string(alpha)<<endl;
		FHW<<"#Max Iterations                  :"<<to_string(cm)<<endl;
		FHW<<"#Convergence Threshold           :"<<to_string(ct)<<endl;
		FHW<<"#####################################################"<<endl;
		string weights= "";
		string transitions = "";
		for (int i = 0; i<BWO.W.size(); i++){
			if (i != BWO.W.size()-1){
				weights+=to_string(BWO.W[i]) + ",";
			}else{
				weights+=to_string(BWO.W[i]) ;
			}
		}
		for(int i =0; i < 2; i++){
			for (int j = 0; j < 2; j++){
				if (i!=1 or j!=1){
					transitions+=to_string(BWO.A[i][j])+",";
				}else{
					transitions+=to_string(BWO.A[i][j]);
					
				}
			}
		}
		string converged;	
		if (BWO.converged){
			converged="True";
		}else{
			converged="False";
		}
		FHW<<"Converged                        :"<<converged<<endl;
		FHW<<"Final LogLikelihood              :"<<to_string(BWO.LL)<<endl;
		FHW<<"Logistic Regression Coefficients :"<<weights<<endl;
		FHW<<"HMM Transition Parameter         :"<<transitions<<endl;
		
		


	}else{
		cout<<"couldn't open: "<<OUT<<endl;
		cout<<"exiting..."<<endl;
	}
}
bool isBed(string OUT){
	for (int i = 0; i < OUT.size(); i++){
		for (int j = i; j < OUT.size(); j++){
			if (OUT.substr(i,j)==".bed"){
				return 1;
			}
		}
	}
	return 0;
}

string getGenes(map<string,map<string,T>> DT, string strand, string chrom, int start, int stop){
	T tree 				= DT[strand][chrom];
	vector<interval> F 	= tree.search_interval(start, stop);
	string  genes 		= "";
	map<string, int> Filter;
	for (int i =0; i <F.size(); i ++){
		genes+=F[i].info+",";
	}
	return genes.substr(0, genes.size()-1);
}
map<string,map<string,T>> makeIntervalTree(map<string, map<string, interval *>> R){
	typedef map<string, map<string, interval *>>::iterator strand_it;
	typedef map<string, interval *>::iterator chrom_it;
	map<string,map<string,T>> DT;
	for (strand_it strand = R.begin(); strand!=R.end(); strand++){
		for (chrom_it chrom = strand->second.begin(); chrom!= strand->second.end(); chrom++){
			DT[strand->first][chrom->first] 	= T(chrom->second);
		}
	}
	return DT;
}

void writeViterbiPaths(string OUT, map<string, state*> results, string refFile, string strand){
	map<string, map<string, interval *>> R;
	map<string,map<string,T>> 	DT;
	if (not refFile.empty() and strand != "."){
		R 								= readRefSeq(refFile);
		DT  							= makeIntervalTree(R);

	}
	bool BED 		= isBed(OUT);
	if (not BED){
		OUT+=".bed";
	}
	ofstream FHW;
	FHW.open(OUT);
	char buff[DTTMSZ];
	string score, RGB;	
	//track name=Strand_-description="FStitch" visibility=2 useScore=2 cgGrades=50 cgColour1=white cgColour2=yellow cgColour3=red height=30
	//chr1    10155   10155   ON=0.000000     100     -       10155   10155   255,0,0 N_A...inf
	FHW<<"track name=FStitch_Annotations " <<getDtTm(buff) << "visibility=2 useScore=2 cgGrades=50 cgColour1=white cgColour2=yellow cgColour3=red height=30\n";	
	typedef map<string,state *>::iterator c_it;
	for (c_it chrom = results.begin(); chrom!=results.end(); chrom++){
		state * C 	= chrom->second;
		int prevStart 		= 0;
		int prevStop 		= 0;
		string prevState 	= "";
		double currProb 	= 0;
		double N 			= 0;
		string genes 		= "";
		while(C){
			if (prevState!=C->ID){
				if (not prevState.empty()){
					if (prevState == "ON"){
						score 	= "100";
						if (strand == "+" or strand == "."){
							RGB 	= "0,0,255";
						}else{
							RGB 	= "255,0,0";
						}
					}else{
						score 	= "500";
						RGB 	= "0,255,0";	
					}
					if (not refFile.empty() and strand != "." and prevState=="ON"){
						genes 		= getGenes(DT, strand, C->chrom, prevStart, prevStop);
					}else{
						genes 		= "";	
					}
					FHW<<C->chrom<<"\t"<<to_string(prevStart)<<"\t"<<to_string(prevStop)<<"\t"<<prevState<<"="<<to_string(currProb/N);
					FHW<<"\t"<<score<<"\t"<<strand<<"\t"<<to_string(prevStart)<<"\t"<<to_string(prevStop)<<"\t"<<RGB<<"\t"<<genes<<endl;
					currProb=0;
					N 		=0;
				}
				prevStart = C->start;
				N++;
				currProb+=C->prob;
			}
			prevStop 	= C->stop;
			prevState 	= C->ID;
			C=C->next;

		}
	}
	


}




