#include <iostream>
#include <fstream>
#include "BaumWelch.h"
#include <time.h>
#include "viterbi.h"
#include "read.h"
#define DTTMFMT "%Y-%m-%d %H:%M:%S "
#define DTTMSZ 21
using namespace std;
static char *getDtTm (char *buff) {
    time_t t = time (0);
    strftime (buff, DTTMSZ, DTTMFMT, localtime (&t));
    return buff;
 }
void writeTrainingFile(string OUT,vector<double> W, vector<vector<double>> A,
	double alpha, double cm, double ct){
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
		for (int i = 0; i<W.size(); i++){
			if (i != W.size()-1){
				weights+=to_string(W[i]) + ",";
			}else{
				weights+=to_string(W[i]) ;
			}
		}
		for(int i =0; i < 2; i++){
			for (int j = 0; j < 2; j++){
				if (i!=1 or j!=1){
					transitions+=to_string(A[i][j])+",";
				}else{
					transitions+=to_string(A[i][j]);
					
				}
			}
		}

		FHW<<"Logistic Regression Coefficients :"<<weights<<endl;
		FHW<<"HMM Transition Parameter         :"<<transitions<<endl;
		
		


	}else{
		cout<<"couldn't open: "<<OUT<<endl;
		cout<<"exiting..."<<endl;
	}
}

void writeViterbiPaths(string OUT, map<string, map<string, vector<segment *> >> S){
	
	char buff[DTTMSZ];
	string score, RGB;	

	typedef map<string, map<string, vector<segment *> >>::iterator it_type;
	typedef   map<string, vector<segment *> >::iterator it_type_2;
	for (it_type s = S.begin(); s!=S.end(); s++){
		ofstream FHW;
		if (s->first == "+"){
			FHW.open(OUT+ ".forward.bed");
		}else{
			FHW.open(OUT+ ".reverse.bed");	
		}
	
		FHW<<"track name=FStitch_Annotations " <<getDtTm(buff) << "visibility=1 useScore=2 cgGrades=50 cgColour1=white cgColour2=yellow cgColour3=red height=30\n";	
		for (it_type_2 c= s->second.begin(); c!=s->second.end(); c++){
			for (int i = 0 ; i < c->second.size(); i++){
				int start 	= c->second[i]->start;
				int stop 	= c->second[i]->stop;
				string state 	= "";
				if (c->second[i]->ID == 1){
					score 	= "100";
					state 	= "ON";
					if (s->first=="+"){
						RGB 	= "0,0,255";
					}else{
						RGB 	= "255,0,0";
					}
				}else{
					state 		= "OFF";
					RGB 		= "0,255,0";	
					score 		= "500";
					
				}

				FHW<<c->first<<"\t"<<to_string(start)<<"\t"<<to_string(stop)<<"\t"<<state;
				FHW<<"\t"<<score<<"\t"<<s->first<<"\t"<<to_string(start)<<"\t"<<to_string(stop)<<"\t"<<RGB<<endl;
					
			}
		}
	}
}




