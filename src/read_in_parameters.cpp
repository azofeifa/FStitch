#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include "read_in_parameters.h"
using namespace std;
paramWrapper::paramWrapper(){
	train = 0, segment = 0, eRNA = 0, EXIT = 0;	
}
void paramWrapper::display(){
	
	cout<<"============================================================="<<endl
		<<"              Fast Read Stitcher Package                     "<<endl
		<<"Last Update: 1/5/2015"<<endl
		<<"Questions? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
	if (segment){
		cout<<"Making Read Intensity Calls\n"<<endl;
		PS.display();
	}else if (train){
		cout<<"Running Parameter Estimation\n"<<endl;
		PT.display();

	}else if (eRNA){
		cout<<"Making BiDirectional Annotations\n"<<endl;
		PE.display();
	}
	cout<<"============================================================="<<endl;
	
	
}

void paramWrapper::help(){


	string header 	= "";
	header+="--------------------------------------------------------------------------------------\n";
	header+="                          Fast Read Stitcher (FStitch)                        \n";
	printf("%s\n",header.c_str() );
	printf("                   ....description of application modules....             \n");
	printf("                                  (critical)                            \n\n");
	
	printf("train     : must be provided immediately following the application call \"FStitch\"\n");
	printf("              this will compute the MLE estimates for the Logistic Regression\n");
	printf("              and Hidden Markov Model transition matrix. Will output a file that\n");
	printf("              will be used by \"segment\"  module\n");
	printf("segment   : must be provided immediately following the application call \"FStitch\"\n");
	printf("              this will compute the viterbi decoded MLE sequence\n");
	printf("              from provided parameter file from \"train\" module\n");
	printf("              will output two bed files that (forward and reverse strand)\n");
	printf("              for viewing in a genome browser\n");
	
	printf("\n\n");
	header="";
	header+="                ....description of non-default parameters....          \n";
	header+="                                 (critical)                             \n";
	printf("%s\n",header.c_str() );
	printf("-i        : /path/to/bedgraph/file\n");
	printf("              this file should be bedgraph formatted\n");
	printf("              chromosome[tab]start[tab]stop[tab]coverage[newline]\n");	
	printf("-j        : /path/to/training/bed/file\n");
	printf("              Note this file should bed formatted\n");
	printf("              chromosome[tab]start[tab]stop[tab]{0,1}[newline]\n");
	printf("              0 is assumed to be a noise interval\n");
	printf("              1 is assumed to be an actively transcribed interval\n");
	printf("              Note that this bed file should be in reference ONLY to -i\n");
	printf("              in this way only forward strand or reverse strand intervals\n");
	printf("              are considered mutually exclusive\n");
	printf("              the output from \"module\" can be used on the strand it was not trained on\n");
	printf("              as the forward and reverse strand presumably arose from the same experiment\n");
	printf("-k        : /path/to/output/from/train/module\n");
	printf("              this file is used when calling the segment module\n");
	printf("              and specifically uses output from the train module\n");
	printf("-o        : /path/to/output/file\n");
	printf("              In the case of the \"train\" module this will output a file\n");
	printf("              with the learned logistic regression coefficients and tranistion matrix\n");
	printf("              In the case of the \"segment\" module this will output two bed files\n");
	printf("              one for the forward and reverse strand respectively respresenting\n");
	printf("              actively transcibed loci from the provided bedgraph file and parameter\n");
	printf("              estimated file\n\n");

	printf("                    ....description of default parameters....          \n");	
	printf("                                (non-critical)                         \n\n");
	printf("-al       : specific to the \"train\" module\n");
	printf("              controls the learning rate of Newton's Method\n");
	printf("              recomended 0.1<al<0.8 default = 0.4\n");
	printf("-ct       : specific to the \"train\" module\n");
	printf("              controls Newton's Method haulting convergence threshold\n");
	printf("              recomended ct <0.001  default = 0.001\n");
	printf("-cm       : specific to the \"train\" module\n");
	printf("              controls Newton's Method haulting max iterations\n");
	printf("              recomended cm < 500  default = 100\n");
	

	printf("\n");

	printf("\nQuestions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu\n" );
		
	printf("--------------------------------------------------------------------------------------------\n");

	

}



paramsTrain::paramsTrain(){
	params["-i"] 	= "";
	params["-j"] 	= "";
	params["-ij"] 	= "";

	params["-k"] 	= "";
	params["-o"] 	= "";
	params["-r"] 	= "";
	params["-v"] 	= "";
	params["-np"] 	= "8";
	params["-cm"] 	= "100";
	params["-reg"] = "1";
	params["-ct"] 	= "0.1";
	params["-al"] 	= "0.4";
	params["-ms"] 	= "20";
	params["-chip"] = "";
}
void paramsTrain::display(){
	cout<<"BedGraph File  : "<<params["-i"]<<endl
		<<"Training File  : "<<params["-j"]<<endl
		<<"Out File       : "<<params["-o"]<<endl
		<<"Processors     : "<<params["-np"]<<endl;
	if (not params["-r"].empty()){
	cout<<"UCSC Gene Table: "<<params["-r"]<<endl;
		}
	if (not params["-v"].empty()){
	cout<<"verbose output : True"<<endl;
	}
	cout<<endl;
	if (params["-cm"]!="100"){
	cout<<"max iterations (user defined): "<<params["-cm"]<<endl;
	}
	if (params["-ct"]!="0.1"){
	cout<<"convergence threshold (user defined): "<<params["-ct"]<<endl;
	}
	if (params["-al"]!="0.4"){
	cout<<"learning rate (user defined): "<<params["-al"]<<endl;
	}
	if (params["-reg"]!="1"){
	cout<<"regularization parameter (user defined): "<<params["-reg"]<<endl;
	}
	if (not params["-chip"].empty()){
	cout<<"ChIP Data      : True"<<endl;
	}





}
paramsSegment::paramsSegment(){
	params["-i"] 	= "";
	params["-j"] 	= "";
	params["-k"] 	= "";
	params["-o"] 	= "";
	params["-np"] 	= "8";
	params["-r"] 	= "";
	params["-v"] 	= "";
	params["-s"] 	= "";
	params["-chip"] = "";
}
void paramsSegment::display(){
	cout<<"BedGraph File    : "<<params["-i"]<<endl
		<<"Training Out File: "<<params["-j"]<<endl
		<<"Out File         : "<<params["-o"]<<endl;
	if (not params["-r"].empty()){
	cout<<"UCSC Gene Table  : "<<params["-r"]<<endl;
		}
	if (not params["-v"].empty()){
	cout<<"verbose output   : True"<<endl;
	}
}
paramsERNA::paramsERNA(){
	params["-i"] 	= "";
	params["-j"] 	= "";
	params["-k"] 	= "";
	params["-o"] 	= "";
	params["-r"] 	= "";
	params["-v"] 	= "";

}
void paramsERNA::display(){
	cout<<"Forward Strand File: "<<params["-i"]<<endl
		<<"Reverse Strand File: "<<params["-j"]<<endl
		<<"Out File           : "<<params["-o"]<<endl;
	if (not params["-r"].empty()){
	cout<<"UCSC Gene Table    : "<<params["-r"]<<endl;
		}
		if (not params["-v"].empty()){
	cout<<"verbose output     : True"<<endl;
		}
}

void fillInOptions(char* argv[],paramWrapper * P){
	bool segment 	= P->segment;
	bool train 		= P->train;
	bool eRNA 		= P->eRNA;
	string F 		= "";
	char * COM 		= "-";
	while (*argv){
		if ((*argv)[0] == COM[0]){
			F 	= string(*argv); 
		}
		transform(F.begin(), F.end(), F.begin(), ::tolower);
		if (not F.empty()) {
			if (segment && P->PS.params.find(F) !=P->PS.params.end()){
				P->PS.params[F]=string(*argv);
			}else if(train && P->PT.params.find(F) !=P->PT.params.end()){
				P->PT.params[F]=string(*argv);
			}else if (eRNA && P->PE.params.find(F) !=P->PE.params.end()){
				P->PE.params[F]=string(*argv);
			}else{
				cout<<"Unknown user option: "<<F<<endl;
			}
		}
		argv++;
	}

}




void readInParameters( char* argv[], paramWrapper * P){	
	string userModParameter = "";
	argv = ++argv;
	if (not *argv){
		cout<<"please specify either: train or segment"<<endl;
		cout<<"use to -h or --help option for quick reference for software usage"<<endl;
		P->EXIT = 1;
	}else{
		string F 		= *argv;
		if (F.substr(0,2)=="-h" or F.substr(0,6)=="--help"){
			P->help();
			P->EXIT 	= 1;
		}
		if (not P->EXIT){
			F 	= *argv;
			if (F.size()==5 and F.substr(0,5) == "train") {
				P->train = 1;
		 	}
			else if (F.size()==7 and  F.substr(0,7) == "segment") {
				P->segment = 1;
			}
			else if (F.size()==4 and  F.substr(0,4) == "eRNA") {
				P->eRNA = 1;		
			}else{
				cout<<"couldn't understand user provided option: "<<" "<<F<<endl;
				cout<<"use to -h or --help option for quick reference for software usage"<<endl;
				P->EXIT = 1;
			}
			if (not P->EXIT){
				argv = ++argv;
				fillInOptions(argv, P);
			}
		}
	}
}


