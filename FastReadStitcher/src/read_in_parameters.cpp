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
paramsTrain::paramsTrain(){
	params["-i"] 	= "";
	params["-j"] 	= "";
	params["-o"] 	= "";
	params["-r"] 	= "";
	params["-v"] 	= "";
	params["-np"] 	= "8";
	params["-cm"] 	= "100";
	params["-ct"] 	= "0.1";
	params["-al"] 	= "1";
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
	if (params["-cm"]!="100"){
	cout<<"max iterations (user defined): "<<params["-cm"]<<endl;
	}
	if (params["-ct"]!="0.1"){
	cout<<"convergence threshold (user defined): "<<params["-ct"]<<endl;
	}
	if (params["-al"]!="1"){
	cout<<"learning rate (user defined): "<<params["-ct"]<<endl;
	}
	if (not params["-chip"].empty()){
	cout<<"ChIP Data      : True"<<endl;
	}





}
paramsSegment::paramsSegment(){
	params["-i"] 	= "";
	params["-j"] 	= "";
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
		P->EXIT = 1;
	}
	if (not P->EXIT){
		string F 	= *argv;
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
			P->EXIT = 1;
		}
		if (not P->EXIT){
			argv = ++argv;
			fillInOptions(argv, P);
		}
	}
}


