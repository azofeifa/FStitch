#include "main_segment.h"
#include "read.h"
#include <omp.h>
#include <map>
#include "interval_tree.h"
#include "write.h"
#include "viterbi.h"
#include "validate.h"

#include <string>
using namespace std;
bool isNeg(string OUT){
	for (int i = 0; i < OUT.size(); i++){
		for (int j = i; j < OUT.size(); j++){
			if (OUT.substr(i,OUT.size()-j)=="neg" or OUT.substr(i,OUT.size()-j) =="-"){
				return 1;
			}
		}
	}
	return 0;
}
bool isPos(string OUT){
	for (int i = 0; i < OUT.size(); i++){
		for (int j = i; j < OUT.size(); j++){
			if (OUT.substr(i,OUT.size()-j)=="pos" or OUT.substr(i,OUT.size()-j) =="+"){
				return 1;
			}
		}
	}
	return 0;
}
bool EXIT(string BedGraphFile, string TrainingOutFile, string outFile, string refFile, string strand, string np){
		bool BGF 	= isFile(BedGraphFile);
		bool TOF 	= isFile(TrainingOutFile);
		bool RF 	= (refFile.empty() or isFile(refFile));
		bool NP 	= isNum(np);
		bool ST;
		if (strand == "+" or strand == "-" or strand.empty()){
			ST = 1;
		}else{
			cout<<"-s parameter: "<<strand<<" not understood"<<endl;
			ST = 0;
			return 1;
		}
		if (not BGF){
			cout<<"FILE (-i): "<<BedGraphFile<<", does not exist"<<endl;
			return 1;
		}
		if (not TOF){
			cout<<"Training File (-j): "<< TrainingOutFile<<", does not exist"<<endl;
			return 1;
		}
		if (not RF){
			cout<<"UCSC Table RefSeq File: "<<refFile<<", does not exist"<<endl;
			return 1;
		}
		if (not NP){
			cout<<"Number of Processors(-np): "<<np<<" is not an integer value"<<endl;
			return 1;
		}
		return 0;
}

int run_main_segment(paramsSegment  PT){
	//=================================================================
	// General Parameters
	string BedGraphFile 		= PT.params["-i"];
	string TrainingOutFile 		= PT.params["-j"];
	string outFile 				= PT.params["-o"];
	string refFile 				= PT.params["-r"];
	string strand 				= PT.params["-s"];
	string np 					= PT.params["-np"];
	bool exit_bool 				= EXIT(BedGraphFile, TrainingOutFile, outFile, refFile, strand, np);
	if (exit_bool){
		cout<<"exiting..."<<endl;
		return 0;
	}
	int num_proc 				= stoi(np);
	if (strand.empty()){
		if (isPos(BedGraphFile)){

			strand = "+";
		}else if(isNeg(BedGraphFile)){
			strand = "-";
		}else{
			strand = ".";
		}
	}
	bool verbose 				= not PT.params["-v"].empty();
	//=================================================================
	//Read in FStich Training out file
	if (verbose){
		cout<<"Reading in FStitch Training Out File: ";
	}
	RTOF RTOF_params 					= readTrainingOutFile(TrainingOutFile);	

	if (RTOF_params.EXIT){
		cout<<"exiting..."<<endl;
		return 0;
	}
	if (verbose){
		cout<<"done"<<endl;
		cout<<flush;
	}
	if (verbose){
		cout<<"Training Set from ChIP Data         : "<<(bool(RTOF_params.ChIP)==1)<<endl;	
	}
	//=================================================================
	//Read in BedGraph File by Chromosome
	if (verbose){
		cout<<"Reading in BedGraph File            : ";
		cout<<flush;
	}
	map<string,contig *> ContigData 	= readBedGraphFileAll(BedGraphFile,num_proc);
	if (ContigData.empty()){
		cout<<"exiting..."<<endl;
		return 0;
	}
	
	if (verbose){
		cout<<"done"<<endl;
		cout<<flush;
	}
	//=================================================================
	//Run Viterbi
	if (verbose){
		cout<<"Running Viterbi                     : ";
		cout<<flush;
	}
	map<string, state*> results 	= runViterbi(ContigData, RTOF_params.W, RTOF_params.A,num_proc, RTOF_params.ChIP);
	if (verbose){
		cout<<"done"<<endl;
		cout<<flush;
	}
	//=================================================================
	//Write Viterbi Paths
	if (verbose){
		cout<<"Writing to IGV                      : ";
		cout<<flush;
	}
	writeViterbiPaths(outFile, results,refFile, strand);
	if (verbose){
		cout<<"done"<<endl;
	}
	
	return 1;

	
	
	

	
	
	
	
}