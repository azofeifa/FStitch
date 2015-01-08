#include "main_train.h"
#include "read.h"
#include <omp.h>
#include <map>
#include "grabTrainingExamples.h"
#include "interval_tree.h"
#include "NewtonsMethod.h"
#include "BaumWelch.h"
#include "write.h"
#include "validate.h"
using namespace std;
bool EXIT(string BedGraphFile, string TrainingFile, string outFile, string np, string mc, string lr, string ct){
	bool BGF 	= isFile(BedGraphFile);
	bool TOF 	= isFile(TrainingFile);
	bool NP 	= isNum(np);
	bool MC 	= isNum(mc);
	bool LR 	= isNum(lr);
	bool MS 	= isNum(ct);
	if (not BGF){
		cout<<"FILE (-i): "<<BedGraphFile<<", does not exist"<<endl;
		return 1;
	}
	if (not TOF){
		cout<<"Training File (-j): "<< TrainingFile<<", does not exist"<<endl;
		return 1;
	}
	if (not NP){
		cout<<"Number of Processors (-np): "<<np<<" is not an integer value"<<endl;
		return 1;
	}
	if (not MC){
		cout<<"Maximum number of learning iterations (-cm): "<<mc<<" is not a number"<<endl;
		return 1;
	}
	if (not LR){
		cout<<"Learning rate (-lr): "<<lr<<" is not a number"<<endl;
		return 1;
	}
	if (not MS){
		cout<<"Convergence threshold (-ct): "<<ct<<" is not a number"<<endl;
		return 1;
	}

	return 0;

}

int run_main_train(paramsTrain PT){
	bool exit_bool 								= EXIT(PT.params["-i"], PT.params["-j"], PT.params["-o"], PT.params["-np"],
		PT.params["-cm"], PT.params["-lr"], PT.params["-ct"]);
	if (exit_bool){
		cout<<"exiting..."<<endl;
		return 0;
	}

	//=================================================================
	// General Parameters

	string BedGraphFile 						= PT.params["-i"];
	string TrainingFile 						= PT.params["-j"];
	string outFile 								= PT.params["-o"];
	int num_proc 								= stoi(PT.params["-np"]);
	bool verbose 								= not PT.params["-v"].empty();
	//=================================================================
	// parameters specific to Baum Welch and Newtons Method
	double max_convergence 						= stof(PT.params["-cm"]);
	double convergence_threshold 				= stof(PT.params["-ct"]);
	double learning_rate 						= stof(PT.params["-al"]);
	int maxSeed 								= stoi(PT.params["-ms"]);
	
	//=================================================================

	if (verbose){
		cout<<"reading training file                     : ";
	}
	//=================================================================
	//READ TRAINING FILE
	readTrainingFileReturn TrainReturn 	= readTrainingFile(TrainingFile);
	if (TrainReturn.EXIT){
		cout<<"exiting..."<<endl;
		return 0;
	}
	map<string, interval *> TrainingIntervals 	= TrainReturn.result;

	if (verbose){
		cout<<"done"<<endl;
		cout<<"making interval tree                      : ";
	}
	//=================================================================
	//INTERVAL TREE FROM TRAINING FILE
	map<string,T> R 							= makeIntervalTree(TrainingIntervals);
	if (verbose){
		cout<<"done"<<endl;
		cout<<"reading bed graph file                    : ";
		cout<<flush;
	}
	//=================================================================
	//READ BEDGRAPH FILE
	map<string,contig *> ContigData 			= readBedGraphFile(BedGraphFile,TrainingIntervals,1);
	if (ContigData.empty()){
		cout<<"exiting..."<<endl;
		return 0;
	}
	if (verbose){
		cout<<"done\n";
		cout<<"grabbing training data from bedgraph file : ";
		cout<<flush;
	
	}
	//=================================================================
	//GET DATA FROM TRAINING INTERVALS 
	run_out RO  								= run_grabTrainingExamples(R, ContigData);
	if (RO.EXIT){
		cout<<"exiting..."<<endl;
		return 0;
	}
	if (verbose){
		cout<<"done\n";
		cout<<"Begin parameter estimation, Newtons Method: ";
		cout<<flush;
	
	}
	//=================================================================
	//NEWTONS METHOD  
	vector<double> W 							= learn(RO.X, RO.Y, 0, learning_rate);
	if (verbose){
		cout<<"done\n";
		cout<<"Parameter estimation, Baum-Welch          : ";
		cout<<flush;
	}
	//=================================================================
	//BAUM WELCH ALG 
	BW_OUT BWO 									= runBW(ContigData, W,max_convergence, convergence_threshold,learning_rate, verbose, num_proc, maxSeed);
	if (verbose){
		cout<<"done\n";
		cout<<"Writing learned parameters                : ";
		cout<<flush;
	
	}
	writeTrainingFile(outFile, BWO, learning_rate, max_convergence, convergence_threshold);
	if (verbose){
		cout<<"done\n";
	}

	
	
}	