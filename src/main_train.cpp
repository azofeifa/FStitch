#include "main_train.h"
#include "read.h"
#include <omp.h>
#include <map>
#include "grabTrainingExamples.h"
#include "NewtonsMethod.h"
#include "BaumWelch.h"
#include "write.h"
#include "validate.h"
using namespace std;

int run_main_train(paramsTrain PT){

	//=================================================================
	// General Parameters

	string BedGraphFile 						= PT.params["-i"];
	string TrainingFile 						= PT.params["-j"];
	string outFile 							= PT.params["-o"];
	int num_proc 								= stoi(PT.params["-np"]);
	bool verbose 								= not PT.params["-v"].empty();
	bool ChIP 									= not PT.params["-chip"].empty();
	//=================================================================
	// parameters specific to Baum Welch and Newtons Method
	int max_convergence 							= stoi(PT.params["-cm"]);
	double convergence_threshold 				= stof(PT.params["-ct"]);
	double learning_rate 						= stof(PT.params["-al"]);
	double reg 										= -stof(PT.params["-reg"]);
	int maxSeed 									= stoi(PT.params["-ms"]);
	
	//=================================================================

	printf("reading in training file.............................");
	cout.flush();
	//=================================================================
	//READ TRAINING FILE
	vector<segment *> training_segments 	= load::load_intervals_of_interest(TrainingFile);
	if (training_segments.empty()){
		printf("Training Segments Not Populated, is (-j) in the correct format?\n");
		printf("exiting....\n");
		return 0;
	}
	map<string, vector<segment *>> GG 		= load::convert_segment_vector(training_segments);

	printf("done\n");
	printf("inserting bedgraph data to training file intervals...");
	cout.flush();
	vector<segment*> integrated_segments 	= load::insert_bedgraph_to_segment_joint(GG, 
			BedGraphFile);
	printf("done\n");
	printf("transforming data to contig information..............");
	cout.flush();
	run_out RO 								= load::convert_to_run_out(integrated_segments);
	printf("done\n");

	
	//=================================================================
	//NEWTONS METHOD  
	printf("running Newton's Method for Logistic Regression......");
	cout.flush();
	vector<double> W 							= learn(RO.X, RO.Y, 0, 
		learning_rate,max_convergence,reg);
	printf("done\n");

	//=================================================================
	//Markov Optimization  
	printf("estimating transition parameters.....................");
	cout.flush();
	vector<vector<double>> A 			= learn_transition_parameters(W , RO.X,  RO.Y);
	printf("done\n");
	cout.flush();
	
	writeTrainingFile(outFile, W,   A, learning_rate,  max_convergence,  convergence_threshold);
	printf("done\n");
	return 1;
	
	
}	