#include "main_segment.h"
#include "read.h"
#include <omp.h>
#include <map>
#include "write.h"
#include "viterbi.h"
#include "validate.h"

#include <string>
using namespace std;

int run_main_segment(paramsSegment  PT){
	// =================================================================
	// General Parameters
	string BedGraphFile 			= PT.params["-i"];
	string TrainingOutFile 		= PT.params["-k"];
	string outFile 				= PT.params["-o"];
	string refFile 				= PT.params["-r"];
	string np 					= PT.params["-np"];
	printf("loading training/parameter out file................");
	cout.flush();
	RTOF RTOF_params 					= readTrainingOutFile(TrainingOutFile);	
	printf("done\n");
	printf("loading bedgraph file(s)...........................");
	cout.flush();
	map<string, segment*> 			G 	= load::load_bedgraphs_total(BedGraphFile);
	printf("done\n");
	if (G.empty()){
		printf("Bedgraph file(s) were not loaded...\n");
		printf("Exiting...\n");
		return 0;
	}

	printf("running viterbi....................................");
	cout.flush();
	map<string, map<string, vector<segment *> >> S  = run_viteribi_across(G , RTOF_params.W , RTOF_params.A );
	printf("done\n");
	printf("writing viterbi....................................");
	cout.flush();
	writeViterbiPaths(outFile, S);
	printf("done\n");
	

	
	
	
	return 1;

	
	
	

	
	
	
	
}
