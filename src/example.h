/*
 * example.h
 *
 *  Created on: Feb 12, 2014
 *      Author: joeyazo
 */

#ifndef EXAMPLE_H_
#define EXAMPLE_H_
#include "readInRefGene.h"

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
#include "baum_welch.h"
#include "grabTrainingData.h"
bool unitTest(string pathToSrc){

	//========================================================
	//	This serves as a unit test library to make sure
	//	everything will run okay, we are going to use a truncated
	//  bedGraph file which will run quickly, binaries will be made
	//	learning algorithms benchmarked etc.
	//	Indeed, this also serves as a perfect run, so please
	//	make sure your files are formatted like the test set.
	//========================================================

	string BedGraphFile		= pathToSrc+"examples/TEST_DMSO2_3.pos.BedGraph";
	string refGeneFileName 	= pathToSrc +"refSeqAnnotations/refSeqGene_hg19.txt";
	string Training_File	= pathToSrc+"training_examples/testTrainingFile.tsv";
	string strand			= "+";

	if (validate(pathToSrc+"src/main.cpp").empty()){
		cout<<"Error Code(1)"<<endl; //Path to src is not set properly
		cout<<"Exiting..."<<endl;
		return 0;

	}
	if (validate(pathToSrc+"examples/TEST_DMSO2_3.pos.BedGraph").empty()){
		cout<<"Error Code(2)"<<endl; //Path to src is not set properly
		cout<<"Exiting..."<<endl;
		return 0;
	}
	if (validate(Training_File).empty()){
		cout<<"Error Code(3)"<<endl; //Path to src is not set properly
		cout<<"Exiting..."<<endl;
		return 0;
	}
	string out_dir = pathToSrc + "segmentation_outputs/";
	if (validate(out_dir).empty()){
		cout<<"Error Code(4)"<<endl; //Path to src is not set properly
		cout<<"Exiting..."<<endl;
		return 0;
	}
	if (validate(refGeneFileName).empty()){
		cout<<"Error Code(5)"<<endl; //hg19 refseq file not in refGene
		cout<<"Exiting"<<endl;
	}


	string linkedFileName		= "linked_" + split(BedGraphFile, "/")[split(BedGraphFile, "/").size()-1] + ".dat";
	linkedFileName				= out_dir + "linked_" + split(BedGraphFile, "/")[split(BedGraphFile, "/").size()-1] + ".dat";
	linkFile(BedGraphFile,strand, refGeneFileName,linkedFileName, 0, 0);
	if (validate(linkedFileName).empty()){
		cout<<"Error Code(6)"<<endl; //Linked binary was not created
		cout<<"Exiting..."<<endl;
		return 0;
	}
	vector<MEMM> struct_array	= readInBinFile(linkedFileName);
	cout<<"Reading Binary...check"<<endl;
	vector<int> codes;
	codes.push_back(1);
	codes.push_back(2);


	if (struct_array.empty()){
		cout<<"Error Code(7)"<<endl; //Was unable to read In binary
		cout<<"Exiting..."<<endl;
	}

	training_data results 	= fetchTrainingData(Training_File, BedGraphFile, struct_array, strand, codes, out_dir, refGeneFileName, 0);
	if (results.X.empty()){
		cout<<"Error Code(8)"<<endl; //Was unable to retrieve training data
		cout<<"Exiting..."<<endl;
	}

	vector<double>	W	= learn(results.X, results.Y,0,1);
	if (W.empty()){
		cout<<"Error Code(9)"<<endl; //Newton-Raphons didn't learn
		cout<<"Exiting..."<<endl;
	}

	cout<<"Raphson-Newton...check"<<endl;


	vector<logistic_equation> emissions(2);

	emissions[0] 			= logistic_equation(1.0, W, 0);
	emissions[1] 			= logistic_equation(1.0, W, 1);

	vector<MEMM> training_struct_array = getBWStruct_array(Training_File,out_dir);



	vector<vector<double> > A_matrix = runBaumWelch(training_struct_array, emissions, codes, 0,5);

	A_matrix[0][0]	= 0.99;
	A_matrix[0][1]	= 0.01;
	A_matrix[1][0]	= 0.01;
	A_matrix[1][1]	= 0.99;

	if (A_matrix.empty()){
		cout<<"Error Code(10)"<<endl;
		cout<<"Exiting..."<<endl;  //Baum Welch didn't work
		return 0;
	}
	cout<<"Baum Welch...check"<<endl;

	struct_array=viterbi(struct_array, A_matrix,  emissions, codes,0);
	if (struct_array.empty()){
		cout<<"Error Code(11)"<<endl; //viterbi didn't run
		cout<<"Exiting..."<<endl;
		return 0;
	}



	string userFileName=out_dir+"TEST_IGV_Segments.bed";

	writeOutFile(struct_array, BedGraphFile, out_dir, userFileName, out_dir);
	if (validate(userFileName).empty()){
		cout<<"Error Code(12)"<<endl;
		cout<<"Exiting..."<<endl;
		return 0;
	}

	cout<<"Path to Segmentation Package: "<<pathToSrc<<endl;

	cout<<"********************** All tests passed **********************"<<endl;
	return 1;
}




#endif /* EXAMPLE_H_ */
