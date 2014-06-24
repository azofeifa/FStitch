//============================================================================
// Name        : main.cpp
// Author      : Joey Azofeifa
// Version     : 1.0
// Copyright   : Your copyright notice
// Description : Main file for running Segmentation Package
//============================================================================
#include <iostream>
#include <fstream>

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
#include "help.h"
#include "MAP_Decoding.h"
#include "baum_welch.h"
#include "eRNA.h"
#include "modelParameters.h"
using namespace std;


int main(int argc, char* argv[]) {

	vector<param_struct> param_array = readInParameters(argv);
	string pathToSrc = "/"+join(split(__FILE__, "/"), split(__FILE__, "/").size()-2,"/");
	//===================================================================
	//	Some Parameters
	bool link 			= 1;
	bool reset 			= 0;
	bool test 			= 0;
	bool write_linked 	= 0;
	bool SHOW			= 0;
	bool TRAIN 			= 0;
	//===================================================================

	//====================================================================
	//	User Parameters
	string fileName 		= "";
	string refGeneFileName 	= "";
	string Training_File	= "";
	double slope_density 	= 0.01; //Logisitic Regression Parameter Density
	double slope_coverage 	= 2; //Logisitic Regression Parameter Coverage
	double weight 			= 1.0; //Weight Coverage Logistic  vs Density Logistic
	double A 				= 0.99; //	Self Transition Loop Parameter
	bool example 			= 0;
	bool strand_bool 		= 0;
	string strand 			= "";
	int dim 				= 3;
	bool exist 				= 0;
	bool viterbi_BOOL 		= 1;
	bool RESET_LEARN		= 0;



	//Model Parameter data structures
	vector<double> W;
	vector<vector<double>> A_matrix;

	string userFileName		= "";
	int BM					= 10;
	double alpha 			= 1.;
	string eRNA_file1="";
	string eRNA_file2="";
	double eRNA_Thresh		= 0.3;
	bool runERNA			= 0;
	vector<int> feature_codes;


	string out_dir = pathToSrc + "segmentation_outputs/";
	vector<string> unknowns;
	//====================================================================


	////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (not param_array.size()){
		runHelp();
		return 1;
	}
	for (int i = 0; i< param_array.size();i++){

		if (param_array[i].flag == "-h"){
			runHelp();
			return 1;
		}

		else if(param_array[i].flag == "-l"){
			TRAIN			= 1;
			Training_File 	= param_array[i].value[0];
			string fileNameBool = validate(Training_File);
			if (fileNameBool.empty()){
				cout<<"user specified training file: "<< Training_File<< " does not exist..."<<endl;
				cout<<"please provide correct path or place in Segmentation/src/ exectuable directory..." << endl;
				cout<<"Exiting..."<<endl;
				return 0;
			}
		}
		else if(param_array[i].flag=="-cd"){
			for (int j =0;j < param_array[i].value.size();j++ ){
				feature_codes.push_back(atoi(param_array[i].value[j].c_str()));
			}
		}


		else if(param_array[i].flag == "-eRNA"){
			if (param_array[i].value.size() < 2){
				cout<<("Please specify two IGV classification files positive and negative strand\n"
						"These files are the output of running ./segment -i file1.pos.BedGrpah and ./segment -i file2.neg.BedGraph\n"
						"and sould be in your specified output directory (default being segmentation outputs)"
						)<<endl;
				cout<<"Exiting..."<<endl;
				return 0;
			}else{
				eRNA_file1=pathToSrc+param_array[i].value[0];
				eRNA_file2=pathToSrc+param_array[i].value[1];

				if (param_array[i].value.size() > 2){
					eRNA_Thresh=atof(param_array[i].value[2].c_str());
				}

				runERNA = 1;


			}

		}
		else if(param_array[i].flag=="-rl"){
			TRAIN			= 1;
			Training_File 	= param_array[i].value[0];
			string fileNameBool = validate(Training_File);
			if (fileNameBool.empty()){
				cout<<"user specified training file: "<< Training_File<< " does not exist..."<<endl;
				cout<<"please provide correct path or place in Segmentation/src/ exectuable directory..." << endl;
				cout<<"Exiting..."<<endl;
				return 0;
			}
			RESET_LEARN=1;

		}

		else if(param_array[i].flag == "-uf"){
			userFileName	= param_array[i].value[0];
		}
		else if(param_array[i].flag =="-t"){
			example = 1;
			refGeneFileName = pathToSrc +"refSeqAnnotations/refSeqGene_hg19.txt";
			out_dir = pathToSrc + "examples/";
			fileName = pathToSrc + "examples/" +"TEST_DMSO2_3.pos.BedGraph";
			cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
			cout<<"Welcome new User! We are running some unit tests to make sure"<<endl;
			cout<<"everything compiled okay. If any error codes arise please look them up"<<endl;
			cout<<"in the ErrorCodeBook.txt"<<endl;
			cout<<"Please direct questions on use to Joey at jgazofeifa [at] gmail [dot] com"<<endl;
			cout<<"Or look up in README"<<endl;
			cout<<"Happy Annotating!"<<endl;
			cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

			unitTest(pathToSrc);
			return 1;

		}
		else if (param_array[i].flag == "-bw"){
			BM	= atoi(param_array[i].value[0].c_str());
		}
		else if (param_array[i].flag == "-i"){
			fileName = param_array[i].value[0];
			string fileNameBool = validate(fileName);
			if (fileNameBool.empty()){
				cout<<"user specified .BedGraph file: "<< fileName<< " does not exist..."<<endl;
				cout<<"please provide correct path or place in Segmentation/src/ exectuable directory..." << endl;
				cout<<"Exiting..."<<endl;
				return 0;
			}
		}


		else if (param_array[i].flag == "-d"){
			dim = atoi(param_array[i].value[0].c_str());
		}

		else if (param_array[i].flag == "-w"){
			for (int j =0; j < param_array[i].value.size(); j++){
				W.push_back(atof(param_array[i].value[j].c_str()));
			}

		}
		else if (param_array[i].flag == "-a"){
			if (atof(param_array[i].value[0].c_str())){
				A = atof(param_array[i].value[0].c_str());
			}
		}
		else if (param_array[i].flag == "-rf"){
			refGeneFileName = param_array[i].value[0];
			refGeneFileName = validate(refGeneFileName);
			if (refGeneFileName.empty()){
				cout<<"Exiting..."<<endl;
				return 0;
			}
		}
		else if (param_array[i].flag == "-o"){
			out_dir = param_array[i].value[0];
			string out_dir_bool = validate(out_dir);
			if (out_dir_bool.empty()){
				cout<<"User specified output directory: "<<out_dir<<", please provide correct path"<<endl;
				cout<<"Exiting..."<<endl;
				return 0;
			}
		}
		else if (param_array[i].flag == "-r"){
			reset=1;
		}
		else if (param_array[i].flag == "-wl"){
			write_linked=1;
		}
		else if (param_array[i].flag == "-v"){
			SHOW=1;
		}else if(param_array[i].flag=="-al"){
			alpha=atof(param_array[i].value[0].c_str());
		}
		else{
			unknowns.push_back(param_array[i].flag);
		}

	}
	if (runERNA){
		cout<<"###############################################"<<endl;
		cout<<"User Parameter Values for LogitGRO"<<endl<<endl;
		cout<<"Running eRNA prediction"<<endl;
		cout<<"File Sense Strand: "<<eRNA_file1<<endl;
		cout<<"File Anti-Sense Strand: "<<eRNA_file2<<endl;
		cout<<"eRNA confidence threshold: "<<eRNA_Thresh<<endl;
		cout<<"Output Directory: "<<out_dir<<endl;
		cout<<"###############################################"<<endl;
		eRNApredictions(eRNA_file1,eRNA_file2,eRNA_Thresh, out_dir);
		return 1;
	}

	if (feature_codes.empty()){
		feature_codes.push_back(1);
		feature_codes.push_back(2);
	}




	if (validate(out_dir).empty()){
		cout<<"Specified output directory does not exist..."<<endl;
		cout<<"Either specify the correct path a directory using -o or run 'segment' within the src directory..."<<endl;
		cout<<"Exiting..."<<endl;
		return 0;
	}
	string BedGraphFile 		= fileName;

	string linkedFileName		= "linked_" + split(BedGraphFile, "/")[split(BedGraphFile, "/").size()-1] + ".dat";
	string fileNamePath 		= join(split(BedGraphFile, "/"),split(BedGraphFile, "/").size()-1, "/");
	string parameterFile= "";


	vector<string> paths ;
	paths.push_back(fileNamePath);
	paths.push_back("./");
	paths.push_back(out_dir);
	paths.push_back("");


	paths.push_back(pathToSrc+"segmentation_outputs/");
	if (not Training_File.empty()){

		parameterFile = split(Training_File, "/")[split(Training_File, "/").size()-1];
		parameterFile= to_string(alpha) + "_" + to_string(BM)+"_" + parameterFile;
		exist	= doesExist(parameterFile, paths, SHOW);
		parameterFile	= search(paths,parameterFile);
		if (RESET_LEARN){
			if( remove( search(paths, parameterFile).c_str() ) != 0 ){
				if (SHOW){
					perror( "Error deleting file" );
				}
			}
			else{
				if(SHOW){
					puts( "Parameter File successfully deleted" );
				}
			}
		}
	}




	if (fileName.empty()){
		cout<<"Please Specify Input BedGraph File..."<<endl;
		cout<<"Exiting..."<<endl;
		return 0;
	}

	if (TRAIN and (Training_File.empty())){
		cout<<"Using Specified learning but did not provide paths using -ON and -OFF to the training dataset..."<<endl;
		cout<<"Warning: Assuming Default Parameters..."<<endl;
		TRAIN 	= 0;
		return 0;
	}

	if (not TRAIN and W.empty()){
		W.push_back(0.24);
		W.push_back(0.005);
	}


	string comp = string(&out_dir[out_dir.size()-1]);
	int pos;
	if (not search(comp, "/", pos)){
		out_dir = out_dir + "/";
	}

	if (SEARCH(fileName, "+") || SEARCH(fileName, "pos")){
		strand = "+";
	}
	else if (SEARCH(fileName, "-") || SEARCH(fileName, "neg")){
		strand = "-";

	}



	cout<<"###############################################"<<endl;
	cout<<"User Parameter Values for LogitGRO"<<endl<<endl;
	cout<<"Input Bed Graph FILE: "<<fileName<<endl;
	if (refGeneFileName.empty()){
		refGeneFileName = pathToSrc +"refSeqAnnotations/refSeqGene_hg19.txt";
		cout<<"RefGeneFile: "<<refGeneFileName<<" (default user did not specify...)"<<endl;
	}
	else{
		cout<<"RefGeneFile: "<<refGeneFileName<<endl;
	}
	if (out_dir.empty()){
		cout<<"Outputting Files to: Current Directory (default user did not specify...)"<<endl;
	}
	else{
		cout<<"Outputting Files to: "<<out_dir<<endl;
	}
	if (strand.empty()){
		cout<<"Strand:	+ (Warning: Couldn't Infer Strand ID from file name of given bed graph, assuming default of +)"<<endl;
	}
	else{
		cout<<"Strand: "<<strand<<endl;
	}

	cout<<"Learning rate:\t"<<alpha<<endl;
	cout<<"Iterations for Baum-Welch:\t"<<BM<<endl;

	map<int, string> feature_disp_map;
	feature_disp_map[1]= "length";
	feature_disp_map[2]= "mean";
	feature_disp_map[3]= "mode";
	feature_disp_map[4]= "max";
	feature_disp_map[5]= "min";
	feature_disp_map[6]= "cov_number";
	feature_disp_map[7]= "density";
	feature_disp_map[8]= "variance";
	feature_disp_map[9]= "GENE";

	cout<<"Features being considered: ";
	for (vector<int>::iterator fc=feature_codes.begin();fc!=feature_codes.end();++fc){
		cout<<feature_disp_map[*fc]<<", ";
	}
	cout<<endl;



	if (not TRAIN){
		cout<<"Self Transition HMM Loop: "<<A<<endl<<endl;
		cout<<"Logistic Regression Parameters: ";
		for (int i = 0 ; i <W.size(); i++){
			cout<<"\t" + to_string(W[i])<<flush;
		}
		cout<<endl;


	}
	cout<<"###############################################"<<endl;
	if (unknowns.size()){
		cout<<"Ignoring User Parameters: ";
		for (int i = 0; i < unknowns.size();i++){
			cout<<unknowns[i]<<" "<<flush;
		}
		cout<<endl;
	}





	//==============================================================
	//	Look for linked_ file; first in the out_dir and then the root
	//	of the input bed graph file
	//==============================================================





	exist = 0;

	if (TRAIN){

		//Check to see if the training file passed by the user is the learning file
		exist	= isTrainingFile(Training_File);
		if (not exist){

			parameterFile = split(Training_File, "/")[split(Training_File, "/").size()-1];
			parameterFile= to_string(alpha) + "_" + to_string(BM)+"_" + parameterFile;
			exist	= doesExist(parameterFile, paths, SHOW);
			parameterFile	= search(paths,parameterFile);
		}else{
			parameterFile=Training_File;
		}



		if (exist){
			TRAIN=0;
			parameters P = readInModelParameters(parameterFile);
			W			= P.W;
			A_matrix 	= P.A_matrix;
			cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
			cout<<"...model parameters learned from before..."<<endl;
			cout<<"Logitistic Regresion Parameters:   ";
			for (vector<double>::iterator w = W.begin(); w!=W.end(); ++w){
				cout<<*w<<",";
			}
			cout<<endl;
			cout<<"Markov Model Transition Matrix"<<endl;
			for (vector<vector<double>>::iterator A_i = A_matrix.begin(); A_i != A_matrix.end(); ++A_i){
				for (vector<double>::iterator A_j = A_i->begin(); A_j!= A_i->end(); ++A_j ){
					cout<<*A_j<<",";
				}
				cout<<endl;
			}
			cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		}
	}



	vector<MEMM> struct_array;

	linkedFileName="linked_" + split(BedGraphFile, "/")[split(BedGraphFile, "/").size()-1] + ".dat";
	if (reset){
		if( remove( search(paths, linkedFileName).c_str() ) != 0 ){
			if (SHOW){
				perror( "Error deleting file" );
			}
		}
		else{
			if(SHOW){
				puts( "Linked File successfully deleted" );
			}
		}
	}
	if (search(paths, linkedFileName).empty()){
		linkFile(BedGraphFile,strand, refGeneFileName,out_dir+linkedFileName, write_linked, SHOW);
	}
	else if(SHOW){
		cout<<"...File "<<(search(paths, linkedFileName))<<" already exists..."<<endl;
	}

	struct_array	= readInBinFile(search(paths, linkedFileName));



	if (not W.empty()){
		dim			= W.size();
	}




	if (SHOW){
		cout<<"...Setting Up Parameters to Run Viterbi..."<<endl;
	}

	string a	= to_string(A); //	Self Transition Loop Parameter
	string d 	= to_string(dim);

	a.erase ( a.find_last_not_of('0') + 1, string::npos );



	training_data results;

	if (TRAIN){
		if (SHOW){
			cout<<"...Grabbing Training Data to learn logistic regression parameters..."<<endl;
		}



		if (struct_array.empty()){
			cout<<"Error (Critical): Linked Class Array was not populated..."<<endl;
			cout<<"Exiting..."<<endl;
			return 0;
		}


		results 	= fetchTrainingData(Training_File, BedGraphFile, struct_array, strand, feature_codes, out_dir,refGeneFileName, SHOW);

		if (results.Y.empty()){
			cout<<"Error (Critical): training data was not gathered..."<<endl;
			cout<<"Exiting..."<<endl;
			return 0;
		}

		W						= learn(results.X, results.Y, SHOW,alpha);
	}

	vector<logistic_equation> emissions(2);

	emissions[0] 			= logistic_equation(weight, W, 0);
	emissions[1] 			= logistic_equation(weight, W, 1);

	int num_of_col = 2;
	int num_of_row = 2;

	if (not exist){

		A_matrix.resize( num_of_col , vector<double>( num_of_row ) );

		A_matrix[0][0] = A;
		A_matrix[0][1] = 1.0-A;
		A_matrix[1][0] = 1.0-A;
		A_matrix[1][1] = A;
	}
	if (TRAIN){
		parameterFile = split(Training_File, "/")[split(Training_File, "/").size()-1];
		parameterFile= to_string(alpha) + "_" + to_string(BM)+"_" + parameterFile;
		vector<MEMM> training_struct_array = getBWStruct_array(Training_File,out_dir);


		A_matrix = runBaumWelch( training_struct_array, emissions, feature_codes, SHOW, BM);
		if (SHOW){
			cout<<"writing learned parameters to: "<<out_dir+parameterFile<<endl;
		}
		writeOutModelParameters(out_dir+parameterFile, A_matrix, W);
	}

	if (struct_array.empty()){
		cout<<"MEMM Struct Array was not populated..."<<endl;
		cout<<"Exiting..."<<endl;
		return 0;
	}

	if (viterbi_BOOL){
		//===================================================================
		// 	Viterbi
		//===================================================================
		if (SHOW){
			cout<<"...about to run viterbi..."<<endl;
		}
		struct_array=viterbi(struct_array, A_matrix,  emissions, feature_codes,SHOW);
	}
	else{
		cout<<"...decoding method was not specified..."<<endl;
		cout<<"Exiting..."<<endl;
		return 0;
	}

	writeOutFile(struct_array, BedGraphFile, out_dir , userFileName, out_dir);












	return 0;
}
