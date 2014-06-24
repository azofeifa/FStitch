/*
 * readInParameters.h
 *
 *  Created on: Jan 23, 2014
 *      Author: joeyazo
 */
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include "MEMM_class.h"
#include "readIn.h"
#include "viterbi.h"
#include "functions.h"
#include "split.h"
#include <numeric>
#include <map>
using namespace std;

#ifndef READINPARAMETERS_H_
#define READINPARAMETERS_H_
class param_struct{
public:
	string flag;
	vector<string> value;
	param_struct(string val){
		flag=val;
	}
};

vector<param_struct> readInParameters( char* argv[]){
	/*
	 * Parse Command Line Parameters
	 */


	map<string, string> PARAMS;

	PARAMS["-i"] 	= "-i";
	PARAMS["-st"] 	= "-st";
	PARAMS["-rf"] 	= "-rf";
	PARAMS["-t"] 	= "-t";
	PARAMS["-l"] 	= "-l";
	PARAMS["-w"] 	= "-w";
	PARAMS["-a"] 	= "-a";
	PARAMS["-al"] 	= "-al";
	PARAMS["-bw"] 	= "-bw";
	PARAMS["-r"] 	= "-r";
	PARAMS["-v"] 	= "-v";
	PARAMS["-o"] 	= "-o";
	PARAMS["-d"] 	= "-d";
	PARAMS["-uf"] 	= "-uf";
	PARAMS["-wl"]	= "-wl";
	PARAMS["-eRNA"]	= "-eRNA";
	PARAMS["-cd"]	= "-cd";
	PARAMS["-rl"]	= "-rl";


	PARAMS["--reset-learning"]				= "-rl";
	PARAMS["--feature-code"]				= "-cd";
	PARAMS["--predict-eRNA"]				= "-eRNA";
	PARAMS["--reset"] 						= "-r";
	PARAMS["--help"] 						= "-h";
	PARAMS["--run-unit-test"] 				= "-t";
	PARAMS["--refseq-file"]	 				= "-rf";
	PARAMS["--alpha-learning-rate"] 		= "-al";
	PARAMS["--baum-welch-iterations"] 		= "-bw";
	PARAMS["--verbose"] 					= "-v";
	PARAMS["--output-directory"] 			= "-o";
	PARAMS["--input-bedgraph-file"] 		= "-i";
	PARAMS["--logistic-regression-weights"] = "-w";
	PARAMS["--feature-dimension"]			= "-d";
	PARAMS["--name-of-IGV-file"]			= "-uf";
	PARAMS["--learn"]						= "-l";
	PARAMS["--self-transition-HMM-weight"]	= "-a";
	PARAMS["--write-linked"]				= "-wl";


	string parameter;
	bool parameter_bool;

	vector<param_struct> param_array;
	bool FOUND = 0;
	int pos;
	while (*argv){
		parameter = *argv;
		parameter_bool= PARAMS.find(parameter)!=PARAMS.end();
		if(FOUND and (not parameter_bool or atoi(parameter.c_str()))){
			param_array[param_array.size()-1].value.push_back(parameter);
		}
		else if (parameter_bool && pos < 2 ) {
			if (PARAMS.find(parameter) != PARAMS.end()){
				param_array.push_back(param_struct(PARAMS[parameter]));
				FOUND  = 1;
			}else{
				param_array.push_back(parameter);
			}


		}
		else{
			FOUND = 0;
		}
		argv = ++argv;
	}


	return param_array;

}



#endif /* READINPARAMETERS_H_ */
