#ifndef read_in_parameters_H
#define read_in_parameters_H

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <map>
using namespace std;


class paramsTrain{
public:
	string ID="train";
	map<string, string> params;
	void display();
	paramsTrain();

};
class paramsSegment{
public:
	string ID= "segment";
	//acceptable user parameters
	map<string, string> params;
	paramsSegment();
	void display();
};
class paramsERNA{
public:
	string ID = "eRNA";
	map<string, string> params;
	string out;
	paramsERNA();
	void display();
};
class paramWrapper{
public:
	bool train 	=0;
	bool segment=0;
	bool eRNA 	=0;
	paramsTrain PT ; 
	paramsSegment PS; 
	paramsERNA PE ; 
	paramWrapper();
	void display();
};


void fillInOptions(char*,paramWrapper);

paramWrapper * readInParameters(char**);

#endif