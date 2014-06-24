/*
 * validate_file.h
 *
 *  Created on: Jan 28, 2014
 *      Author: joeyazo
 */

#ifndef VALIDATE_FILE_H_
#define VALIDATE_FILE_H_

#include <iostream>
#include <fstream>
#include <string>
using namespace std;
string validate(string fileName){

	ifstream ifile(fileName);
	if (not ifile){
		return "";
	}
	else{
		return fileName;
	}
}







#endif /* VALIDATE_FILE_H_ */
