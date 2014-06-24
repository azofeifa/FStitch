/*
 * search_for_files.h
 *
 *  Created on: Jan 28, 2014
 *      Author: joeyazo
 */

#ifndef SEARCH_FOR_FILES_H_
#define SEARCH_FOR_FILES_H_
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

string search(vector<string> paths, string filename){
	for (int i = 0; i < paths.size(); i++){
		ifstream ifile(paths[i] + filename);
		if (ifile){
			return  paths[i] + filename;
		}
	}
	return "";


}




#endif /* SEARCH_FOR_FILES_H_ */
