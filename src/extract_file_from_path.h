/*
 * extract_file.h
 *
 *  Created on: Jan 28, 2014
 *      Author: joeyazo
 */

#ifndef EXTRACT_FILE_H_
#define EXTRACT_FILE_H_
#include <vector>
#include <string>
#include "split.h"
using namespace std;
class split_struct{
public:
	string file;
	string path;
	split_struct(string val, string val2){
		file = val;
		path = val2;

	}
};


split_struct extract(string path){
	vector<string> info = split(path, "/");
	string file = info[info.size()-1];
	string new_path ="";
	for (int i = 0 ; i < info.size()-1; i++){
		new_path = new_path + info[i] + "/";
	}
	return split_struct(file, new_path);
}








#endif /* EXTRACT_FILE_H_ */
