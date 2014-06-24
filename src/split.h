/*
 * split.h
 *
 *  Created on: Jan 25, 2014
 *      Author: joeyazo
 */

#ifndef SPLIT_H_
#define SPLIT_H_
#include <iostream>
#include <vector>
#include <string>
using namespace std;


int search(string str, string delimiter, int &pos_ptr){
	const char * string_array = str.c_str();
	const char * del_array = delimiter.c_str();
	int N			= str.size();
	int pos 		= 0;
	int N_DEL		= delimiter.size();

	for (int i =0; i < N; i++){
		if ((string_array[i]==del_array[pos]) and pos < N_DEL){
			pos++;
		}
		else if (pos == N_DEL){
			pos_ptr = i;
			return 1;
		}else if(pos >0 && string_array[i]!=del_array[pos]){
			pos=0;
		}
	}
	if (pos == N_DEL){
		pos_ptr = str.size();
		return 1;
	}
	return 0;

}
bool SEARCH(string str, string delimiter){
	const char * string_array 	= str.c_str();
	const char * del_array 		= delimiter.c_str();
	int N			= str.size();
	int pos 		= 0;
	int N_DEL		= delimiter.size();
	for (int i =0; i < N; i++){
		if ((string_array[i]==del_array[pos]) and pos < N_DEL){
			pos++;
		}
		else if (pos == N_DEL){
			return 1;
		}
		else if(pos >0 && string_array[i]!=del_array[pos]){
			pos=0;
		}
	}
	if (pos == N_DEL){
		return 1;
	}
	return 0;

}

string replace(string str, string delimiter, string replacement){
	const char * string_array = str.c_str();
	const char * del_array = delimiter.c_str();
	int N			= str.size();
	int pos 		= 0;
	int N_DEL		= delimiter.size();
	int beg 		= 0;
	int end			= 0;
	int i			= 0;
	for (i =0; i < N; i++){
		if ((string_array[i]==del_array[pos]) and pos < N_DEL){
			if (pos == 0){
				beg=i;
			}
			pos++;
		}
		else if (pos == N_DEL){
			end = i;
			break;
		}
		else{
			pos = 0;
		}
	}
	if (pos==N_DEL){
		end=i;
	}
	return str.substr(0,beg) + replacement + str.substr(end,N);

}


vector<string> split(string str, string delimiter){
	int N			= str.size();
	string curr 	= "";
	vector<string> info;
	const char * string_array = str.c_str();
	const char * del_array = delimiter.c_str();
	for (int i =0; i < N; i++){
		if ((string_array[i]!=*del_array)){
			curr	= curr+string_array[i];
		}
		else if (!curr.empty()){
			info.push_back(curr);
			curr="";
		}
	}
	if (!curr.empty()){
		info.push_back(curr);
	}
	return info;
}

string join(vector<string> string_vector, int threshold = 0, string delimiter = ""){
	if (not threshold){
		threshold =string_vector.size();
	}
	string ss = "";
	for (int i=0; i < threshold; i++){
		ss = ss  + string_vector[i]+delimiter;
	}
	return ss;
}







#endif /* SPLIT_H_ */
