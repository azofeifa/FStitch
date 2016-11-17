#include "validate.h"
#include <string>
#include <locale>
#include <fstream>
#include <iostream>

using namespace std;

bool isNum(string num){
	for (int i=0; i < num.size(); i++){
		if ((not isdigit(num[i])) and num.substr(i,1) != "." ){
			return false;
		}
	}
	return true;
}
bool isFile(string FILE){
	ifstream f(FILE);
	if (f.good()) {
		f.close();
		return true;
	}else{
		f.close();
		return false;
	}   
}
