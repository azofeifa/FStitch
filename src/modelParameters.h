/*
 * modelParameters.h
 *
 *  Created on: May 10, 2014
 *      Author: joeyazo
 */

#ifndef MODELPARAMETERS_H_
#define MODELPARAMETERS_H_

bool doesExist(string parameterFile,vector<string> paths, bool SHOW){

	bool exists 	= search(paths,parameterFile).empty();
	if (SHOW){
		if (not exists){
			cout<<parameterFile<<" exists about to readin"<<endl;
		}else{
			cout<<parameterFile<<" doesn't exist...learning paramteres"<<endl;
		}
	}
	return not exists;
}

void writeOutModelParameters(string parameterFile, vector< vector<double> > A_matrix, vector<double> W){
	ofstream writeFileHandle;
	writeFileHandle.open(parameterFile);
	writeFileHandle<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	for (int i =0; i < A_matrix.size(); i++){
		for (int j = 0; j < A_matrix.size(); j++){
			writeFileHandle<<to_string(A_matrix[i][j]) + ",";
		}
		writeFileHandle<<endl;
	}
	writeFileHandle<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	for (int i = 0 ; i < W.size(); i++){
		writeFileHandle<<to_string(W[i]) + ",";
	}
	writeFileHandle<<endl;
}


class parameters{
public:
	vector< vector<double>> A_matrix;
	vector<double>	W;
	parameters(vector<vector<double>> A, vector<double> w){
		A_matrix =	 A;
		W		=	 w;
	}
	parameters(){};
};

parameters readInModelParameters(string parameterFile){
	ifstream fileHandle;
	fileHandle.open(parameterFile);
	string line;
	int i=0;
	int j = 0;
	int k =0;
	vector<vector<double>>A;
	vector<double>	W;
	vector<string> lineInfo;
	parameters P;
	if (fileHandle.is_open()){
		while ( getline (fileHandle,line) ){
			if (not SEARCH(line, "~") and i == 0){
				return P;
			}
			else if (SEARCH(line, "~")){
				i++;
			} else if (not SEARCH(line, "~") and i == 1){
				//transition matrix
				lineInfo	= split(line, ",");
				A.push_back(vector<double>(lineInfo.size()));
				k=0;
				for (vector<string>::iterator val = lineInfo.begin(); val!= lineInfo.end(); ++val){
					A[j][k]	= atof(val->c_str());
					k++;
				}
				j++;
			}else if (not SEARCH(line,"~") and i == 2){
				lineInfo	= split(line, ",");

				for (vector<string>::iterator val = lineInfo.begin(); val!= lineInfo.end(); ++val){
					W.push_back(atof(val->c_str()));
				}
			}
		}
	}else{
		cout<<"couldn't open: "<<parameterFile<<endl;
	}

	P.A_matrix	= A;
	P.W			= W;
	return P;
}

bool isTrainingFile(string file){
	ifstream fileHandle;
	fileHandle.open(file);
	string line;
	if (fileHandle.is_open()){
		while ( getline (fileHandle,line) ){
			if (SEARCH(line, "~")){
				return 1;
			}else{
				return 0;
			}
		}
	}
	return 0;
}






#endif /* MODELPARAMETERS_H_ */
