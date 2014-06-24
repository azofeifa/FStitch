/*
 * SeperateByStrand.h
 *
 *  Created on: Feb 28, 2014
 *      Author: joeyazo
 */

#ifndef SEPERATEBYSTRAND_H_
#define SEPERATEBYSTRAND_H_
using namespace std;
vector<string> SeperateByStrand(string BedGraphFile,string out_dir){
	ifstream myfile(BedGraphFile);

	string filePos	= replace(BedGraphFile, ".BedGraph", ".pos.BedGraph");
	string fileNeg	= replace(BedGraphFile, ".BedGraph", ".neg.BedGraph");
	filePos			= out_dir+split(filePos, "/")[split(filePos, "/").size()-1];
	fileNeg			= out_dir+split(fileNeg, "/")[split(fileNeg, "/").size()-1];


	ofstream posHandle;
	ofstream negHandle;
	posHandle.open(filePos);
	negHandle.open(fileNeg);

	cout<<"...Seperating by Strand..."<<endl;
	vector<string> files(2);
	files[0] 	= filePos;
	files[1] 	= fileNeg;
	string line;
	vector<string> lineInfo;
	int cov;
	string chrom, start, stop;
	if (myfile.is_open()){
		while ( getline (myfile,line) ){
			lineInfo	= split(line, "\t");

			cov 	= atoi(lineInfo[3].c_str());
			chrom 	= lineInfo[0];
			start 	= lineInfo[1];
			stop 	= lineInfo[2];

			if (cov < 0){
				negHandle<<(chrom+"\t" + start+ "\t"+ stop+ "\t" + to_string(abs(cov))+"\n");
			}
			else{
				posHandle<<(chrom+"\t" + start+ "\t"+ stop+ "\t" + to_string(abs(cov))+"\n");
			}
		}
	}
	return files;
}





#endif /* SEPERATEBYSTRAND_H_ */
