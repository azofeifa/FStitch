/*
 * writeOutFIle.h
 *
 *  Created on: Jan 23, 2014
 *      Author: joeyazo
 */

#ifndef WRITEOUTFILE_H_
#define WRITEOUTFILE_H_
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include "MEMM_class.h"
#include "readIn.h"
#include "viterbi.h"
#include "functions.h"
#include "extract_file_from_path.h"
#include <numeric>
using namespace std;

void writeOutFile(vector<MEMM> struct_array, string fileName, string path, string userName, string out_dir){
	string bed_delimiter = "Bed";
	string line;
	string BEGIN_POS;
	string BEGIN_CHROM;
	int pos;
	int prev_call = -1;
	string annotation;
	string RGB;
	bool FOUND;



	FOUND 	= search(fileName, bed_delimiter,pos);
	if (FOUND)
	{
		string strand = struct_array[0].strand;
		if (strand.empty()){
			strand = "+";
		}
		if (userName.empty() ){
			fileName = replace(fileName, ".BedGraph", "") + "_segs_IGV.bed";
			fileName = path+extract(fileName).file;
		}else{
			fileName	= userName;
		}
		string header = ("track name=Strand_" +strand+
				"description=\"LogitGRO\" visibility=2 useScore=2 cgGrades=50 cgColour1=white cgColour2=yellow cgColour3=red height=30\n");
		ofstream writeFileHandle;
		writeFileHandle.open(fileName);
		writeFileHandle<<header;
		string INFO;
		double score=0.;
		double N_score=0.;
		double cov_total	= 0.;
		double density	= 0.;
		for (int i = 0; i < struct_array.size() ; i++){
			strand = struct_array[i].strand;
			if (strand.empty()){
				strand = "+";
			}

			if (prev_call == -1){
				BEGIN_CHROM = struct_array[i].chrom;
				BEGIN_POS = to_string(static_cast <long long> (struct_array[i].begin));
			}
			else if(prev_call != struct_array[i].SegCall or BEGIN_CHROM != struct_array[i].chrom ){
				if (prev_call){
					annotation = "ON="+ to_string(score / N_score);
					if (strand =="+"){
						RGB = "0,0,255";
					}else{
						RGB = "255,0,0";
					}
					INFO = "100";
				}
				else{
					annotation  = "OFF=" + to_string(score / N_score);
					RGB = "0,255,0";
					INFO = "500";
				}
				density 	= cov_total / (struct_array[i-1].begin - atof(BEGIN_POS.c_str()));
			       
				if (not  (BEGIN_CHROM != struct_array[i].chrom )){
					line = (BEGIN_CHROM +"\t" + BEGIN_POS + "\t"+ to_string(static_cast <long long> (struct_array[i-1].begin))
							+ "\t" + annotation + "\t" + INFO + "\t" +strand + "\t" + BEGIN_POS + "\t" + to_string(static_cast <long long> (struct_array[i-1].begin))
							+ "\t" + RGB +"\t" + struct_array[i].annotation+  "..." + to_string(density) +"\n");
				}else{
					line = (BEGIN_CHROM +"\t" + BEGIN_POS + "\t"+ to_string(static_cast <long long> (struct_array[i-1].begin))
												+ "\t" + annotation + "\t" + INFO + "\t" +strand + "\t" + BEGIN_POS + "\t" +  to_string(static_cast <long long> (struct_array[i-1].begin))
												+ "\t" + RGB +"\t" + struct_array[i].annotation+ "..." + to_string(density) + "\n");
				}
				writeFileHandle<<line;
				if ((BEGIN_CHROM != struct_array[i].chrom )){
					BEGIN_POS="0";
				}else{
					BEGIN_POS = to_string(static_cast <long long> (struct_array[i-1].begin));
				}
				BEGIN_CHROM = static_cast <string> (struct_array[i].chrom);

				score =0.;
				N_score=0.;
				cov_total=0.;
			}
			score+=struct_array[i].score;
			N_score+=1;
			cov_total+=struct_array[i].cov_number;
			prev_call =  struct_array[i].SegCall;
		}
		if (prev_call){
			annotation = "ON";
			if (strand =="+"){
				RGB = "0,0,255";
			}else{
				RGB = "255,0,0";
			}
			INFO = "100";
		}
		else{
			annotation  = "OFF";
			RGB = "0,255,0";
			INFO = "500";
		}
		line = (BEGIN_CHROM +"\t" + BEGIN_POS + "\t"+ to_string(static_cast <long long> (struct_array[struct_array.size()-1].end))
									+ "\t" + annotation + "\t" + INFO + "\t" +struct_array[struct_array.size()-1].strand + "\t" + BEGIN_POS + "\t" +  to_string(static_cast <long long> (struct_array[struct_array.size()-1].end))
									+ "\t" + RGB +"\t" + struct_array[struct_array.size()-1].annotation+ "\n");

		writeFileHandle<<line;



		writeFileHandle.close();

	}


}



#endif /* WRITEOUTFILE_H_ */
