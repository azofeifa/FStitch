/*
 * link.h
 *
 *  Created on: Jan 25, 2014
 *      Author: joeyazo
 */
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "split.h"
#include "MEMM_class.h"
#include <algorithm>
#include <time.h>
#include "readInRefGene.h"
#include "window.h"
#include <map>
using namespace std;
#ifndef LINK_H_
#define LINK_H_

double maximum(vector<double> array){
	double * curr_max = NULL;
	for (int i =0; i<array.size(); i++){
		if (curr_max == NULL){
			curr_max = &array[i];
		}
		else if (*curr_max < array[i]){
			curr_max = &array[i];
		}
	}
	return *curr_max;
}

double minimum(vector<double>array){
	double curr_min	= array[0];
	for (vector<double>::iterator a = array.begin(); a!=array.end(); ++a){
		if (*a < curr_min){
			curr_min=*a;
		}
	}
	return curr_min;

}


double mean(vector<double> array){
	double N 	= array.size();
	double sum	= 0;
	for (int i =0;i< N; i++){
		sum+=array[i];
	}
	return sum / N;
}

double median(vector<double> array){
	double N 	= array.size();
	int j 		= (int) N/2;
	return array[j];
}
double mode(vector<double> array){
	map <double, int> values;
	int N 	= array.size();
	for (int i =0;i< N; i++){
		if (not (values.find(array[i]) == values.end())){
			values[array[i]] = 0;
		}
		values[array[i]]+=1;
	}
	int curr_max	 	= 0;
	double curr_num		= 0;
	for(map<double,int>::iterator it = values.begin(); it != values.end(); ++it) {
	  if (it->second > curr_max){
		  curr_max 	= it->second;
		  curr_num	= it->first;
	  }
	  else if (it-> second == curr_max and it->first > curr_num){
		  curr_num	= (curr_num+  it->first) /2.0;
	  }
	}
	return curr_num;


}



double sample_variance(vector<double> array, double mean){
	double sum = 0;
	for (int i =0; i < array.size(); i++){
		sum = sum + pow((array[i] - mean),2);
	}
	return (sum / (float)array.size());
}
double sum_array(vector<double> array){
	double SUM 	= 0;
	for (int i =0; i<array.size();i++){
		SUM+=array[i];
	}
	return SUM;
}
int max_array(vector<int> array){
	int MAX = 0;
	for (int i =0; i<array.size();i++){
		if (array[i] > MAX){
			MAX 	= array[i];
		}
	}
	return MAX;
}


int linkFile(string fileName, string input_strand, string refGeneFileName,string linkedFileName, bool write_linked, bool SHOW){
	clock_t time;
	int f;
	time = clock();

	ifstream myfile(fileName);
	tree * refGeneStructure = new tree();

	refGeneStructure  = readRefGene(refGeneFileName);

	if (SHOW){
		cout<<"...Linking File..."<<endl;
	}
	if (refGeneStructure->strand_ptr != NULL){
		if (myfile.is_open()){
			string line, chrom;
			string token;
			int start, stop, coverage;
			int prev = 0;
			vector<string> lineInfo;
			int pos;
			int query_pos;
			//=======================================
			//	Binary File Handle
			//=======================================
			string BinfileName = linkedFileName;
			string rplc = ".dat";
			ofstream writeFileHandleBinary;

			writeFileHandleBinary.open(BinfileName, ios::binary);
			string TextFileName = replace(BinfileName, rplc, ".bed");
			ofstream TextFileHandle;
			if (write_linked){
				//=======================================
				//	Text File Handle
				//=======================================
				string header = "chrom\tstart\tstop\tlength\tGENEID\tmean\tvariance\tmedian\tmode\tmax\tdensity\tmin\n";
				TextFileHandle.open(TextFileName);
				TextFileHandle<<header;
			}

			double sum;
			int cov;
			int N;
			vector<double> coverages;
			string GENE_ID;
			string length;
			string begin_chrom="";
			int begin = 0;
			double mean;
			string out_line;
			MEMM MEMM_obj;
			double v2,v3,v4,v5,v6,v7,v8;
			int SHOW_CT 	= 0;
			vector< vector<double> > 	window;
			int window_size	= 500;
			vector<MEMM> struct_array;
			int ERROR = 0;
			while ( getline(myfile,line) ){
				lineInfo	= split(line, "\t");

				//===========================================================//
				//	Index of lineInfo					  					 //
				//	0 = chromosome						  					 //
				// 	1 = start							  					 //
				//	2 = stop							  					 //
				//	3 = coverage across span			  					 //
				// 	Search = refGeneStructure, strand, chromosome, query_pos //
				//===========================================================//
				if (not (lineInfo.size() == 4) and ERROR < 1){ //Assume this is a header
					ERROR++;
				}
				else if (not (lineInfo.size() == 4)){
					cout<<"File: "<<fileName<<" is not formatted properly..."<<endl;
					cout<<"must be in tab seperated format..."<<endl;
					cout<<"Exiting..."<<endl;
					return 0;
				}
				else{
					start = atoi(lineInfo[1].c_str());
					stop = atoi(lineInfo[2].c_str());
					if (begin_chrom != lineInfo[0] and SHOW){
						cout<<"Linking..."<< lineInfo[0]<<"\t"<<flush;
						SHOW_CT+=1;
					}
					if (SHOW_CT > 3){
						cout<<endl;
						SHOW_CT=0;
					}


					if (prev){

						if (prev != start or begin_chrom != lineInfo[0]){
							query_pos = (prev + begin) / 2;
							GENE_ID = search(refGeneStructure, input_strand, lineInfo[0], query_pos);
							if (GENE_ID.empty()){
								GENE_ID = "N_A";
							}

							length = to_string(static_cast <long long> (prev - begin));
							mean = (float)sum / (float)N;
							v2	= median(coverages);
							v3	= mode(coverages);
							v4	= sample_variance(coverages, mean);
							v5	= maximum(coverages);
							v6	= sum_array(coverages);
							v7	= (sum / (prev - begin));
							v8	= minimum(coverages);
							if (write_linked){
								string v22 	= to_string(v2); //Logisitic Regression Parameter Density
								v22.erase ( v22.find_last_not_of('0') + 1, string::npos );
								string v32 	= to_string(v3); //Logisitic Regression Parameter Density
								v32.erase ( v32.find_last_not_of('0') + 1, string::npos );
								string v42 	= to_string(v4); //Logisitic Regression Parameter Density
								v42.erase ( v42.find_last_not_of('0') + 1, string::npos );
								string v52 	= to_string(v5); //Logisitic Regression Parameter Density
								v52.erase ( v52.find_last_not_of('0') + 1, string::npos );
								string v12 	= to_string(mean); //Logisitic Regression Parameter Density
								v12.erase ( v12.find_last_not_of('0') + 1, string::npos );
								string v72	= to_string(v7);
								v72.erase ( v72.find_last_not_of('0') + 1, string::npos );
								string v82	= to_string(v8);
								v82.erase ( v82.find_last_not_of('0') + 1, string::npos );




								//==================================================
								//	Writing to linked data to text file for
								//	other programs (python...) to interpret and parse
								//	This Linked data is the continuous segment with
								//	no interruptions in read connection
								//==================================================
								out_line = (begin_chrom + "\t" +  to_string(static_cast <long long> (begin)) + "\t" +
										to_string(static_cast <long long>  (prev)) + "\t" + length + "\t" + GENE_ID + "\t" + v12  + "\t" + v42 + "\t" + v22 + "\t" + v32 + "\t" +
										v52 + "\t" + v72 + "\t" + v82+ "\n");
								TextFileHandle<<out_line;
							}

							MEMM_obj = MEMM(input_strand , begin_chrom, begin, prev,  GENE_ID, (prev - begin), mean , maximum(coverages), sample_variance(coverages,mean));
							MEMM_obj.setStats(mean, v2, v3, v5, v4, v6, v7, v8);

							//==================================================
							//	Writing to Binary File for faster processing...
							writeFileHandleBinary.write((char*)&MEMM_obj, sizeof(MEMM_obj));
							sum	= 0;
							N 	= 0;


							length = to_string(static_cast <long long> ((start - prev) * -1));

							coverages.clear();
							if (write_linked){
								//==================================================
								//	Writing to linked data to text file for
								//	other programs (python...) to interpret and parse
								//	This is gap between when a continuous segment of
								// 	reads end and the next segment beings
								//	signified by the fact that its a negative length
								//==================================================
								out_line = lineInfo[0] + "\t" + to_string(static_cast <long long> (prev)) + "\t" + to_string(static_cast <long long>  (start)) + "\t" + length + "\t" + GENE_ID + "\t0\t0\t0\t0\t0\t0\n";
								TextFileHandle<<out_line;

							}

							MEMM_obj = MEMM(input_strand , begin_chrom, prev, start,  GENE_ID,((start - prev) * -1),0,0,0);
							//==================================================
							//	Writing to Binary File for faster processing...
							writeFileHandleBinary.write((char*)&MEMM_obj, sizeof(MEMM_obj));

							begin = start;
							begin_chrom = lineInfo[0];
						}
					}
					else{
						begin = start;
						begin_chrom = lineInfo[0];
					}
					prev = stop;

					N++;
					cov = atof(lineInfo[3].c_str());
					sum = sum+cov;
					coverages.push_back(cov);
					lineInfo.clear();
				}


			}



			refGeneStructure->clear(refGeneStructure);



			writeFileHandleBinary.close();
			TextFileHandle.close();


		}
		else
		{
			cout<<"Unable to open bed graph file: "<<fileName<<"...try specifying full path/spelling error correction"<<endl;
		}

	}
	if (SHOW){
		cout<<"\n...linked bed graph file:\t"<<((float)clock()-time)/CLOCKS_PER_SEC<<" seconds..."<<endl;
	}
	return 1;
}

#endif /* LINK_H_ */
