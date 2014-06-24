//============================================================================
// Name        : grab_segments.cpp
// Author      : Joey Azofeifa
// Version     : 1.0
// Copyright   : Your copyright notice
// Description : Main file for running Segmentation Package
//============================================================================
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <random>
#include "split.h"
#include <map>
using namespace std;

#ifndef GRAB_SEGMENTS_H_
#define GRAB_SEGMENTS_H_
#include "split.h"
map<string, int> chromProgress(){
	string pathToSrc = "/"+join(split(__FILE__, "/"), split(__FILE__, "/").size()-2,"/");

	string refFileName = pathToSrc+"refSeqAnnotations/chromosomeLengths.txt";
	ifstream segmentFH(refFileName);
	string line;
	vector<string> line_array;
	map <string, int> chromosomes;
	string chrom;
	int pos;
	if (segmentFH.is_open()){
		while ( getline(segmentFH,line) ){
			line_array = split(line, "\t");
			chrom 	= "chr" + line_array[0];
			pos		= atoi(line_array[1].c_str());
			if (chromosomes.find(chrom) == chromosomes.end()){
				if (chromosomes[chrom] < pos){
					chromosomes[chrom]	= pos;
				}
			}
			else{
				chromosomes[chrom] 	 	= pos;
			}
		}
	}
	else{
		cout<<"Couldnt Open File; Progress Bar: Disabled"<<endl;
		return chromosomes;
	}
	return chromosomes;
}

class segment{
public:
	vector<double> sense_coordinates;
	vector<double> sense_coverages;
	vector<double> anti_coverages;
	vector<double> anti_coordinates;

	string chrom;
	string strand;
	string GENE_ID;
	void add_sense(double coord, double cov){
		sense_coordinates.push_back(coord);
		sense_coverages.push_back(cov);
	}
	void add_anti(double coord, double cov){
		anti_coordinates.push_back(coord);
		anti_coverages.push_back(cov);
	}


	void display_info(){
		if (sense_coordinates.size()){
			cout<<"Segment:\t"<<sense_coordinates[0]<<"-"<<sense_coordinates[sense_coordinates.size()-1]<<"\tSize:\t"<<sense_coordinates.size()<<endl;
		}
	}
	segment(string s, string c,string g){
		strand	= s;
		chrom	= c;
		GENE_ID	= g;
	};
};


vector<segment> grab(string segmentFile, string sortedBedSumFile, bool test, bool SHOW){

	map<string, int> chromosomes = chromProgress();
	string line;
	ifstream segmentFH(segmentFile);
	ifstream BedSumFH(sortedBedSumFile);
	vector<string> lineInfo;
	vector<segment> segment_array;
	int start, stop;
	string chrom, strand;
	int start_BD, stop_BD,coverage_BD;
	string chrom_BD, strand_BD;
	string call;
	string GENE_ID;
	bool header = 1;
	int i =0;
	int fileStart, fileStop;
	fileStart = BedSumFH.tellg();
	int ct=0;
	string prev_chrom = "";
	int progress_counter = 0;
	double rand;
	double prev_percent=0.0;
	uniform_real_distribution<double> distribution(0.0,1.0);
	default_random_engine generator;

	if (segmentFH.is_open()){
		while ( getline(segmentFH,line) ){

			if (test){
				if (ct > 1){
					break;
				}
			}

			if (!header){
				lineInfo= split(line, "\t");
				chrom	= lineInfo[0];
				strand	= lineInfo[5];
				start	= atoi(lineInfo[1].c_str());
				stop	= atoi(lineInfo[2].c_str());

				if (chrom!=prev_chrom){
					if (prev_percent and SHOW){
						cout<<" * "<<flush;
					}
					if (SHOW){
						cout<<"\n...Collecting " <<chrom<<" Segments..."<<endl;
					}
					prev_chrom = chrom;
					prev_percent = 0.0;
				}

				GENE_ID	= lineInfo[9];

				if (lineInfo[3] == "ON"){
					lineInfo.clear();

					segment 	current(strand, chrom, GENE_ID);
					getline(BedSumFH, line);
					lineInfo	= split(line, "\t");

					start_BD	= atoi(lineInfo[1].c_str());
					chrom_BD	= lineInfo[0];
					strand_BD	= lineInfo[5];
					progress_counter+=1;

					if ((( (double)start_BD / chromosomes[chrom_BD]) - 0.1) > prev_percent and SHOW){
						prev_percent = ((double) start_BD / chromosomes[chrom_BD]);
						cout<<" * "<<flush;
					}
					while (chrom_BD != chrom){
						getline(BedSumFH, line);
						lineInfo	= split(line, "\t");
						chrom_BD	= lineInfo[0];
						lineInfo.clear();
					}

					while (start_BD < start){
						getline(BedSumFH, line);
						lineInfo	= split(line, "\t");
						start_BD	= atoi(lineInfo[1].c_str());
						chrom_BD	= lineInfo[0];
						lineInfo.clear();
					}
					while (1){
						getline(BedSumFH, line);
						lineInfo 	= split(line, "\t");
						start_BD	= atoi(lineInfo[1].c_str());
						chrom_BD	= lineInfo[0];
						coverage_BD	= atoi(lineInfo[4].c_str());

						if (lineInfo[5] == strand_BD ){
							current.add_sense(start_BD, coverage_BD);
						}
						else if (lineInfo[5] != strand_BD){
							current.add_anti(start_BD, coverage_BD);
						}

						if (start_BD > stop or chrom_BD != chrom){
							lineInfo.clear();
							break;
						}
						lineInfo.clear();
					}
					if (distribution(generator) < 0.5){
						ct+=1;
					}

					segment_array.push_back(current);
				}
				lineInfo.clear();
			}
			else{
				header = 0;
			}
		}
	}
	else{
		cout<<"Could not open: "<<segmentFile<<endl;
	}
	return segment_array;
}


#endif /* GRAB_SEGMENTS_H_ */
