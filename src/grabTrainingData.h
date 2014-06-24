/*
 * grabTrainingData.h
 *
 *  Created on: Feb 24, 2014
 *      Author: joeyazo
 */

#ifndef GRABTRAININGDATA_H_
#define GRABTRAININGDATA_H_
#include <map>
#include "split.h"
#include "readIn.h"
#include <random>
using namespace std;
class training_data{
public:
	vector< vector<double> >X;
	vector<int> 			Y;

	training_data(vector< vector<double> > V1, vector<int> V2){
		X	= V1;
		Y	= V2;
	}
	training_data(){};

};

//================================================
//	Possible Feature Set (Available in MEMM)
//	int begin;
//	int end;
//	int length;
//	int SegCall;
//	char strand[10];
//	bool GENE;
//	char annotation[150];
//	char chrom[10];
//	double median;
//	double mean;
//	double mode;
//	double max;
//	double variance;
//================================================

//================================================
// 	Make Data Structures for Training Data Set
//	Basically a map of maps
class TR_INT{
public:


	int start, end;
	string strand;
	string chrom;


	TR_INT(){};
	TR_INT(int st, int sp){
		start=st;
		end=sp;

	};
	TR_INT(int st, int sp, string str,string chr){
		start	= st;
		end		= sp;
		strand	= str;
		chrom	= chr;
	};

	void SHOW(){
		cout<<"Start: "<< start<< " Stop: "<<end<<endl;
	};
};


class intervals{
public:
	map<string, map<string, vector<TR_INT>> > ON_intervals;
	map<string, map<string, vector<TR_INT>> > OFF_intervals;
	string ReferncedBedGraphFile;
	intervals(map<string, map<string, vector<TR_INT>> > ON, map<string, map<string, vector<TR_INT>> > OFF, string v1 ){
		ON_intervals			= ON;
		OFF_intervals			= OFF;
		ReferncedBedGraphFile	= v1;
	}
	intervals(){};
};


map<string,TR_INT> makeRefStruct(string refGeneFileName){
	ifstream myfile(refGeneFileName);
	string line;
	vector<string> lineInfo;
	map<string,TR_INT> refStruct;
	string strand,chrom;
	if (myfile.is_open()){
		while ( getline(myfile,line) ){
			lineInfo=split(line, "\t");
			strand 	= lineInfo[3];
			chrom	= lineInfo[2];
			refStruct[lineInfo[1]]=TR_INT(atoi(lineInfo[4].c_str()),atoi(lineInfo[5].c_str()),strand, chrom);
			refStruct[lineInfo[12]]=TR_INT(atoi(lineInfo[4].c_str()),atoi(lineInfo[5].c_str()),strand, chrom);
		}
	}
	else{
		cout<<"Unable to open reference Gene File: "<<refGeneFileName<<"...try specifying full path/spelling error correction"<<endl;
	}
	return refStruct;

}



intervals makeIntervalStucture(string Training_FILE)
{
	string pathToSrc = "/"+join(split(__FILE__, "/"), split(__FILE__, "/").size()-2,"/");
	map<string, map<string, vector<TR_INT> > > ON;
	map<string, map<string, vector<TR_INT> >> OFF;
	string line;
	vector<string> lineInfo;
	ifstream TrainingFH(Training_FILE);
	string ReferncedBedGraphFile;
	bool HEADER = 1;
	string strand, chrom;
	int start, stop;
	bool REF=0;
	string GENE;
	TR_INT T;
	map<string,TR_INT> refStruct;
	intervals I;
	if (TrainingFH.is_open()){
		while ( getline(TrainingFH,line) ){

			lineInfo 	= split(line, "\t");
			if (HEADER and lineInfo.size() == 2){
				ReferncedBedGraphFile = lineInfo[1];

				if (ReferncedBedGraphFile == "hg19"){
					string refGeneFileName 			= pathToSrc +"refSeqAnnotations/refSeqGene_hg19.txt";
					refStruct	= makeRefStruct(refGeneFileName);
					REF = 1;
				}


				HEADER = 0;
			}
			else if(HEADER){
				cout<<"Header was not formattted properly in training file..."<<endl;
				cout<<"Exiting..."<<endl;
				return I;
			}
			else if (not HEADER and not REF and lineInfo.size()!=5){
				cout<<"Intervals were not formatted properly in training file..."<<endl;
				cout<<"Offending line: "<<line<<endl;
				cout<<"Exiting..."<<endl;
				return I;
			}

			else if(REF){
				GENE 	= lineInfo[0];
				if (refStruct.find(GENE)!=refStruct.end()){
					T	= refStruct[GENE];
					if (lineInfo[1]=="1"){
						ON[T.strand][T.chrom].push_back(T);
					}else{
						OFF[T.strand][T.chrom].push_back(T);
					}
				}
				else
				{
					cout<<"Couldnt Find: "<<GENE<<endl;
				}
			}

			else{
				strand 	= lineInfo[0];
				chrom	= lineInfo[1];
				start	= atoi(lineInfo[2].c_str());
				stop	= atoi(lineInfo[3].c_str());
				TR_INT	T(start, stop);
				if (lineInfo[4]=="1"){
					ON[strand][chrom].push_back(T);
				}else{
					OFF[strand][chrom].push_back(T);
				}


			}
		}
	}
	else{
		cout<<"Couldn't Open: "<<Training_FILE<<endl;
		cout<<"Exiting..."<<endl;
	}
	TrainingFH.close();
	return intervals(ON, OFF,ReferncedBedGraphFile);

}


vector<vector<double>> standardize(vector<vector<double>> X){
	int N 	= X.size();
	int K 	= X[0].size();
	//calculate the mean / sample variance of each feature dimension
	vector<double> mean(K);
	vector<double> variance(K);
	for (int j = 0; j < K; j++){
		mean[j]		=0;
		variance[j]	= 0;
	}
	for (int i = 0; i < N; i++){
		for (int j = 0; j < K; j++){
			mean[j]+=X[i][j];
		}
	}
	for (int j = 0; j < K; j++){
		mean[j]= mean[j] / N;
	}
	for (int i = 0; i < N; i++){
		for (int j = 0; j < K; j++){
			variance[j]+=pow((X[i][j]-mean[j]), 2);
		}
	}
	for (int i = 0; i < N; i++){
		for (int j = 0; j < K; j++){
			X[i][j]	= (X[i][j]-mean[j]) / sqrt(variance[j]);
		}
	}
	return X;
}

training_data randSample(vector<vector<double>> X, vector<int> Y, int ON_COUNT, int OFF_COUNT, bool SHOW){
	int threshold = 5000;
	double rand;
	uniform_real_distribution<double> distribution(0.0,1.0);
	default_random_engine generator;
	vector<vector<double>> newX;
	vector<int> newY;
	int ON_CT 	= 0;
	int OFF_CT 	= 0;
	int i = 0;
	for (i = 0; i < Y.size(); i++){
		if (distribution(generator) < (OFF_COUNT / (double)(ON_COUNT + OFF_COUNT)) and Y[i] and (ON_CT < threshold)){
				newX.push_back(X[i]);

				newY.push_back(Y[i]);
				ON_CT++;
			}
		if (distribution(generator) < (ON_COUNT / (double)(ON_COUNT + OFF_COUNT)) and (Y[i]==0) and (OFF_CT < threshold)){

			if (X[i][1]<200){
				newX.push_back(X[i]);
				newY.push_back(Y[i]);
			}


			OFF_CT++;
		}
		if (ON_CT == threshold and OFF_CT == threshold){
				break;
		}
	}
//newX 	= standardize(newX);
	if (SHOW){
		cout<<"There are "<<ON_CT<<" ON training examples"<<endl;
		cout<<"There are "<<OFF_CT<<" OFF training examples"<<endl;
	}
	return training_data(newX, newY);

}

training_data fetchTrainingData(string training_FILE, string bedGraphFileName, vector<MEMM> struct_array, string strand, vector<int> dim,string out_dir, string refGeneFileName, bool SHOW){


	intervals Intervals 											= makeIntervalStucture(training_FILE);

	if (SHOW){
		cout<<"Training File interval structure made"<<endl;
	}

	map<string, map<string, vector<TR_INT>> > ON 					= Intervals.ON_intervals;
	map<string, map<string, vector<TR_INT>> > OFF 					= Intervals.OFF_intervals;
	string refbedGraphFileName=Intervals.ReferncedBedGraphFile;

	//bedGraphFileName= /Users/joeyazo/Desktop/Lab/gro_seq/src/files/bed_files/DMSO2_3.pos.BedGraph

	vector<string> linkedFileNames;

	//Search to make sure binary exstis for bedGraphFileName
	string linkedFileName		= "linked_" + split(refbedGraphFileName, "/")[split(refbedGraphFileName, "/").size()-1] + ".dat";
	string fileNamePath 		= join(split(refbedGraphFileName, "/"),split(refbedGraphFileName, "/").size()-1, "/");

	vector<string> paths ;
	paths.push_back(fileNamePath);
	paths.push_back("");
	paths.push_back(out_dir);
	string pathToSrc = "/"+join(split(__FILE__, "/"), split(__FILE__, "/").size()-2,"/");
	paths.push_back(pathToSrc+"segmentation_outputs/");


	//Search for the opposing strand
	string opposingStrandFileName;
	string option1, option2;
	string input_strand;
	if (SEARCH(refbedGraphFileName, "+") || SEARCH(refbedGraphFileName, "pos")){
		if (SEARCH(refbedGraphFileName, "+") ){
			option1 = replace(refbedGraphFileName, "+", "-");
			option2 = replace(refbedGraphFileName, "+", "neg");

		}else{
			option1 = replace(refbedGraphFileName, "pos", "-");
			option2 = replace(refbedGraphFileName, "pos", "neg");
		}

		input_strand='-';
	}
	else if (SEARCH(refbedGraphFileName, "-") || SEARCH(refbedGraphFileName, "neg")){
		if (SEARCH(refbedGraphFileName, "-") ){
			option1 = replace(refbedGraphFileName, "-", "+");
			option2 = replace(refbedGraphFileName, "-", "pos");

		}else{
			option1 = replace(refbedGraphFileName, "neg", "+");
			option2 = replace(refbedGraphFileName, "neg", "pos");
		}
		input_strand='+';
	}






	string linkedOpposingStandName;

	if (not search(paths, option1).empty() and not option1.empty() ){
		linkedOpposingStandName	= "linked_" + split(option1, "/")[split(option1, "/").size()-1] + ".dat";
	}
	else if (not search(paths, option2).empty() and not option2.empty()){
		linkedOpposingStandName	= "linked_" + split(option2, "/")[split(option2, "/").size()-1] + ".dat";
	}
	if (not linkedOpposingStandName.empty() and search(paths, linkedOpposingStandName).empty() and not linkedOpposingStandName.empty()){
		cout<<"Linked File: "<<linkedOpposingStandName<<" doesn't exist...making it"<<endl;
		linkFile(bedGraphFileName, input_strand,  refGeneFileName, out_dir+linkedOpposingStandName,0,SHOW);
	}






	typedef map<string, map<string, vector<TR_INT>> > it_type;
	typedef map<string, vector<TR_INT>> it_type2;
	int ON_CT=0;
	int OFF_CT=0;
	for (it_type::iterator it = ON.begin(); it!=ON.end(); it++){
		for (it_type2::iterator it2 = ON[it->first].begin(); it2!=ON[it->first].end(); it2++){
			ON_CT+=ON[it->first][it2->first].size();



		}
	}
	for (it_type::iterator it = OFF.begin(); it!=OFF.end(); it++){
		for (it_type2::iterator it2 = OFF[it->first].begin(); it2!=OFF[it->first].end(); it2++){
			OFF_CT+=OFF[it->first][it2->first].size();
		}
	}
	if (not linkedOpposingStandName.empty()){
		linkedFileNames.push_back(search(paths, linkedOpposingStandName));
	}



	if (search(paths, linkedFileName).empty()){
		cout<<"Couldnt find: "<<linkedFileName<<endl;
		cout<<"This will error..."<<endl;
	}

	linkedFileNames.push_back(search(paths, linkedFileName));



	if (SHOW){
		cout<<"...Made Interval Data Structure..."<<endl;
		cout<<"...Number of Contigs: "<<struct_array.size()<<endl;
		cout<<"Number of ON intervals: "<<ON_CT<<endl;
		cout<<"Number of OFF intervals: "<<OFF_CT<<endl;

	}
	string ReferencedBedGraphFile 									= Intervals.ReferncedBedGraphFile;
	string ReferenceStrand 											= "";






	//=======================================
	//	Need to make On and OFF interval
	//=======================================

	vector<vector<double>> 	X;
	vector<int> 			Y;
	int start,end;
	int chrom_ct 	= 0;
	int ON_COUNT 	= 0;
	int OFF_COUNT 	= 0;
	int N			= struct_array.size();
	string curr_chrom, curr_strand;
	vector<double> features;
	vector<double> featuresPast;
	vector<double> featuresFuture;
	string prev_chrom	= "";
	int CHROM_TOTAL		= ON["-"].size() + ON["+"].size();
	for (vector<string>::iterator file =linkedFileNames.begin(); file!=linkedFileNames.end(); ++file ){
		ifstream binaryFile(*file, ios::binary|ios::in);
		binaryFile.seekg(0, ios::beg);
		MEMM memmObj;
		cout<<*file<<endl;
		if (binaryFile.is_open()){
			int j			= 0;
			int k			= 0;
			int chrom_ct 	= 0;
			while (1){
				binaryFile.read((char*)&memmObj, sizeof(MEMM));
				start 			= memmObj.begin;
				end				= memmObj.end;
				curr_strand 	= string(memmObj.strand);
				curr_chrom		= string(memmObj.chrom);
				if (curr_strand.empty()){
					curr_strand = "+";
				}

				if (binaryFile.eof() or ON[curr_strand].empty()){
					break;
				}
				while ((j<ON[curr_strand][curr_chrom].size()) && (ON[curr_strand][curr_chrom][j].end < start)){
					j++;
				}
				while (k<OFF[curr_strand][curr_chrom].size() && OFF[curr_strand][curr_chrom][k].end < start){
					k++;
				}
				if (j < ON[curr_strand][curr_chrom].size() and ON[curr_strand][curr_chrom][j].start < start
							and end < ON[curr_strand][curr_chrom][j].end)
				{
										features = memmObj.getFeatures(dim);
					X.push_back(features);
					Y.push_back(1);
					ON_COUNT++;
				}
				if (k < OFF[curr_strand][curr_chrom].size() and
						 OFF[curr_strand][curr_chrom][k].start < start and end < OFF[curr_strand][curr_chrom][k].end){
					features = memmObj.getFeatures(dim);
					X.push_back(features);
					Y.push_back(0);
					OFF_COUNT++;
				}
				if (prev_chrom!=curr_chrom){
					if (j or k){
						chrom_ct++;
					}
					if (chrom_ct==CHROM_TOTAL){
						break;
					}
					j=0;
					k=0;
				}
				prev_chrom = curr_chrom;
			}
			binaryFile.close();

		}

		else{
			cout<<"Couldn't Open Binary File: "<<*file<<endl;
		}
	}

	//return training_data(X,Y);
	return randSample(X,Y, ON_COUNT, OFF_COUNT,SHOW);
}



#endif /* GRABTRAININGDATA_H_ */
