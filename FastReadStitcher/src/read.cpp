#include "read.h"
#include <string>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <map>
#include <cmath>
#include "BaumWelch.h"
#include "split.h"
#include "validate.h"
#include <stdexcept>
using namespace std;

string getChrom(string line){
	const char * tab = "\t";
	for (int i = 0; i < line.size(); i++){
		if (line[i]==*tab){
			return line.substr(0,i);
		}
	}
	return "";
}
class file_stats{
public:
	vector<int> start_stop;
	map<int, string> relate;
	file_stats(vector<int> Start_stop, map<int, string> Relate){
		start_stop=Start_stop;
		relate=Relate;
	}
};
file_stats getFHstats(string FILE){
	ifstream FH(FILE);
	vector<int> start_stop;
	map<int, string> relate;
	string line;
	string chrom;
	string prev="";
	int i=0; 
	if (FH){
		while (getline(FH,line)){
			chrom = getChrom(line);
			if (prev!=chrom){
				start_stop.push_back(FH.tellg());
				relate[i] 	= chrom;
				i++;
			}
			prev=chrom;
		}

	}else{
		cout<<"\""<<FILE<<"\""<<" doesn't exist, exiting..."<<endl;
	}
	start_stop.push_back(FH.tellg());
	FH.close();
	return file_stats(start_stop, relate);
}
contig::contig(){}
void contig::setStats(int st, int sp, double l, double r, double len, float C, string CHROM){
	start 	= st;
	stop 	= sp;
	left 	= l;
	right 	= r;
	cov 	= C ;
	length 	= len;
	chrom 	= CHROM;
}
void contig::display(){
	cout<<chrom<<":"<<start<<"-"<<stop<<endl;
}
vector<double> contig::getVect(bool ChIP){
	vector<double> x;
	if (not ChIP){
		x.push_back(1);
		x.push_back(log((left+right)/2));
		x.push_back(log(length));
		x.push_back(log(cov/ (left+right+length)));
	}else{
		x.push_back(1);
		x.push_back(cov / (length+ ((left+right)/2) ));
		x.push_back(cov);
	}
	return x;
}
class contigOut{
public:
	bool EXIT;
	contig * result;
};


contigOut makeContig(string FILE, int start, int stop){
	contigOut CO;
	ifstream FH(FILE);
	FH.seekg(start);
	string line, chrom;
	vector<string> lineArray;
	vector<contig> contigs;
	int st, sp;
	float cov;
	int prevStart= 0;
	int prevStop = 0;
	int p = 0;
	double l = 0;
	double r = 0;
	bool begin=1;
	contig * C;
	C 		= new contig;
	contig * root 	= C;
	float coverage 	= 0;

	while (FH.tellg()<stop){
		getline(FH,line);
		lineArray 	= splitter(line, "\t");
		if (lineArray.size()!=4){
			cout<<endl;
			cout<<"couldn't parse line: "+line+"\n";
			CO.EXIT=true;
			return CO;
		}else{

			chrom 	= lineArray[0];
			if (not isNum(lineArray[1]) or not  isNum(lineArray[2]) or not isNum(lineArray[3]) ){
				cout<<endl;
				cout<<"Line: "<<line<<endl;
				cout<<"could not convert coordinates or coverage value to number"<<endl;
				CO.EXIT=true;
				return CO;
			}
			sp 		= stoi(lineArray[2]);
			st 		= stoi(lineArray[1]);
			cov 	= stof(lineArray[3]);
			if (begin){
				begin 		= 0;
				prevStop	= sp;
				l 			= sp-prevStop;
				prevStart 	= 0;
			}else if(not begin and (st-p)>2 ){
				r 			= st - p;
				C->setStats(prevStart-l,p+r, l, r, p-prevStart, coverage, chrom);
				C->next 	= new contig;
				C 			= C->next;
				prevStart 	= st;
				l 			= prevStart - p;

				coverage=0;			
			}

			p 		= sp;
			coverage+=cov;
		}
		start++;
	}
	C->next 	= NULL;
	
	FH.close();
	CO.EXIT 	= false;
	CO.result 	= root;
	return CO;
}

interval::interval(int st, int sp , string INFO){
	start 	= st;
	stop 	= sp;
	info 	= INFO;
	next 	= NULL;
}
interval::interval(){}

readTrainingFileReturn::readTrainingFileReturn(){}

readTrainingFileReturn readTrainingFile(string FILE){
	readTrainingFileReturn RETURN;
	map<string, interval *> 	R;
	map<string, interval *> 	roots;
	
	vector<int> start_stop(3);
	ifstream FH(FILE);
	string line;
	vector<string>lineArray;
	if (FH){
		while (getline(FH, line)){
			lineArray 			= splitter(line, "\t");
			if (lineArray.size()!=4){
				cout<<"Line: "<<line<<", in training file is not formatted properly, tab delimited"<<endl;
				RETURN.EXIT 	= true;
				return RETURN;
			}

			if (not (lineArray[3]=="0" or lineArray[3]=="1")){
				cout<<endl;
				cout<<"Line: "<<line<<", must have training label as either 0 or 1"<<endl;
				RETURN.EXIT 	= true;
				return RETURN;
			}
			if (not isNum(lineArray[1]) or not isNum(lineArray[2])){
				cout<<"Line: "<<line<<", coordinates must be numbers"<<endl;
				RETURN.EXIT 	= true;
				return RETURN;
			}
			if (R.find(lineArray[0])==R.end()){
				R[lineArray[0]] 		= new interval(stoi(lineArray[1] ), stoi(lineArray[2]), lineArray[3]);
				roots[lineArray[0]] 	= R[lineArray[0]];
			}else{
				R[lineArray[0]]->next 	= new interval(stoi(lineArray[1] ), stoi(lineArray[2]), lineArray[3]);
				R[lineArray[0]] 		= R[lineArray[0]]->next;
			}

		}
	}else{
		cout<<"couldn't open: "<<FILE<<endl;
		cout<<"exiting..."<<endl;
	}
	typedef map<string, interval *>::iterator r_it;
	for (r_it I = R.begin(); I!=R.end(); I++){
		R[I->first] 	= roots[I->first];
	}

	FH.close();
	RETURN.EXIT 	= false;
	RETURN.result 	= R;
	return RETURN;

}
map<string,contig *> readBedGraphFile(string FILE, map<string, interval *> T, bool verbose){
	file_stats fs 					= getFHstats(FILE);
	vector<int> start_stop 			= fs.start_stop;

	map<int,string> relate 			= fs.relate;
	map<string,contig *> 	D;
	vector<contig *> M(start_stop.size());
	int nthreads 					= omp_get_max_threads();
	bool abort = false;
	#pragma omp parallel for
	for(int n=1; n<start_stop.size(); ++n)
	{
		#pragma omp flush (abort)
		if (T.find(relate[n-1])!=T.end()){
			contigOut CO 	= makeContig(FILE, start_stop[n-1],start_stop[n]);
			if (CO.EXIT){
				D.clear();
				abort = true;
			}
			M[n-1] 			= CO.result;
			D[relate[n-1]] 	= M[n-1];
		}
		
	}
	if (abort){
		D.clear();
	}
	return D;
}
map<string,contig *> readBedGraphFileAll(string FILE,int np){
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(np); // Use 4 threads for all consecutive parallel regions
	
	file_stats fs 					= getFHstats(FILE);
	vector<int> start_stop 			= fs.start_stop;
	map<int,string> relate 			= fs.relate;
	map<string,contig *> 	D;
	vector<contig *> M(start_stop.size());
	bool abort = false;
	#pragma omp parallel for
	for(int n=1; n<start_stop.size(); ++n)
	{
		#pragma omp flush (abort)
		contigOut CO 	= makeContig(FILE, start_stop[n-1],start_stop[n]);
		if (CO.EXIT){
			D.clear();
			abort = true;
		}
		M[n-1] 			= CO.result;
		D[relate[n-1]] 	= M[n-1];
		
	}
	if (abort){
		D.clear();
	}
	return D;
}


map<string, map<string, interval *>> readRefSeq(string FILE){
	map<string, map<string, interval *>> R;	
	map<string, map<string, interval *>> roots;	
	
	R["+"] =	map<string, interval *>();
	R["-"] = 	map<string, interval *>();
	roots["+"] =	map<string, interval *>();
	roots["-"] = 	map<string, interval *>();
	
	if (not FILE.empty()){
		ifstream FH(FILE);
		string line, chrom, strand, ID;
		int 	start, stop; 
		vector<string> lineArray;
		if (FH){
			while (getline(FH,line)){
				lineArray 	= splitter(line, "\t");
				chrom 		= lineArray[2], strand = lineArray[3], ID=lineArray[1];
				start 		= stoi(lineArray[4]), stop=stoi(lineArray[5]);
				if (R[strand].find(chrom)==R[strand].end()){
					R[strand][chrom] 		= new interval(start, stop, ID);
					roots[strand][chrom] 	= R[strand][chrom];
				}else{
					R[strand][chrom]->next 	= new interval(start, stop, ID);
					R[strand][chrom] 		= R[strand][chrom]->next;
				}
			}
		}else{
			cout<<"couldn't open: "<<FILE<<"\nExiting..."<<endl;
			return R;
		}
		FH.close();
	}
	typedef map<string, map<string, interval *>>::iterator st_it;
	typedef map<string, interval *>::iterator chrom_it;
	
	for (st_it strand = R.begin(); strand != R.end(); strand++ ){
		for (chrom_it chrom = strand->second.begin(); chrom!= strand->second.end(); chrom++){
			R[strand->first][chrom->first] 	= roots[strand->first][chrom->first];
		}
	}
	
	return R;
}
RTOF::RTOF(vector<double> w, vector<vector<double>> a, bool CH){
		W=w,A=a;
		EXIT=false;
		ChIP=CH;
}
RTOF::RTOF(){
	EXIT=true;
}

RTOF readTrainingOutFile(string FILE){
	ifstream FH(FILE);
	string line;
	vector<double> W;
	vector<vector<double>> A;
	vector<string> lineArray;
	bool begin 	= 1;
	bool ChIP 	= 0;
	if (FH){
		while (getline(FH,line)){
			if (begin and ("#" != line.substr(0,1) ) ){
				RTOF ROOT;
				cout<<"This is not an output training file\nfrom the fast read stitcher"<<endl;
				return ROOT;
			}

			begin = false;
			if ("#ChIP" == line.substr(0,5)){
				lineArray 		= splitter(line, ":");
				ChIP 			= (lineArray[1]=="1");
			}
			if ("#" != line.substr(0,1)){
				lineArray 		= splitter(line, ":");
				if (lineArray[0].substr(0,1)=="L"){
					lineArray 	= splitter(lineArray[1], ",");
					for (int i = 0; i < lineArray.size(); i++){
						W.push_back(stof(lineArray[i]));
					}
				}else if(lineArray[0].substr(0,1)=="H"){
					lineArray 	= splitter(lineArray[1], ",");
					int k 		= 0;
					for (int i  = 0; i < 2;i++){
						vector<double> row;
						for (int j=0;j<2;j++){
							row.push_back(stof(lineArray[k]));
							k++;
						}
						A.push_back(row);
					}
				}
			}
		}

	}else{
		cout<<"\""<<FILE<<"\""<<" doesn't exist, exiting..."<<endl;
	}
	FH.close();
	return RTOF(W,A, ChIP);
}




