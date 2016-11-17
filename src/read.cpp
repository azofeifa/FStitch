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
	int LAST 	= 0;
	if (FH){
		while (getline(FH,line)){
			if (!line.empty() && line[line.size() - 1] == '\r'){
				line.erase(line.size() - 1);
			}
			
			chrom = getChrom(line);
			if (prev!=chrom){
				start_stop.push_back(FH.tellg());
				relate[i] 	= chrom;
				i++;
			}
			prev=chrom;
			LAST 	= FH.tellg();
		}

	}else{
		cout<<"\""<<FILE<<"\""<<" doesn't exist, exiting..."<<endl;
	}
	start_stop.push_back(LAST);
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
		if (!line.empty() && line[line.size() - 1] == '\r'){
			line.erase(line.size() - 1);
		}
			
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
				prevStart 	= st;
			}else if(not begin and (st-p)>2 ){
				r 			= st - p;
				C->setStats(prevStart,p, l, r, p-prevStart, coverage, chrom);
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
			if (!line.empty() && line[line.size() - 1] == '\r'){
				line.erase(line.size() - 1);
			}
			lineArray 			= splitter2(line, "\t");
			if (lineArray.size()!=4){
				cout<<"Line: "<<line<<", in training file is not formatted properly, tab delimited"<<endl;
				RETURN.EXIT 	= true;
				return RETURN;
			}
			if (not (lineArray[3]=="0" or lineArray[3]=="1")){
				printf("Line: %s, must contain either 0 or 1 as training input\n",line.c_str() );
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


//========================================================================
//The very very important segment class

segment::segment(string chr, int st, int sp){
	chrom	= chr;
	start	= st;
	stop	= sp;
	N 		= 0;
	fN 		= 0;
	rN 		= 0;
	minX=st, maxX=sp;
	counts 	= 1;
	XN 		= 0;
	ID 		= 0;
	strand 	= ".";
	chrom_ID= 0;

}
segment::segment(string chr, int st, int sp, int i){
	chrom	= chr;
	start	= st;
	stop	= sp;
	fN 		= 0;
	rN 		= 0;
	N 		= 0;
	minX=st, maxX=sp;
	counts 	= 1;
	XN 		= 0;
	ID 		= i;
	strand 	= ".";
	chrom_ID= 0;

}

segment::segment(string chr, int st, int sp, int i, string STR){
	chrom	= chr;
	start	= st;
	stop	= sp;
	N 		= 0;
	fN 		= 0;
	rN 		= 0;
	minX=st, maxX=sp;
	counts 	= 0;
	XN 		= 0;
	ID 		= i;
	strand 	= STR;
	chrom_ID= 0;

}

segment::segment(){
	N 		= 0;
	fN 		= 0;
	rN 		= 0;
	counts 	= 1;
	XN 		= 0;
	ID 		= 0;
	strand 	= ".";
	chrom_ID= 0;
}

string segment::write_out(){
	string text 	= ("#" + chrom + ":" + to_string(start) + "-" 
		+ to_string(stop) + "," + to_string(int(N))+ "\n");
	return text;
}


void segment::add2(int strand, vector<double> x){
	if (strand == 1){
		forward.push_back(x);
	}else if (strand==-1){
		reverse.push_back(x);
	}
}



void segment::bin(double delta, double scale, bool erase){
	X 				= new double*[3];
	SCALE 			= scale;
	int BINS;
	BINS 		= (maxX-minX)/delta;
	start = minX, stop=maxX;
	for (int j = 0 ; j < 3;j++){
		X[j] 		= new double[BINS];
	}
	N 				= 0;
	fN = 0, rN = 0;
	XN 				= BINS;
	//===================
	//populate bin ranges
	X[0][0] 		= double(minX);
	X[1][0]=0,X[2][0]=0;
	



	for (int i = 1; i < BINS; i++){
		X[0][i] 	= X[0][i-1] + delta;
		X[1][i] 	= 0;
		X[2][i] 	= 0;
	}
	// ===================
	//insert forward strand
	int j 	=0;
	//printf("start: %d , stop: %d , bins: %d ,delta: %f, forward: %d, reverse: %d\n", start, stop, BINS, delta, forward.size(), reverse.size() );
	for (int i = 0 ; i < forward.size(); i++){
		while (j < BINS and X[0][j] <=forward[i][0]){
			j++;
		}
		if (j < BINS and forward[i][0]<= X[0][j]){
			X[1][j-1]+=forward[i][1];
			N+=forward[i][1];
			fN+=forward[i][1];
		}
	}
	j 	=0;
	//===================
	//insert reverse strand
	for (int i = 0 ; i < reverse.size(); i++){
		while (j < BINS and X[0][j] <=reverse[i][0]){
			j++;
		}
		if (j < BINS and reverse[i][0]<= X[0][j]){
			X[2][j-1]+=reverse[i][1];
			N+=reverse[i][1];
			rN+=reverse[i][1];
		}
	}
	//===================
	//scale data down for numerical stability
	if (scale){
		for (int i = 0; i < BINS; i ++ ){

			X[0][i] 	= (X[0][i]-minX)/scale;
			// X[1][i]/=delta;
			// X[2][i]/=delta;
		}
	}
	//we also want to get rid of those data points that we don't need
	//i.e. the ones where there is no data coverage values on either the 
	//forward or reverse strands

	int realN 		= 0;
	for (int i = 0; i < BINS;i++){
		if (X[1][i]>0 or X[2][i]>0){
			realN++;
		}
	}
	if (erase){
		double ** newX 	= new double*[3];
		for (int j=0; j<3;j++){
			newX[j] 	= new double[realN];
		}
		j = 0;
		for (int i = 0; i < BINS; i ++){
			if (X[1][i]>0 or X[2][i]>0){
				newX[0][j] 	= X[0][i];
				newX[1][j] 	= X[1][i];
				newX[2][j] 	= X[2][i];
				j++;
			}
		}
		if (realN!=j){
			printf("WHAT? %d,%d\n", j, realN);
		}
		//clear previous memory
		for (int i = 0; i < 3; i ++){
			delete X[i];
		}
		delete X;
		X 				= newX;
		XN 				= realN;
	}
	if (scale){
		if (not centers.empty()){
			for (int i = 0; i < centers.size(); i++){
				centers[i]=(centers[i]-minX)/scale;			
			}
		}
		if (not fitted_bidirs.empty() ){
			for (int fb = 0; fb < fitted_bidirs.size(); fb++){
				int center 	= fitted_bidirs[fb][0];
				int std 	= fitted_bidirs[fb][1]*0.5 + (1. /  fitted_bidirs[fb][2]);
				int a 		= center - std*3;
				int b 		= center + std*3;
				


				fitted_bidirs[fb][0] = (fitted_bidirs[fb][0] - minX)/scale;
				fitted_bidirs[fb][1] /= scale;
				fitted_bidirs[fb][2] *= scale; 
			}
		}

		maxX 			= (maxX-minX)/scale;
		minX 			=0;
	}
	double S=0;
	for (int i = 0; i < XN; i++){
		S+=X[1][i];
	}
	forward.clear();
	reverse.clear();
}

vector<vector<double>> segment::get_contig_info(int strand, 
		vector<vector<double>> & st_sp){
	vector<vector<double>> X;
	double prev 	= 0;
	double st 		= 0;
	double S 		= 0, M=0;
	double sp 		= 0;
	double l,gap;
	vector<vector<double>> current;
	if (strand == 1){
		current 	= forward;
	}else{
		current 	= reverse;
	}
	for (int i = 0 ; i < current.size(); i++){
		if ( abs(current[i][1] - prev) > 1.0){
			if (st != 0){
				sp 	= prev;
				l 		= sp-st, gap 	= prev-current[i][1];
				vector<double> x 		= {1, gap, l, S/l,M, sqrt(pow(l,2) + pow(gap,2))  };
				X.push_back(x);
				vector<double> y 		= {st, current[i][1]};
				st_sp.push_back(y);
			}
			st 		= current[i][1];
			S 		= 0;
			M 		= 0;
		}
		prev 	= current[i][2];
		S  		+= current[i][3];
		M 		= max(current[i][3],M);
	}
	return X;
}




//================================================================================================
//interval tree code

node::node(){};

node::node(vector<segment * > segments ){
	center 	= (double(segments[0]->start)  + double(segments[segments.size()-1]->stop)) / 2.;
	vector<segment * > Left;
	vector<segment * > Right;
	left=NULL, right=NULL;
	for (int i = 0 ; i < segments.size(); i++){
		if (segments[i]->stop < center){
			Left.push_back(segments[i]);
		}
		else if (segments[i]->start > center){
			Right.push_back(segments[i]);
		}
		else{
			current.push_back(segments[i]);
		}
	}
	if (Left.size() > 0){
		left 	= new node(Left);
	}
	if (Right.size() > 0){
		right 	= new node(Right);
	}
}
void node::insert_coverage(vector<double> x, int s){
	for (int i = 0 ; i < current.size(); i++){
		if (x[0] > current[i]->start and  x[0] < current[i]->stop  ){
			if (s==1){
				current[i]->forward.push_back(x);
			}else{
				current[i]->reverse.push_back(x);	
			}
		}
	}	

	if (x[0] >= center and right != NULL ){
		right->insert_coverage(x, s);
	}
	if (x[0] <= center and left !=NULL){
		left->insert_coverage(x,  s);
	}
}
void node::searchInterval(int start, int stop, vector<int>& finds ){
	for (int i = 0 ; i < current.size(); i++){
		if (stop > current[i]->start and  start < current[i]->stop  ){
			finds.push_back(1);
		}
	}	
	if (start >= center and right != NULL ){
		right->searchInterval(start, stop, finds);
	}
	if (stop <= center and left !=NULL){
		left->searchInterval(start, stop, finds);
	}	
}

void node::retrieve_nodes(vector<segment*> & saves){
	for (int i = 0; i < current.size(); i++){
		saves.push_back(current[i]);
	}
	if (right!= NULL){
		right->retrieve_nodes(saves);
	}
	if (left != NULL){
		left->retrieve_nodes(saves);		
	}
}


//================================================================================================
//LOADING from file functions...need to clean this up...


map<string, segment*> load::load_bedgraphs_total(string forward_strand, 
	string reverse_strand, string joint_bedgraph){

	map<string, segment*> 	G;
	vector<segment*> segments;
	vector<string> FILES;
	if (forward_strand.empty() and reverse_strand.empty()){
		FILES 	= {joint_bedgraph};
	}else if (not forward_strand.empty() and not reverse_strand.empty()){
		FILES 	= {forward_strand, reverse_strand};
	}else{
		return G;
	}
	
	string line, chrom;
	int start, stop;
	double coverage;
	vector<string> lineArray;
	string prevChrom="";
	segment * S =NULL;
	bool INSERT 	= false;
	bool EXIT 		= false;
	int line_number = 0;
	for (int u = 0 ; u < FILES.size(); u++){
		ifstream FH(FILES[u]) ;
		if (not FH ){
			printf("couln't open FILE %s\n", FILES[u].c_str());
		}
		if (EXIT){
			break;
		}
		while (getline(FH, line)){
			lineArray=splitter(line, "\t");
			if (lineArray.size()!=4){
				EXIT 	= true;
				printf("\nLine number %d  in file %s was not formatted properly\nPlease see manual\n",line_number, FILES[u].c_str() );
			}
			line_number++;
			chrom=lineArray[0], start=stoi(lineArray[1]), stop=stoi(lineArray[2]), coverage=(stof(lineArray[3]));
			if (chrom != prevChrom   )  {
				if (chrom.size() < 6 and u==0){
					G[chrom] 	= new segment(chrom, start, stop );
					INSERT 		= true;
				}else if(chrom.size() > 6){
					INSERT 		= false;
				}
			}
			if (INSERT){
				vector<double> x(4);
				double center 	= (stop + start) /2.;

				x[0]=center, x[1]=start, x[2]=stop, x[3] = abs(coverage);

				if (u==0){
					if (coverage > 0){
						G[chrom]->add2(1, x);
					}else{
						G[chrom]->add2(-1, x);						
					}
				}else{
					G[chrom]->add2(-1, x);	
				}
			}
			prevChrom=chrom;

		}
	}
	
	return G;
}
vector<segment*> load::load_intervals_of_interest(string FILE){
	ifstream FH(FILE);


	vector<segment *> G;
	int ct 	= 1;
	map<string, vector<segment * > > GS;
	map<int, string> IDS_first;
	int T 	= 0;
	bool EXIT 		= false;
	if (FH){
		string line, chrom;
		int start, stop, label;
		int 	i = 0;
		vector<string> lineArray;
		string strand; 
		bool PASSED 	= true;

		while(getline(FH, line)){
			lineArray=splitter(line, "\t");
			if (lineArray[0].substr(0,1)!="#" and lineArray.size()>2){
				try{
					chrom=lineArray[0], start=max(stoi(lineArray[1]), 0), stop=stoi(lineArray[2]);
					label 	= stoi(lineArray[3]);
					
				}
				catch(exception& e){
					printf("\n\nIssue with file %s at line %d\nPlease consult manual on file format\n\n",FILE.c_str(), i );
					EXIT=true;
					GS.clear();
					break;
				}
				if (start < stop){
					
					segment * S 	= new segment(chrom, start, stop,label,strand);
					GS[S->chrom].push_back(S);
				
					i++;
				}
			}
		}
	}else{
		printf("couldn't open %s for reading\n", FILE.c_str() );
		EXIT 	= true;
	}
	if (not EXIT){
		typedef map<string, vector<segment * > >::iterator it_type;
		for (it_type c 	= GS.begin(); c!=GS.end(); c++){
			vector<segment *> m_segs;
			m_segs 	= c->second;	
			for (int i = 0 ; i < m_segs.size(); i++){
				G.push_back(m_segs[i]);
			}
		}
	}else{
		G.clear();
	}
	return G;
}

vector<segment* > load::insert_bedgraph_to_segment_joint(map<string, vector<segment *> > A , 
	string bedgraph ){
	
	
	
	map<string, node> NT;
	typedef map<string, vector<segment *> >::iterator it_type_5;
	for(it_type_5 c = A.begin(); c != A.end(); c++) {
		NT[c->first] 	= node(c->second);
	}
	int start, stop, N, j;
	double coverage;
	N 	= 0,j 	= 0;
	int strand;
	int o_st, o_sp;
	vector<string> lineArray;
	string chrom, prevchrom, line;
	vector<segment *> segments;
	double center;
	vector<string> FILES;
	FILES 	= {bedgraph};
	string FILE;
	for (int i =0; i < FILES.size(); i++){
		FILE=FILES[i];
		ifstream FH(FILE);
		if (FH){
			prevchrom="";
			while (getline(FH, line)){
				lineArray 	= splitter2(line, "\t");
				if (lineArray.size()==4){
					chrom 		= lineArray[0];
					start=stoi(lineArray[1]),stop=stoi(lineArray[2]), coverage = stod(lineArray[3]);
					
					center 	= (stop + start) /2.;
					if (NT.find(chrom)!=NT.end()){
						vector<double> x(4);
						x[0]=center, x[1]=start, x[2]=stop, x[3] = abs(coverage);
						NT[chrom].insert_coverage(x, 1);
						
					}
				}
				else{
					printf("\n***error in line: %s, not bedgraph formatted\n", line.c_str() );
					segments.clear();
					return segments;
				}
			}
			FH.close();
			
		}else{
			cout<<"could not open forward bedgraph file: "<<FILE<<endl;
			segments.clear();
			return segments;
		}
	}
	//now we want to get all the intervals and make a vector<segment *> again...
	vector<segment *>NS;
	typedef map<string, node>::iterator it_type_6;
	for (it_type_6 c = NT.begin(); c!=NT.end(); c++){
		c->second.retrieve_nodes(NS);
	}

	return NS;
}
map<string, vector<segment *> > load::convert_segment_vector(vector<segment *> FSI){
	map<string, vector<segment *>> GG;
	for (int s = 0 ; s < FSI.size(); s++){
		GG[FSI[s]->chrom].push_back(FSI[s]);
	}
	return GG;
}
run_out::run_out(vector< vector<double> >x , vector<int> y){
	X=x,Y=y;
};
run_out::run_out(){};
run_out load::convert_to_run_out(vector<segment *> FSI){
	vector<vector<double>> X;
	vector<int> Y;
	for (int s = 0 ; s < FSI.size(); s++){
		vector<vector<double>> xy;
		vector<vector<double>> d 	= FSI[s]->get_contig_info(1, xy);
		for (int i = 0 ; i < d.size(); i++){
			X.push_back(d[i]);
			Y.push_back(FSI[s]->ID);
		}
	}
	run_out RO(X,Y);
	return RO;
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
		if (T.find(relate[n-1])!=T.end() ){
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
				if (!line.empty() && line[line.size() - 1] == '\r'){
					line.erase(line.size() - 1);
				}
			
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
			if (!line.empty() && line[line.size() - 1] == '\r'){
				line.erase(line.size() - 1);
			}
			
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




