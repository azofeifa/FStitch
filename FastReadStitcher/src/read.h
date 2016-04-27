#ifndef read_H
#define read_H
#include <string>
#include <vector>
#include <map>
using namespace std;
class contig{
public:
	int start, stop;
	double left, right, length;
	float cov;
	string chrom;
	contig * next;
	contig();
	void setStats(int , int , double , double , double, float , string);
	void display();
	vector<double> getVect(bool);
};
class interval{
public:
	int start;
	int stop;
	string info;
	interval * next = NULL;
	interval(int, int, string);
	interval();	
};
class RTOF{
public:
	vector<double> W;
	vector<vector<double>> A;
	bool ChIP;
	bool EXIT;
	RTOF(vector<double>,vector<vector<double>>, bool);
	RTOF();

};

class readTrainingFileReturn{
public:
	bool EXIT;
	map<string, interval *> result;
	readTrainingFileReturn();
	
};

class segment{
public:
	string chrom; 
	int start, stop, ID, chrom_ID;
	double minX, maxX;
	vector< vector<double> > forward;
	vector< vector<double> > reverse;
	string strand ;
	int counts;
	vector<double> centers;
	vector<vector<double>> parameters; //for bootstrapping
	map<int, vector<double> > variances;
	segment(string, int , int);
	segment(string, int , int, int);
	segment(string, int , int, int,string);

	segment();
	string write_out();
	void bin(double, double, bool);
	void add2(int , vector<double> );
	vector<vector<double>> get_contig_info(int ,vector<vector<double>> &);
	double N;
	double fN;
	double rN;
	double XN;
	double ** X;
	double SCALE;
	vector<vector<double> > bidirectional_bounds;
	vector<segment *> bidirectional_data;
	vector<int>  bidir_counts; //used for optimization of BIC?
	vector<int> bidirectional_N;
	vector<vector<double>> fitted_bidirs; //mu, si, l,pi
};

class node{
public:
	double center;
	int start, stop;
	node * left;
	node * right;
	vector<segment * > current;
	void retrieve_nodes(vector<segment * >&);
	void insert_coverage(vector<double>, int);
	node();
	node(vector<segment *>);
	void searchInterval(int, int, vector<int> &) ;
};

readTrainingFileReturn readTrainingFile(string);
map<string,contig *> readBedGraphFile(string, map<string, interval *>, bool);
map<string, map<string, interval *>> readRefSeq(string);
RTOF readTrainingOutFile(string);
map<string,contig *> readBedGraphFileAll(string,int);
class run_out{
public:
	vector< vector<double> > X ;
	vector<int> Y;
	bool EXIT;
	run_out(vector< vector<double> >, vector<int> );
	run_out();
};



namespace load{
	map<string, segment*> load_bedgraphs_total(string, 
		string, string);
	vector<segment *> load_intervals_of_interest(string);
	vector<segment* > insert_bedgraph_to_segment_joint(map<string, vector<segment *> >  , 
		 string );
	map<string, vector<segment *> > convert_segment_vector(vector<segment *> );

	run_out convert_to_run_out(vector<segment *>);


}

#endif