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


readTrainingFileReturn readTrainingFile(string);
map<string,contig *> readBedGraphFile(string, map<string, interval *>, bool);
map<string, map<string, interval *>> readRefSeq(string);
RTOF readTrainingOutFile(string);
map<string,contig *> readBedGraphFileAll(string,int);
#endif