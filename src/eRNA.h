/*
 * eRNA.h
 *
 *  Created on: May 7, 2014
 *      Author: joeyazo
 */

#ifndef ERNA_H_
#define ERNA_H_
class eInterval{
public:
	int start, stop;
	double score;
	string chrom;
	eInterval(int st, int sp, double sc){
		start	= st;
		stop	= sp;
		score	= sc;
	}
	eInterval(int st, int sp, double sc, string chr){
		start	= st;
		stop	= sp;
		score	= sc;
		chrom	= chr;
	}


};


map<string, vector<eInterval>> makeIntervals(string file){
	map<string, vector<eInterval>> intervals;
	ifstream myfile(file);
	string line;
	vector<string> lineInfo;
	bool header=1;
	string chrom, state;
	int start, stop;
	double score;
	if (myfile.is_open()){
		while ( getline(myfile,line) ){
			if (not header){
				lineInfo=split(line,"\t");
				if (lineInfo.size() <4){
					cout<<"File not formatted properly..."<<endl;
					cout<<"Are you sure this file was the output of LogitGRO?"<<endl;
					break;
				}
				chrom	= lineInfo[0];
				start 	= atoi(lineInfo[1].c_str());
				stop 	= atoi(lineInfo[2].c_str());
				if (split(lineInfo[3], "=").size() >1){
					state	= split(lineInfo[3], "=")[0];
					score	= atof(split(lineInfo[3], "=")[1].c_str());
					if (state == "ON"){
						intervals[chrom].push_back(eInterval(start, stop,score));
					}
				}

			}else{
				header=0;
			}

		}
	}
	else{
		cout<<"file: "<<file<<" doesn't exist..."<<endl;
		cout<<"please specify full path..."<<endl;
	}

	return intervals;
}

int min(int A, int B){
	if (A < B){
		return A;
	}
	return B;
}
int max(int A, int B){
	if (A > B){
		return A;
	}
	return B;
}

class normal{
public:
	double mean, variance;
	normal(double m, double v){
		mean=m;
		variance=v;
	}
	double cdf(double x){
		return 0.5 * (1.0 + erf((x - mean) / sqrt(2*variance)));
	}


};

vector<eInterval> getOverlaps(map<string, vector<eInterval>> A_I, map<string, vector<eInterval>> B_I, double thresh){
	vector<eInterval> overlaps;
	vector<eInterval> A,B;
	int i,j;
	normal NRV(0.6, 0.63);
	typedef map<string, vector<eInterval>>::iterator chrom_it;
	typedef vector<eInterval>::iterator A_it;
	int o_st, o_sp;
	double prob;
	for (chrom_it chrom = A_I.begin(); chrom!=A_I.end(); ++chrom ){
		A=chrom->second;
		if (B_I.find(chrom->first)!=B_I.end()){
			B=B_I[chrom->first];
			j = 0;
			for (A_it Ae = A.begin(); Ae!=A.end(); ++Ae){
				while (j+1 < B.size() and abs(B[j+1].stop - Ae->start) < abs(B[j].stop - Ae->start)){
					j++;
				}
				while (j < B.size() and B[j].start < Ae->stop){
					o_st 	= max(B[j].start, Ae->start);
					o_sp 	= min(B[j].stop, Ae->stop);
					prob	= min(1-NRV.cdf((o_sp - o_st) / 1000.0), NRV.cdf((o_sp - o_st) / 1000.0)) * 2.;
					if (prob > thresh){

						if(B[j].stop < Ae->start){
							overlaps.push_back(eInterval(B[j].stop, Ae->start, prob, chrom->first));
						}
						else{
							overlaps.push_back(eInterval(o_st, o_sp, prob, chrom->first));
						}
					}
					j++;
				}
				if (j -1 > 0){
					j--;
				}


			}



		}
	}
	return overlaps;


}

void writeOut(vector<eInterval> overlaps, string fileName){
	string header = ("track name='eRNA Predictions"
					"description=\"FStitch eRNA Predictions\" visibility=2 useScore=2 cgGrades=50 cgColour1=white cgColour2=yellow cgColour3=red height=30\n");
	ofstream writeFileHandle;
	writeFileHandle.open(fileName);
	writeFileHandle<<header;
	string line;
	typedef vector<eInterval>::iterator e_it;
	for (e_it o = overlaps.begin(); o!=overlaps.end(); ++o){
		line = (o->chrom + "\t" + to_string(o->start) + "\t" + to_string(o->stop) + "\teRNA=" + to_string(o->score) + "\t100\t.\t" +
				 to_string(o->start) + "\t" + to_string(o->stop) + "\t0,0,255\tN_A\n");
		writeFileHandle<<line;

	}
	writeFileHandle.close();


}


int eRNApredictions(string file1, string file2, double thresh, string out_dir){
	map<string, vector<eInterval>> A = makeIntervals(file1);
	map<string, vector<eInterval>> B = makeIntervals(file2);


	vector<eInterval> overlaps;
	if (A.empty() or B.empty()){
		cout<<"eRNA predictions not made"<<endl;
		cout<<"Exiting..."<<endl;
		return 0;
	}else{
		overlaps = getOverlaps(A,B, thresh);
		cout<<"There are "<<overlaps.size()<<" eRNA predictions"<<endl;
	}
	writeOut(overlaps, out_dir + "eRNA_Predictions.bed");

	return 1;
}




#endif /* ERNA_H_ */
