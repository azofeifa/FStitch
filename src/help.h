/*
 * help.h
 *
 *  Created on: May 3, 2014
 *      Author: joeyazo
 */

#ifndef HELP_H_
#define HELP_H_
//string fileName 		= "";
//string refGeneFileName 	= "";
//string Training_File	= "";
//double slope_density 	= 0.01; //Logisitic Regression Parameter Density
//double slope_coverage 	= 2; //Logisitic Regression Parameter Coverage
//double weight 			= 1.0; //Weight Coverage Logistic  vs Density Logistic
//double A 				= 0.99; //	Self Transition Loop Parameter
//bool example 			= 0;
//bool strand_bool 		= 0;
//string strand 			= "";
//int dim 				= 3;
//bool viterbi_BOOL 		= 1;
//vector<double> W;
//string userFileName		= "";
//int BM					= 20;
//double alpha 			= 1.;


void runHelp(){
	cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	cout<<"\t\t\tLogitGRO Parameter Use"<<endl;
	cout<<"\tPlease see software manual for more detailed exposition"<<endl<<endl;


	cout<<"-i\t...\t--input-bedgraph-file"<<endl;
	cout<<("\tPath to BedGraphFile...must either be positive or negative strand"
				"\n\tAnd file name must specify strand by having either .pos.\n"
				"\tor .+. to specify positive strand"
				"\n\tor .neg. or .-. to specify negative strand")<<endl<<endl;
	cout<<"-st\t...\t--strand"<<endl;
	cout<<("\tIf you did not place a .pos. or .+. (etc.) in bedgraph file name\n"
				"\tyou can specify the strand either + or - with this flag")<<endl<<endl;

	cout<<"-rf\t...\t--refseq-file"<<endl;
	cout<<("\teither hg18 or hg19, default is hg19\n"
					"\tsee refseqannotations directory for these files")<<endl<<endl;


	cout<<"-t\t...\t--run-unit-test"<<endl;
	cout<<("\truns unit tests")<<endl<<endl;


	cout<<"-o\t...\t--output-directory"<<endl;
	cout<<("\tSpecify path to where you want the IGV segments outputted to.\n"
			"\tDefault is segmentation_outputs/ LogitGRO directory")<<endl<<endl;

	cout<<"-uf\t...\t--name-of-IGV-file"<<endl;
	cout<<("\tSpecficy the name that you want to call IGV output file\n"
			"\tDefault appends _IGV.bed to whatever the input bedgraph file was called")<<endl<<endl;
	cout<<"-l\t...\t--learn"<<endl;
	cout<<("\tPath to training file; see training_files/ for examples on use"
			"\n\tOr manuel of types of acceptable training data formats")<<endl<<endl;

	cout<<"-d\t...\t--feature-dimension"<<endl;
	cout<<("\tNumber of Features you wish to consider;\n"
			"\tif 2 then you only consider the bias term and contig length\n"
			"\tif 3 then you consider, bias, contig length and average coverage\n"
			"\tetc. please software manuel for full exposition")<<endl<<endl;


	cout<<"-al\t...\t--alpha-learning-rate"<<endl;
	cout<<"\tLearning rate for raphson newton method\n\tdefault is 0.5"<<endl<<endl;

	cout<<"-bw\t...\t--baum-welch-iterations"<<endl;
	cout<<"\tMax Number of iterations for Baum-Welch\n\tDefault is 15"<<endl<<endl;

	cout<<"-w\t...\t--logistic-regression-weights"<<endl;
	cout<<("\tIf -learn is not specified you can manually put in the weights\n"
			"\tof the logistic regression\n"
			"\tIn order it is, bias term weight, weight for "
			"\n\tcontig lenght, weight for coverage mean")<<endl<<endl;

	cout<<"-a\t...\t--self-transition-HMM-weight"<<endl;
	cout<<("\tIf -learn is not specified you can manually specficy the markov transition parameter\n"
					"\tOne number makes it a symetric markov process where\n"
					"\tthe self transition loop is specified here\n"
					"\tdefault is 0.9 for both Inactive and Active classes")<<endl<<endl;
	cout<<"-cd\t...\t--feature-code"<<endl;
	cout<<"\tThis is what features you want to consider..."<<endl;
	cout<<"\t1\tContig Length"<<endl;
	cout<<"\t2\tMean Coverage"<<endl;
	cout<<"\t3\tTotal Coverage"<<endl;
	cout<<"\t4\tMax Coverage"<<endl;
	cout<<"\t5\tVariance Coverage"<<endl;








	cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;




}





#endif /* HELP_H_ */
