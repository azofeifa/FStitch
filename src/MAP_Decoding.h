/*
 * MAP_Decoding.h
 *
 *  Created on: Mar 4, 2014
 *      Author: joeyazo
 */

#ifndef MAP_DECODING_H_
#define MAP_DECODING_H_

using namespace std;


vector< vector<double> > runForward(vector<MEMM> data_array, vector<vector<double> > A, vector<logistic_equation> emission_eqs, vector<int> feature_dimension, bool SHOW){
	vector<double> xs;
	int K 		= A.size(); //number of states
	int N 		= data_array.size(); //number of time points
	int t 		= 0;
	double pi	= 1.0 / emission_eqs.size();
	double sum	= 0;
	double tsum	= 0;

	int f,i,j;

	vector< vector<double> > alpha(N);
	for (int t = 0; t < N; t++){
		alpha[t] = vector<double>(K);
	}

	while (t < N){
		xs 	= data_array[t].getFeatures(feature_dimension);
		if (t==0)
		{
			for (i = 0; i< K; i++){
				alpha[t][i] 	= pi*B(xs, emission_eqs, i);
			}
		}
		else{
			tsum=0;
			for (i = 0; i< K; i++){
				sum				= B(xs, emission_eqs, i);
				for (j = 0; j < K; j++){
					sum+=(alpha[t-1][j] * A[j][i]);
				}
				alpha[t][i] 	= sum;
				tsum+=sum;
			}
			//Normalize
			for (i = 0; i< K; i++){
				alpha[t][i]	= alpha[t][i] / tsum;
			}
		}
		t++;
	}
	return alpha;
}
vector< vector<double> > runBackward(vector<MEMM> data_array, vector<vector<double> > A, vector<logistic_equation> emission_eqs, vector<int> feature_dimension, bool SHOW){
	vector<double> xs;
	int K 		= A.size(); //number of states
	int N 		= data_array.size(); //number of time points
	int t 		= N-1;
	double pi	= 1.0 / emission_eqs.size();
	double sum	= 0;
	double tsum	= 0;

	int f,i,j;

	vector< vector<double> > beta(N);
	for (int t = 0; t < N; t++){
		beta[t] = vector<double>(K);
	}

	while (t > -1){
		if (t==(N-1))
		{
			xs 	= data_array[t].getFeatures(feature_dimension);

			for (i = 0; i< K; i++){
				beta[t][i] 	= pi*B(xs, emission_eqs, i);
			}
		}
		else{
			xs 	= data_array[t+1].getFeatures(feature_dimension);

			tsum	= 0;
			for (i=0;i<K;i++){
				sum = 0;
				for (j=0;j<K;j++){
					sum+=(B(xs, emission_eqs,j) * A[i][j] * beta[t+1][j]);
				}
				beta[t][i]=sum;
				tsum+=sum;
			}
			for (i=0;i<K;i++){
				beta[t][i]	= beta[t][i] / tsum;
			}
		}
		t--;
	}
	return beta;
}
vector< vector< double> > runForwardAndBackward(vector<MEMM> data_array, vector<vector<double> > A, vector<logistic_equation> emission_eqs, vector<int> feature_dimension, bool SHOW){
	clock_t time;
	time = clock();
	vector< vector<double> > alpha_trellis 	= runForward(data_array, A, emission_eqs,feature_dimension, SHOW);
	vector< vector<double> > beta_trellis 	= runBackward(data_array, A, emission_eqs, feature_dimension, SHOW);
	int N 									= data_array.size(); //number of time points
	int K 									= A.size();
	int t,i;
	vector< vector< double> > probs(N);
	for (t = 0; t < N; t++){
		probs[t]=vector<double>(K);
	}
	for (t = 0; t < N; t++){
		for (i = 0; i < K; i++){
			probs[t][i] = alpha_trellis[t][i] * beta_trellis[t][i];
		}
	}
	return probs;

}



#endif /* MAP_DECODING_H_ */
