/*
 * window.h
 *
 *  Created on: Mar 10, 2014
 *      Author: joeyazo
 */

#ifndef WINDOW_H_
#define WINDOW_H_

double calc_density(vector<double> window){
	double D;
	for (int i = 0; i <window.size();i++){
		D+=window[i];
	}
	double length 	= window[window.size()-1] -window[0];

	return D / length;
}




#endif /* WINDOW_H_ */
