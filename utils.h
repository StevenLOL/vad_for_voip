/*
 * utils.h
 *
 *  Created on: 2015/11/20
 *      Author: aihara
 */

#ifndef UTILS_H_
#define UTILS_H_

template<typename X> X *makeVector(int winsize, X init){
	X *vec = new X[winsize];
	for (int i = 0; i< winsize; i++){
		vec[i] = init;
	}
	return vec;
}

template<typename X> X sumVector(int winsize, X *vec){
	X s = 0;
	for(int i = 0; i< winsize; i++){
		s += vec[i];
	}
	return s;
}

template<typename X> void powerVector(int winsize, X *vec){
	for(int i = 0; i< winsize; i++){
		vec[i] = vec[i] * vec[i];
	}
}

template<typename X> void resetVector(int winsize, X *vec, X init){
	for(int i = 0; i< winsize; i++){
		vec[i] = init;
	}
}

#endif /* UTILS_H_ */
