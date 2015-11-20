/*
 * utils.h
 *
 *  Created on: 2015/11/20
 *      Author: aihara
 */

#ifndef UTILS_H_
#define UTILS_H_


float *makeZeroVector(int winsize){
	float *vec = new float[winsize];
		for (int i = 0; i< winsize; i++){
			vec[i] = 0.0;
		}
		return vec;
}

int *makeZeroIntVector(int winsize){
	int *vec = new int[winsize];
		for (int i = 0; i< winsize; i++){
			vec[i] = 0;
		}
		return vec;
}

float *makeVector(int winsize, float init){
	float *vec = new float[winsize];
		for (int i = 0; i< winsize; i++){
			vec[i] = init;
		}
		return vec;
}


float sumVector(int winsize, float *vec){
	float s = 0.0;
	for(int i = 0; i< winsize; i++){
		s += vec[i];
	}
	return s;
}

void powerVector(int winsize, float *vec){
	for(int i = 0; i< winsize; i++){
		vec[i] = vec[i] * vec[i];
	}
}

void resetIntVector(int winsize, int *vec, int init){
	for(int i = 0; i< winsize; i++){
		vec[i] = init;
	}
}

void resetVector(int winsize, float *vec, float init){
	for(int i = 0; i< winsize; i++){
		vec[i] = init;
	}
}

#endif /* UTILS_H_ */
