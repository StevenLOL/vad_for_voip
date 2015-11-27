/*
 * MmseBasedNpe.h
 *
 *  Created on: 2015/11/27
 *      Author: aihara
 */

#ifndef MMSEBASEDNPE_H_
#define MMSEBASEDNPE_H_

#include "utils.h"
#include "string.h"
#include "math.h"

class MmseBasedNpe {
public:
	MmseBasedNpe(int size, double *noiseProfile);
	virtual ~MmseBasedNpe();
	void process(double *amp);
	void updateNoiseProfile(double *noise);
private:
	int fftsize;
	double* PH1mean;
	double alphaPH1mean;
	double alphaPSD;
	double q;
	double priorFact;
	double xiOptDb;
	double xiOpt;
	double logGLRFact;
	double GLRexp;

	double* noisePow;
	double* noisyPer;
	double* snrPost1;
	double* estimate;
	double* GLR;
	double* PH1;

};

#endif /* MMSEBASEDNPE_H_ */
