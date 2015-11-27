/*
 * MmseBasedNpe.cpp
 *
 *  Created on: 2015/11/27
 *      Author: aihara
 */

#include <MmseBasedNpe.h>

MmseBasedNpe::MmseBasedNpe(int size, double *noiseProfile) {
	// TODO Auto-generated constructor stub
	fftsize = size;
	alphaPH1mean = 0.9;
	alphaPSD = 0.8;

	q = 0.5; // a priori probability of speech presence:
	priorFact = q / (1 - q);
	xiOptDb = 15; // optimal fixed a priori SNR for SPP estimation
	xiOpt = pow(10, (xiOptDb / 10));
	logGLRFact = log(1.0 / (1 + xiOpt));
	GLRexp = xiOpt / (1 + xiOpt);

	PH1mean = makeVector(fftsize, 0.5);
	noisePow = new double[fftsize];
	memcpy(noisePow, noiseProfile, sizeof(double) * fftsize);
	noisyPer = new double[fftsize];
	snrPost1 = new double[fftsize];
	estimate = new double[fftsize];
	GLR = new double[fftsize];
	PH1 = new double[fftsize];
}

MmseBasedNpe::~MmseBasedNpe() {
	// TODO Auto-generated destructor stub
	delete[] noisePow;
	delete[] noisyPer;
	delete[] snrPost1;
	delete[] estimate;
	delete[] GLR;
}


void MmseBasedNpe::process(double* amp) {
	int i = 0;
	double tmp;
	for(i = 0; i< fftsize; i++){
		noisyPer[i] = amp[i] * amp[i];
		snrPost1[i] = noisyPer[i] / noisePow[i];

		tmp = logGLRFact + GLRexp *snrPost1[i];
		if(tmp > 200.0){
			tmp = 200.0;
		}
		GLR[i] = priorFact * exp(tmp);
		PH1 = GLR[i] / (1.0 + GLR[i]);
		PH1mean[i] = alphaPH1mean * PH1mean[i] + (1.0-alphaPH1mean) * PH1[i];
		if(PH1mean[i] > 0.99){
			if (PH1[i] > 0.99){
				PH1[i] = 0.99;
			}
		}
		estimate[i] = PH1[i] * noisePow[i] + (1.0 - alphaPSD) * noisyPer[i];
		noisePow[i] = alphaPSD *noisePow[i] + (1.0 - alphaPSD)*estimate[i];
	}
}

void MmseBasedNpe::updateNoiseProfile(double *noise){
	memcpy(noise, noisePow, sizeof(double) * fftsize);
}

