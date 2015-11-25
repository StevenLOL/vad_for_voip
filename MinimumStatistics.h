/*
 * MinimumStatistics.h
 *
 *  Created on: 2015/11/19
 *      Author: Shunsuke Aihara(aihara@argmax.jp)
 */

#ifndef MINIMUMSTATISTICS_H_
#define MINIMUMSTATISTICS_H_

#include "utils.h"
#include "string.h"
#include "math.h"

class MinimumStatistics {
public:
	MinimumStatistics(int size, int samplingrate, double *noiseProfile);
	virtual ~MinimumStatistics();
	void process(double *amp);
	void updateNoiseProfile(double *noise);
private:
	int counter;
	int fftsize;
	double snrexp;
	double av;
	double alpha_c_lambda;
	int U;
	int V;
	int D;
	double M;
	double M2;
	double H;
	int subwc;
	int ibuf;
	int *lmin_flag_lambda;
	double alpha_max;
	double beta_max;
	double qeqimin;
	double clear_max;
	double *actmin_lambda;
	double *actmin_lambda_sub;
	double *Pmin_u_lambda;
	double **actbuf;
	double *P_lambda;
	double *sn2_lambda;
	double *eP_lambda;
	double *eP2_lambda;
	double *power;

	double *alpha_lambda_hat;
	double *Qeq_lambda_inverse;
	double *Bmin_lambda;
	double *Bmin_lambda_sub;
	int *k_mod;
};

#endif /* MINIMUMSTATISTICS_H_ */
