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
	MinimumStatistics(int winsize, int samplingrate, float *noiseProfile);
	virtual ~MinimumStatistics();
	void process(float *amp);
	void updateNoiseProfile(float *noise);
private:
	int counter;
	int windowsize;
	float snrexp;
	float av;
	float alpha_c_lambda;
	int U;
	int V;
	int D;
	float M;
	float M2;
	float H;
	int subwc;
	int ibuf;
	int *lmin_flag_lambda;
	float alpha_max;
	float beta_max;
	float qeqimin;
	float clear_max;
	float *actmin_lambda;
	float *actmin_lambda_sub;
	float *Pmin_u_lambda;
	float **actbuf;
	float *P_lambda;
	float *sn2_lambda;
	float *eP_lambda;
	float *eP2_lambda;
	float *power;

	float *alpha_lambda_hat;
	float *Qeq_lambda_inverse;
	float *Bmin_lambda;
	float *Bmin_lambda_sub;
	int *k_mod;
};

#endif /* MINIMUMSTATISTICS_H_ */
