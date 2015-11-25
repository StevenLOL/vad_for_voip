/*
 * MinimumStatistics.cpp
 *
 *  Created on: 2015/11/19
 *      Author: Shunsuke Aihara(aihara@argmax.jp)
 */

#include "MinimumStatistics.h"

MinimumStatistics::MinimumStatistics(int size, int samplingrate, double *noise_profile) {
	// TODO Auto-generated constructor stub
	int i = 0;
	counter = 0;
	fftsize = size;
	double frametime = (double)fftsize / (double)samplingrate;
	snrexp = (-1.0)*frametime/0.064;
	av = 2.12;
	alpha_c_lambda = 0.7;
	U = 8;
	V = 12;
	D = U * V;
	M = 0.87;
	M2 = 0.64;
	H = 3.5;
	subwc = V - 1;
	ibuf  = 0;
	lmin_flag_lambda = makeVector(fftsize, int(0));
	alpha_max = 0.96;
	beta_max  = 0.8;
	float qeqmax = 14.0;
	qeqimin = 1/qeqmax;
	clear_max = 4294836225; // max of uint32_t
	actmin_lambda = makeVector(fftsize, clear_max);
	actmin_lambda_sub = makeVector(fftsize, clear_max);
	Pmin_u_lambda = makeVector(fftsize, clear_max);
	actbuf = new double*[U];
	for (i = 0;i < U; i++){
		actbuf[i] = makeVector(fftsize, clear_max);
	}
	P_lambda = new double[fftsize];
	memcpy(P_lambda, noise_profile, sizeof(double) * fftsize);
	sn2_lambda = new double[fftsize];
	memcpy(sn2_lambda, noise_profile, sizeof(double) * fftsize);
	eP_lambda = new double[fftsize];
	memcpy(eP_lambda, noise_profile, sizeof(double) * fftsize);
	eP2_lambda = new double[fftsize];
	memcpy(eP2_lambda, noise_profile, sizeof(double) * fftsize);
	powerVector(fftsize, eP2_lambda);
	Pmin_u_lambda = new double[fftsize];
	memcpy(Pmin_u_lambda, noise_profile, sizeof(double) * fftsize);


	power = new double[fftsize];

	alpha_lambda_hat = new double[fftsize];
	Qeq_lambda_inverse = new double[fftsize];
	Bmin_lambda = new double[fftsize];
	Bmin_lambda_sub = new double[fftsize];
	k_mod = new int[fftsize];
}

MinimumStatistics::~MinimumStatistics() {
	// TODO Auto-generated destructor stub
	int i = 0;
	for (i = 0;i < U; i++){
		delete[] actbuf[i];
	}
	delete[] actbuf;
	delete[] lmin_flag_lambda;
	delete[] actmin_lambda;
	delete[] actmin_lambda_sub;
	delete[] Pmin_u_lambda;
	delete[] P_lambda;
	delete[] sn2_lambda;
	delete[] eP_lambda;
	delete[] eP2_lambda;

	delete[] power;

	delete[] alpha_lambda_hat;
	delete[] Qeq_lambda_inverse;
	delete[] Bmin_lambda;
	delete[] Bmin_lambda_sub;
	delete[] k_mod;
}

void MinimumStatistics::process(double *amp){
	int i = 0;
	int j = 0;
	for (i = 0; i< fftsize; i++){
		power[i] = amp[i] * amp[i];
	}
	// eq9
	double tmp = (sumVector(fftsize, P_lambda)/sumVector(fftsize, power) - 1);
	double alpha_c_lambda_tilde = 1.0 / (tmp * tmp + 1.0);

	// eq10
	if (alpha_c_lambda_tilde > 0.7){
		tmp = alpha_c_lambda_tilde;
	}else{
		tmp = 0.7;
	}
	alpha_c_lambda = alpha_c_lambda * 0.7 + 0.3 * tmp;

	// eq11
	for(i = 0; i < fftsize; i++){
		tmp = (P_lambda[i] / sn2_lambda[i] - 1.0);
		alpha_lambda_hat[i] = (alpha_max / alpha_c_lambda) / (tmp * tmp + 1);
	}

	// eq12
	double snr = sumVector(fftsize, P_lambda) / sumVector(fftsize, sn2_lambda);
	tmp = powf(snr, snrexp);
	if (tmp > 0.3){
		tmp = 0.3;
	}
	for(i = 0; i < fftsize; i++){
		if (alpha_lambda_hat[i] < tmp){
			alpha_lambda_hat[i] = tmp;
		}
	}

	// eq4 smoothed periodgram
	for(i = 0; i < fftsize; i++){
		P_lambda[i] = alpha_lambda_hat[i] * P_lambda[i] + (1.0 - alpha_lambda_hat[i]) * power[i];
	}

	for(i = 0; i < fftsize; i++){
		// eq20
		tmp = alpha_lambda_hat[i] * alpha_lambda_hat[i];
		if (tmp > beta_max){
			tmp = beta_max;
		}
		eP_lambda[i] = tmp * eP_lambda[i] + (1.0 - tmp) * P_lambda[i];
		eP2_lambda[i] = tmp * eP2_lambda[i] + (1.0 - tmp) * P_lambda[i] * P_lambda[i];

		// eq22
		tmp = eP2_lambda[i] - eP_lambda[i] * eP_lambda[i];

		// eq23
		tmp = tmp / (sn2_lambda[i] * sn2_lambda[i] * 2.0);
		if (tmp > 0.5){
			tmp = 0.5;
		}
		if (tmp < qeqimin/(counter + 1)){
			tmp = qeqimin/(counter + 1);
		}
		Qeq_lambda_inverse[i] = tmp;
	}
	double eQ_lambda = sumVector(fftsize, Qeq_lambda_inverse) / fftsize;
	double Bc_lambda = 1.0 + av * sqrtf(eQ_lambda);

	// eq 16
	for(i = 0; i < fftsize; i++){
		// for overall window of length D
		tmp = (1.0 / Qeq_lambda_inverse[i] - 2 * M) / (1.0 - M);
		Bmin_lambda[i] = 1.0 + (D - 1) * 2.0 / tmp;
		// for subwindow U of length V
		tmp = (1.0 / Qeq_lambda_inverse[i] - 2 * M2) / (1.0 - M2);
		Bmin_lambda_sub[i] = 1.0 + (V - 1) * 2.0 / tmp;
	}

	// calc actmin,
	resetVector(fftsize, k_mod, 0); // reset to 0
	for(i = 0; i < fftsize; i++){
		// if (P * Bmin * Bc < actmin)
		tmp = P_lambda[i] * Bmin_lambda[i] * Bc_lambda;
		if (actmin_lambda[i] > tmp){
			actmin_lambda[i] = tmp;
			actmin_lambda_sub[i] = P_lambda[i] * Bmin_lambda_sub[i] * Bc_lambda;
			k_mod[i] = 1;
		}
	}

	if(0 < subwc && subwc < V-1){
		for(i = 0; i < fftsize; i++){
			// sample is NOT the fisrt or the last; middle of buffer allows a local minimum
			if (lmin_flag_lambda[i] + k_mod[i] >= 1){
				lmin_flag_lambda[i] = 1;
			}else{
				lmin_flag_lambda[i] = 0;
			}
			if ( Pmin_u_lambda[i] > actmin_lambda_sub[i]){
				Pmin_u_lambda[i] = actmin_lambda_sub[i];
			}
		}
		memcpy(sn2_lambda, Pmin_u_lambda, sizeof(double) * fftsize);
		subwc++;
	}else if(subwc >= V - 1){
		// store actmin_lamnda, note actbuf is NOT cyclic!
		ibuf = ibuf % U;
		memcpy(actbuf[ibuf], actmin_lambda, sizeof(double) * fftsize);
		ibuf++;

		// calc noise_slope_max
		double noise_slope_max;
		if(eQ_lambda < 0.03){
			noise_slope_max = 8.0;
		}else if(eQ_lambda < 0.05){
			noise_slope_max = 4.0;
		}else if(eQ_lambda < 0.06){
			noise_slope_max = 2.0;
		}else{
			noise_slope_max = 1.2;
		}

		// sample IS the last; end of buffer lets finishing subwindow process and a buffer switch
		for(i = 0; i < fftsize; i++){
			if (lmin_flag_lambda[i] - k_mod[i] >= 0){
				lmin_flag_lambda[i] = 1;
			}else{
				lmin_flag_lambda[i] = 0;
			}
			// find Pmin_u, the minimum of the last U stored value of actmin
			Pmin_u_lambda[i] = clear_max;
			for(j = 0; j < U; j++){
				if (Pmin_u_lambda[i] > actbuf[j][i]){
					Pmin_u_lambda[i] = actbuf[j][i];
				}
			}
			// replace all previously stored values of actmin by actminsub
			if (lmin_flag_lambda[i] && actmin_lambda_sub[i] < noise_slope_max * Pmin_u_lambda[i] && Pmin_u_lambda[i] < actmin_lambda_sub[i]){
				Pmin_u_lambda[i] = actmin_lambda_sub[i];
				for(j = 0; j < U; j++){
					actbuf[j][i] = Pmin_u_lambda[i];
				}
			}

		}
		resetVector(fftsize, lmin_flag_lambda, 0);
		resetVector(fftsize, actmin_lambda, clear_max);
		subwc = 0;
	}else{
		subwc++;
	}
	counter++;
}

void MinimumStatistics::updateNoiseProfile(double *noise){
	memcpy(noise, sn2_lambda, sizeof(double) * fftsize);
}
