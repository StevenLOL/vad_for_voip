/*
 * LTSD.cpp
 *
 *  Created on: 2015/11/18
 *      Author: Shunsuke Aihara (aihara@argmax.jp)
 */

#include "LTSD.h"

//#include <android/log.h>
//#define LOGE(...) ((void)__android_log_print(ANDROID_LOG_VERBOSE, "vaddsp-jni", __VA_ARGS__))


LTSD::LTSD(int winsize, int samprate, int order, double e0, double e1, double lambda0, double lambda1){
	windowsize = winsize;
	fftsize = winsize / 2;
	samplingrate = samprate;
	m_order = order;
	m_e0 = e0;
	m_e1 = e1;
	m_lambda0 = lambda0;
	m_lambda1 = lambda1;
	fft_in = new double[windowsize];
	ltse = new double[fftsize];
	noise_profile = new double[fftsize];
	for(int i=0; i< fftsize; i++){
		noise_profile[i] = 0.0;
	}
	fft_out = new double[windowsize];
    for(int i=0; i< fftsize; i++){
      fft_out[i] = 0.0;
    }
	estimated = false;
	createWindow();
	fftreal = new ffft::FFTReal<double>(fftsize);
}

LTSD::~LTSD() {
	for (std::deque<short*>::iterator its = signal_history.begin(); its != signal_history.end(); its++){
		delete[] (*its);
	}
	for (std::deque<double*>::iterator ita = amp_history.begin(); ita != amp_history.end(); ita++){
		delete[] (*ita);
	}
	if (window != NULL){
		delete[] window;
	}
	delete[] fft_in;
	delete[] fft_out;
	delete[] ltse;
	delete[] noise_profile;
	delete fftreal;

	if (mmse != NULL){
		delete mmse;
	}
}

bool LTSD::process(char *input){
    short *signal = (short *)input;
	for(int i=0; i<windowsize; i++){
		fft_in[i]=(double(signal[i]) / 32767.0) * window[i];
		fft_out[i] = 0.0;
	}
	fftreal->do_fft(fft_out, fft_in);
	double *amp = new double[fftsize];
	for(int i=0; i<fftsize; i++) {
		if (!std::isinf(fft_out[i]) && !std::isnan(fft_out[i])) {
			amp[i] = fabs(fft_out[i]);
		}
	}

	short* sig = new short[windowsize];
	memcpy(sig, signal, sizeof(short) * windowsize);
	signal_history.push_back(sig);
	amp_history.push_back(amp);

	if (signal_history.size() > m_order){
		if(!estimated){
			createNoiseProfile();
			estimated = true;
			mmse = new MmseBasedNpe(fftsize, noise_profile);
		}
		if (mmse != NULL){
			mmse->process(amp);
			mmse->updateNoiseProfile(noise_profile);
		}
		//履歴長が指定以上なので、先頭を削除してltsd判定
		delete[] signal_history[0];
		signal_history.pop_front();
		delete[] amp_history[0];
		amp_history.pop_front();

		return isSignal();
	}else{
		return false;
	}
}

bool LTSD::isSignal(){
	calcLTSE();
	double ltsd = calcLTSD();
	double e = calcPower();
    double e2 = calcNoisePower();
	double sn = fabs(e - e2);
    double lamb = (m_lambda0 - m_lambda1) / (m_e0 / m_e1) * e2 + m_lambda0 -
                  (m_lambda0 - m_lambda1) / (1.0 - (m_e1 / m_e0));

    LOGE("signal: %f, noise: %f, lambda:%f", e, e2, lamb);


	if (e2 < m_e0){
		if(ltsd > m_lambda0){
			return true;
		}else{
			return false;
		}
	}else if (e2 > m_e1){
		if(ltsd > m_lambda1){
			return true;
		}else{
			return false;
		}
	}else {
        if (ltsd > lamb) {
            return true;
        } else {
            return false;
        }
    }
}

double LTSD::calcPower(){
	double* amp = amp_history.at(amp_history.size() - 1);
	double sum = 0.0;
	for(int i = 0; i < fftsize; i++){
		sum += amp[i] * amp[i];
	}
	return 10 * log10((sum / fftsize) / (1.0e-6 * 1.0e-6));
}

double LTSD::calcNoisePower(){
    double s = 0.0;
    for(int i = 0; i < fftsize; i++){
        s += noise_profile[i];
    }
    return 10 * log10((s / fftsize) / (1.0e-6 * 1.0e-6));
}

char* LTSD::getSignal(){
	if (signal_history.size() != m_order){
		return NULL;
	}else{
		short* src = signal_history.at(signal_history.size() - 1);
		short* dest = new short[windowsize];
		memcpy(dest, src, sizeof(short) * windowsize);
		return (char *)dest;
	}
}

void LTSD::calcLTSE(){
	int i = 0;
	double amp;
	for(i=0;i < fftsize; i++){
		ltse[i] = 0.0;
	}
	for (std::deque<double*>::iterator ita = amp_history.begin(); ita != amp_history.end(); ita++){
		for(i=0;i < fftsize; i++){
			amp = (*ita)[i];
			if (ltse[i] < amp){
				ltse[i] = amp;
			}
		}
	}
	for(i=0;i < fftsize; i++){
		ltse[i] = ltse[i] * ltse[i];
	}
}

double LTSD::calcLTSD(){
	double sum = 0.0;
	for(int i = 0; i < fftsize; i++){
		sum += ltse[i] / noise_profile[i];
	}
	return 10 * log10(sum / fftsize);
}

void LTSD::createNoiseProfile(){
	int i = 0;
	double s = (double)amp_history.size();
	for (std::deque<double*>::iterator ita = amp_history.begin(); ita != amp_history.end(); ita++){
        double *x = (*ita);
		for(i=0;i < fftsize; i++){
			noise_profile[i] += x[i];
		}
	}
    double sum = 0.0;
	for(i=0;i < fftsize; i++){
        noise_profile[i] = pow(noise_profile[i] / s, 2);
	}
}

void LTSD::createWindow(){
	window = new double[windowsize];
	if (windowsize == 1){
		window[0] = 1.0;
	}else{
		double n = windowsize -1;
		double coef = M_PI * 2 / double(n);
		for (int i = 0; i < n; i++){
			window[i] = 0.54 - 0.46 * cos(coef * double(i));
		}
	}
}

void LTSD::updateParams(double e0, double e1, double lambda0, double lambda1){
 	m_e0 = e0;
	m_e1 = e1;
	m_lambda0 = lambda0;
	m_lambda1 = lambda1; 
}

