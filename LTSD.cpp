/*
 * LTSD.cpp
 *
 *  Created on: 2015/11/18
 *      Author: Shunsuke Aihara (aihara@argmax.jp)
 */

#include "LTSD.h"

LTSD::LTSD(int winsize, int samprate, int order, float e0, float e1, float lambda0, float lambda1){
	windowsize = winsize;
	samplingrate = samprate;
	m_order = order;
	m_e0 = e0;
	m_e1 = e1;
	m_lambda0 = lambda0;
	m_lambda1 = lambda1;
	fft_in = new float[windowsize];
	ltse = new float[windowsize];
	noise_profile = new float[windowsize];
	for(int i=0; i< windowsize; i++){
		noise_profile[i] = 0.0;
	}
	fft_out = new float[windowsize];
	estimated = false;
	createWindow();
	fftreal = new ffft::FFTReal<float>(windowsize);
}

LTSD::~LTSD() {
	for (std::deque<unsigned char*>::iterator its = signal_history.begin(); its != signal_history.end(); its++){
		delete[] (*its);
	}
	for (std::deque<float*>::iterator ita = amp_history.begin(); ita != amp_history.end(); ita++){
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

	if (ms != NULL){
		delete ms;
	}
}

bool LTSD::process(unsigned char* signal){

	for(int i=0; i<windowsize; i++){
			fft_in[i]=((float(signal[i] - 127)) / 128.0) * window[i];
	}
	fftreal->do_fft(fft_in, fft_out);
	float *amp = new float[windowsize];
	for(int i=0; i<windowsize; i++){
		amp[i] = fabsf(fft_out[i]);
	}
	unsigned char* sig = new unsigned char[windowsize];
	memcpy(sig, signal, sizeof(unsigned char) * windowsize);
	signal_history.push_back(sig);
	amp_history.push_back(amp);

	if (signal_history.size() > m_order){
		if(!estimated){
			// 先頭orderフレームをノイズプロファイルに使う。本当はもっと長く取りたいけどとりあえずやってみる
			// MS法で逐次ノイズプロファイル更新を実装するほうが良さそう
			createNoiseProfile();
			estimated = true;
			ms = new MinimumStatistics(windowsize, samplingrate, noise_profile);
		}

		if (ms != NULL){
			ms->process(amp);
			ms->updateNoiseProfile(noise_profile);
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
	float ltsd = ltsd();
	float e = calcPower();
	if (e < m_e0){
		if(ltsd > m_lambda0){
			return true;
		}else{
			return false;
		}
	}else if (e > m_e1){
		if(ltsd > m_lambda1){
			return true;
		}else{
			return false;
		}
	}else{
		float lamb = (m_lambda0-m_lambda1)/(m_e0/m_e1)*e + m_lambda0 - (m_lambda0-m_lambda1)/(1.0-(m_e0/m_e1));
		if (ltsd > lamb){
			return true;
		}else{
			return false;
		}
	}
}

float LTSD::calcPower(){
	float* amp = amp_history.at(amp_history.size() - 1);
	float sum = 0.0;
	for(int i = 0; i < windowsize; i++){
		sum += amp[i] ^ 2;
	}
	return 10 * log10(sum / windowsize);
}

unsigned char* LTSD::getSignal(){
	if (signal_history.size() != m_order){
		return NULL;
	}else{
		unsigned char* src = signal_history.at(signal_history.size() - 1);
		unsigned char* dest = new unsigned char[windowsize];
		memcpy(dest, src, sizeof(unsigned char) * windowsize);
		return dest;
	}
}

void LTSD::calcLTSE(){
	int i = 0;
	float amp;
	for(i=0;i < windowsize; i++){
		ltse[i] = 0.0;
	}
	for (std::deque<double*>::iterator ita = amp_history.begin(); ita != amp_history.end(); ita++){
		for(i=0;i < windowsize; i++){
			amp = (*ita)[i];
			if (ltse[i] < amp){
				ltse[i] = amp;
			}
		}
	}
	for(i=0;i < windowsize; i++){
		ltse[i] = ltse[i] * ltse[i];
	}
}

float LTSD::calcLTSD(){
	float sum = 0.0;
	for(int i = 0; i < windowsize; i++){
		sum += ltse[i] / noise_profile[i];
	}
	return 10 * log10(sum / windowsize);
}

void LTSD::createNoiseProfile(){
	int i = 0;
	int s = amp_history.size();
	for (std::deque<float*>::iterator ita = amp_history.begin(); ita != amp_history.end(); ita++){
		for(i=0;i < windowsize; i++){
			noise_profile[i] += (*ita)[i];
		}
	}
	for(i=0;i < windowsize; i++){
		noise_profile[i] = powf(noise_profile[i] / s, 2);
	}
}

void LTSD::createWindow(){
	window = new double[windowsize];
	if (windowsize == 1){
		window[0] = 1.0;
	}else{
		double n = windowsize -1;
		double coef = M_PI * 2 / float(n);
		for (int i = 0; i < n; i++){
			window[i] = 0.54 - 0.46 * cos(coef * float(n));
		}
	}
}
