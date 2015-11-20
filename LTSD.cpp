/*
 * LTSD.cpp
 *
 *  Created on: 2015/11/18
 *      Author: Shunsuke Aihara (aihara@argmax.jp)
 */

#include "LTSD.h"

LTSD::LTSD(int windowsize, int samplingrate, int order, float e0, float e1, float lambda0, float lambda1){
	m_windowsize = windowsize;
	m_samplingrate = samplingrate;
	m_order = order;
	m_e0 = e0;
	m_e1 = e1;
	m_lambda0 = lambda0;
	m_lambda1 = lambda1;
	m_in = new float[m_windowsize];
	m_ltse = new float[m_windowsize];
	m_noiseProfile = new float[m_windowsize];
	for(int i=0; i< m_windowsize; i++){
		m_noiseProfile[i] = 0.0;
	}
	m_out = new float[m_windowsize];
	m_estimated = false;
	createWindow();
	m_fftreal = new ffft::FFTReal<float>(m_windowsize);
}

LTSD::~LTSD() {
	for (std::deque<unsigned char*>::iterator its = m_signalHistory.begin(); its != m_signalHistory.end(); its++){
		delete[] (*its);
	}
	for (std::deque<float*>::iterator ita = m_ampHistory.begin(); ita != m_ampHistory.end(); ita++){
		delete[] (*ita);
	}
	if (m_window != NULL){
		delete[] m_window;
	}
	delete[] m_in;
	delete[] m_out;
	delete[] m_ltse;
	delete[] m_noiseProfile;
	delete m_fftreal;

	if (ms != NULL){
		delete ms;
	}
}

bool LTSD::process(unsigned char* signal){

	for(int i=0; i<m_windowsize; i++){
			m_in[i]=((float(signal[i] - 127)) / 128.0) * m_window[i];
	}
	m_fftreal->do_fft(m_in, m_out);
	float *amp = new float[m_windowsize];
	for(int i=0; i<m_windowsize; i++){
		amp[i] = fabsf(m_out[i]);
	}
	unsigned char* sig = new unsigned char[m_windowsize];
	memcpy(sig, signal, sizeof(unsigned char) * m_windowsize);
	m_signalHistory.push_back(sig);
	m_ampHistory.push_back(amp);

	if (m_signalHistory.size() > m_order){
		if(!m_estimated){
			// 先頭orderフレームをノイズプロファイルに使う。本当はもっと長く取りたいけどとりあえずやってみる
			// MS法で逐次ノイズプロファイル更新を実装するほうが良さそう
			createNoiseProfile();
			m_estimated = true;
			ms = new MinimumStatistics(m_windowsize, m_samplingrate, m_noiseProfile);
		}

		if (ms != NULL){
			ms->process(amp);
			ms->updateNoiseProfile(m_noiseProfile);
		}
		//履歴長が指定以上なので、先頭を削除してltsd判定
		delete[] m_signalHistory[0];
		m_signalHistory.pop_front();
		delete[] m_ampHistory[0];
		m_ampHistory.pop_front();

		return isSignal();
	}else{
		return false;
	}
}

bool LTSD::isSignal(){
	ltse();
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
	float* amp = m_ampHistory.at(m_ampHistory.size() - 1);
	float sum = 0.0;
	for(int i = 0; i < m_windowsize; i++){
		sum += amp[i] ^ 2;
	}
	return 10 * log10(sum / m_windowsize);
}

unsigned char* LTSD::getSignal(){
	if (m_signalHistory.size() != m_order){
		return NULL;
	}else{
		unsigned char* src = m_signalHistory.at(m_signalHistory.size() - 1);
		unsigned char* dest = new unsigned char[m_windowsize];
		memcpy(dest, src, sizeof(unsigned char) * m_windowsize);
		return dest;
	}
}

void LTSD::ltse(){
	int i = 0;
	float amp;
	for(i=0;i < m_windowsize; i++){
		m_ltse[i] = 0.0;
	}
	for (std::deque<double*>::iterator ita = m_ampHistory.begin(); ita != m_ampHistory.end(); ita++){
		for(i=0;i < m_windowsize; i++){
			amp = (*ita)[i];
			if (m_ltse[i] < amp){
				m_ltse[i] = amp;
			}
		}
	}
	for(i=0;i < m_windowsize; i++){
		m_ltse[i] = m_ltse[i] * m_ltse[i];
	}
}

float LTSD::ltsd(){
	float sum = 0.0;
	for(int i = 0; i < m_windowsize; i++){
		sum += m_ltse[i] / m_noiseProfile[i];
	}
	return 10 * log10(sum / m_windowsize);
}

void LTSD::createNoiseProfile(){
	int i = 0;
	int s = m_ampHistory.size();
	for (std::deque<float*>::iterator ita = m_ampHistory.begin(); ita != m_ampHistory.end(); ita++){
		for(i=0;i < m_windowsize; i++){
			m_noiseProfile[i] += (*ita)[i];
		}
	}
	for(i=0;i < m_windowsize; i++){
		m_noiseProfile[i] = powf(m_noiseProfile[i] / s, 2);
	}
}

void LTSD::createWindow(){
	m_window = new double[m_windowsize];
	if (m_windowsize == 1){
		m_window[0] = 1.0;
	}else{
		double n = m_windowsize -1;
		double coef = M_PI * 2 / float(n);
		for (int i = 0; i < n; i++){
			m_window[i] = 0.54 - 0.46 * cos(coef * float(n));
		}
	}
}
