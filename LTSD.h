/*
 * LTSD.h

 *
 *  Created on: 2015/11/18
 *      Author: Shunsuke Aihara (aihara@argmax.jp)
 */

#ifndef LTSD_H_
#define LTSD_H_

#include <vector>
#include <deque>
#include "string.h"
#include "math.h"
#include "ffft/FFTReal.h"
#include "MinimumStatistics.h"
class LTSD {
public:
	// windowsizeは2の冪乗(256以上)で必ず偶数サイズにすること、orderは奇数(5, 7, 11)が望ましい
	LTSD(int winsize, int samprate, int order = 5,float e0 = 200.0, float e1 = 300.0, float lambda0 = 40.0, float lambda1 = 50.0);
	virtual ~LTSD();
	bool process(short *signal);
	short* getSignal(); // 取得したsignalは必ず利用後deleteすること
private:
	void createWindow();
	void calcLTSE();
	float calcLTSD();
	bool isSignal();
	void createNoiseProfile();
	float calcPower();
	int windowsize;
	int fftsize;
	int samplingrate;
	int m_order;
	float m_e0;
	float m_e1;
	float m_lambda0;
	float m_lambda1;
	float* window;
	float* ltse;
	float* noise_profile;
	bool estimated;

	ffft::FFTReal<float> *fftreal;
	MinimumStatistics *ms;
	float *fft_in;
	float* fft_out;
	std::deque<float*> amp_history;
	std::deque<short*> signal_history;
};

#endif /* LTSD_H_ */
