/*
 * LTSD.h

 *
 *  Created on: 2015/11/18
 *      Author: Shunsuke Aihara (aihara@argmax.jp)
 */

#ifndef LTSD_H_
#define LTSD_H_

#include <string.h>
#include <vector>
#include <deque>
#include "math.h"
#include "ffft/FFTReal.h"

class LTSD {
public:
	// windowsizeは2の冪乗(256以上)、orderは奇数(5, 7, 11)が望ましい
	LTSD(int windowsize, int samplingrate, int order = 5,float e0 = 200.0, float e1 = 300.0, float lambda0 = 40.0, float lambda1 = 50.0);
	virtual ~LTSD();
	bool process(unsigned char *signal);
	unsigned char* getSignal(); // 取得したsignalは必ず利用後deleteすること
private:
	void createWindow();
	void ltse();
	float ltsd();
	bool isSignal();
	void createNoiseProfile();
	float calcPower();
	int m_windowsize;
	int m_samplingrate;
	int m_order;
	float m_e0;
	float m_e1;
	float m_lambda0;
	float m_lambda1;
	float* m_window;
	float* m_ltse;
	float* m_noiseProfile;
	bool m_estimated;

	ffft::FFTReal<float> *m_fftreal;
	MinimumStatistics *ms;
	float *m_in;
	float* m_out;
	std::deque<float*> m_ampHistory;
	std::deque<unsigned char*> m_signalHistory;
};

#endif /* LTSD_H_ */
