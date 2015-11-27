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
#include <cmath>
#include "math.h"
#include "utils.h"
#include "ffft/FFTReal.h"
//#include "MinimumStatistics.h"
#include "MmseBasedNpe.h"

class LTSD {
public:
	// windowsizeは2の冪乗(256以上)で必ず偶数サイズにすること。
	// この値は、short型で換算した場合のサンプル数(つまりbyte配列の長さ/2)。
	// orderは奇数(5, 7, 11)が望ましい
	// sampling rateは、録音時のサンプリングレートを用いる
	LTSD(int winsize, int samprate, int order = 7,double e0 = -30.0, double e1 = -20.0, double lambda0 = 20.0, double lambda1 = 10.0);
	virtual ~LTSD();
	bool process(char *input);
	char* getSignal(); // 取得したsignalは必ず利用後deleteすること。byteで返しているが、アーキテクチャのエンディアンでのshort型の配列になっている
private:
	void createWindow();
	void calcLTSE();
	double calcLTSD();
	bool isSignal();
	void createNoiseProfile();
	double calcPower();
	double calcNoisePower();
	int windowsize;
	int fftsize;
	int samplingrate;
	int m_order;
	double m_e0;
	double m_e1;
	double m_lambda0;
	double m_lambda1;
	double* window;
	double* ltse;
	double* noise_profile;
	bool estimated;

	ffft::FFTReal<double> *fftreal;
	MmseBasedNpe *mmse;
	double *fft_in;
	double* fft_out;
	std::deque<double*> amp_history;
	std::deque<short*> signal_history;
};

#endif /* LTSD_H_ */
