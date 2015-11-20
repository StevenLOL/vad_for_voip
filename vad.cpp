//============================================================================
// Name        : vad.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "LTSD.h"
using namespace std;

int main() {
	LTSD ltsd(256, 8000, 7);
	unsigned char signal[256 * 10000];
	for (int i = 0; i < 256 * 10000; i++){
		signal[i] = rand() % 256;
	}
	unsigned char *sig = &signal[0];

	for(int i = 0; i< 10000; i++){
		bool is_signal = ltsd.process(sig);
		(*sig)+=256;
		cout << is_signal << endl;
		if(is_signal){
			unsigned char *hoge = ltsd.getSignal();
			delete[] hoge;
		}
	}
	return 0;
}
