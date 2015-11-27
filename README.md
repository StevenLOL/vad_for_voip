# vad_for_voip
maybe not working



## sample code

```vad.cpp
#include <iostream>
#include "LTSD.h"
using namespace std;

int main() {
	LTSD ltsd(256, 8000, 7);
	short signal[256 * 10000];
	for (int i = 0; i < 256 * 10000; i++){
		signal[i] = (rand() * rand()) % (32767 * 2 + 1) - 32767;
	}
	char *sig = (char *)&signal[0];

	for(int i = 0; i< 10000; i++){
		bool is_signal = ltsd.process(sig);
		(*sig)+=512;
		cout << is_signal << endl;
		if(is_signal){
			char *hoge = ltsd.getSignal();
			delete[] hoge;
		}
	}
	return 0;
}
```