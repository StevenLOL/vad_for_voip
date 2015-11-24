CC :=g++
CFLAGS :=  -I.
LDFLAGS :=


all: ltsd.o ms.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o vad  ms.o ltsd.o vad.cpp

ltsd.o: ltsd.cpp
	$(CC) $(CFLAGS) -o ltsd.o -c LTSD.cpp

ms.o: MinimumStatistics.cpp
	$(CC) $(CFLAGS) -o ms.o -c MinimumStatistics.cpp

clean:
	rm -rf *.o
	rm -rf vad
