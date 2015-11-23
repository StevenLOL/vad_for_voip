CC :=g++
CFLAGS :=  -I.
LDFLAGS :=

ltsd.o: ltsd.cpp
	$(CC) $(CFLAGS) -o ltsd.o -c LTSD.cpp

ms.o: MinimumStatistics.cpp
	$(CC) $(CFLAGS) -o ms.o -c MinimumStatistics.cpp

all: ltsd.o ms.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o vad  ms.o ltsd.o vad.cpp
