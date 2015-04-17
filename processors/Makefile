CC = g++
RM-F = rm -f
CPPFLAGS += -O3 -MMD -g -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -Wno-deprecated
.PHONY : all deps objs clean veryclean rebuild

all : preprocessor postprocessor dinucSim

preprocessor: preprocessor.o
	$(CC) -o preprocessor preprocessor.o

postprocessor: postprocessor.o
	$(CC) -o postprocessor postprocessor.o
	
dinucSim: dinucSim_cc_0.7_1Jan2014.cpp
	$(CC) -std=c++0x -o dinucSim dinucSim_cc_0.7_1Jan2014.cpp

clean:
	rm *.o *.d preprocessor postprocessor dinucSim
