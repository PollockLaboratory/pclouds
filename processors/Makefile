CC = g++
RM-F = rm -f
CPPFLAGS += -O3 -MMD -g -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -Wno-deprecated
.PHONY : all deps objs clean veryclean rebuild

all : preprocessor postprocessor

preprocessor: preprocessor.o
	$(CC) -o preprocessor preprocessor.o

postprocessor: postprocessor.o
	$(CC) -o postprocessor postprocessor.o

clean:
	rm *.o *.d preprocessor postprocessor
