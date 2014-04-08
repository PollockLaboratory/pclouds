CC = g++
RM-F = rm -f
CPPFLAGS += -O3 -MMD -g -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -Wno-deprecated
.PHONY : all deps objs clean veryclean rebuild

all : preprocessor postprocessor

preprocessor: preprocessor.o
	$(CC) -o ../../bin/preprocessor preprocessor.o -I.

postprocessor: postprocessor.o
	$(CC) -o ../../bin/postprocessor postprocessor.o -I.

clean:
	rm *.o *.d ../../bin/preprocessor ../../bin/postprocessor
