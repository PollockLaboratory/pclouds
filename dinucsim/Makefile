EXECUTABLE = ../../bin/DinucSimulator
SOURCE = $(wildcard *.c) $(wildcard *.cpp)
OBJS = $(patsubst %.c, %.o, $(patsubst %.cpp, %.o, $(SOURCE)))
DEPS = $(patsubst %.o, %.d, $(OBJS))
MISSING_DEPS = $(filter-out $(wildcard $(DEPS)),$(DEPS))
MISSING_DEPS_SOURCES = $(wildcard $(patsubst %.d, %.c, $(MISSING_DEPS))\
$(patsubst %.d, %.cpp, $(MISSING_DEPS)))

CC = /usr/bin/g++  
RM-F = rm -f
CPPFLAGS += -O3 -MMD -g -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -Wno-deprecated -I /usr/include -fpermissive 

.PHONY : all deps objs clean veryclean rebuild

all : $(EXECUTABLE)

deps : $(DEPS)

objs : $(OBJS)
	
clean: 
	rm ../../bin/DinucSimulator
	@$(RM-F) *.o
	@$(RM-F) *.d
	
ifneq ($(MISSING_DEPS),)
$(MISSING_DEPS) :
	@$(RM-F) $(patsubst %.d, %.o, $@)
endif

-include $(DEPS)

$(EXECUTABLE) : $(OBJS)
	$(CC) -o $(EXECUTABLE) $(OBJS)

