EXECUTABLE = pclouds.exe

# List of sources to include
include sources.mk

OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
DEPS = $(patsubst %.o, %.d, $(OBJECTS))
MISSING_DEPS = $(filter-out $(wildcard $(DEPS)), $(DEPS))

# Implicit rule for cpp files
#  �$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c� 
# http://www.gnu.org/software/make/manual/html_node/Using-Implicit.html#Using-Implicit

# C++ compiler
CXX = g++ 

# C PreProcessor flags
CPPFLAGS = -MMD #-fprofile-arcs -ftest-coverage
# Special variable for implicit rule generation
# -MMD generates (-M) the dependency files (*.d) without 
# the system header files (-MM) and does not stop compilation at the 
# preprocessor state (-MMD).
# -fprofile-arcs -ftest-coverage for coverage checking

# C++ compiler flags
CXXFLAGS = -g #-pg -fprofile-arcs -ftest-coverage
# Special variable for implicit rule generation
# -g for debugging symbols
# -pg -fprofile-arcs -ftest-coverage for coverage checking

# GCC linker flags
LDFLAGS = -g #-pg -fprofile-arcs -ftest-coverage
# Special variable for implicit rule generation
# -g for debugging symbols
# -pg -fprofile-arcs -ftest-coverage for coverage checking



# This protects against any files called 'all', 'clean', etc. 
.PHONY : all clean rebuild 

all : $(EXECUTABLE)

clean: 
	$(RM) *.o
	$(RM) *.d
	$(RM) $(EXECUTABLE)

rebuild: clean all

# This protects against having a missing dependency file and remakes the object 
ifneq ($(MISSING_DEPS),)
$(MISSING_DEPS) :
	$(RM) $(patsubst %.d, %.o, $@)
endif


$(EXECUTABLE) : $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $(EXECUTABLE) $(OBJECTS)
