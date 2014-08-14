$(CXX) = g++
VPATH = src include
CFLAGS = -O3 -Wall -std=c++11
CPPFLAGS =  -I include
source = damping_fn.cpp main.cpp mc_fitting.cpp \
	molecule_funcs.cpp potential_funcs.cpp
objects = $(source:.cpp=.o)
####### MAIN EXECUTABLE ########
Fitting_MC : $(objects)
	$(CXX) -o $@ $(objects)
#### MAKING MOLECULE.H COMPULSORY FOR ALL FILES
#$(objects) : molecule.h
#### COMPILING OBJECTS WITH SPECIFIC OPTIONS
%.o : %.cpp molecule.h
	$(CXX) -c $(CFLAGS) $(CPPFLAGS) $< -o $@	
.PHONY: clean
clean:
	rm -rf $(objects) Fitting_MC
