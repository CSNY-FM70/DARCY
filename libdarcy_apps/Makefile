EXE=$(addsuffix .exe, $(basename $(wildcard *.cpp)))
CPPC=g++
CPPC_FLAGS=-g -Wall -Wextra -pedantic -Ofast -march=native -std=c++11
EIGEN_DIR=$(HOME)/Software/include/Eigen

BASE_DIR=..
LIBDARCY_INC_DIR=$(BASE_DIR)/libdarcy/include
LIBDARCY_APPS_INC_DIR=$(BASE_DIR)/libdarcy_apps/include

apps : $(EXE)

all :	$(EXE)
	$(doc)

%.exe : %.cpp
	$(CPPC) $(CPPC_FLAGS) -I$(EIGEN_DIR) -I$(LIBDARCY_INC_DIR) -I$(LIBDARCY_APPS_INC_DIR) $< -o $@

doc:
	cd doc && $(MAKE)

clean :
	cd doc && $(MAKE) clean
	rm -fr $(EXE)
	rm -fr *.csv

cleandoc:
	cd doc && $(MAKE) clean

.PHONY: all doc clean cleandoc apps
