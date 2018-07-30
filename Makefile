#modied from htslib makefile
FLAGS=-O3

CXXFLAGS += $(FLAGS)

CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)

all: ngsRelate

# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif

.PHONY: misc clean test

-include $(OBJ:.o=.d)

ifdef HTSSRC
%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp -I$(HTS_INCDIR) -D__WITH_BCF__
	$(CXX) -MM $(CXXFLAGS) $*.cpp -I$(HTS_INCDIR) -D__WITH_BCF__ >$*.d

ngsRelate: $(OBJ)
	$(CXX) $(FLAGS)  -o ngsRelate *.o $(HTS_LIBDIR) -lz -lm -lbz2 -llzma -lpthread -lcurl -D__WITH_BCF__
else
%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -D__WITH_BCF__ $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -D__WITH_BCF__ $*.cpp >$*.d

ngsRelate: $(OBJ)
	$(CXX) $(FLAGS)  -o ngsRelate *.o -lz -lpthread -lhts  -lbz2 -llzma 
endif

clean:
	rm  -f *.o *.d ngsRelate *~
