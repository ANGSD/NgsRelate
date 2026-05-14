#modied from htslib makefile
FLAGS=-O3

CXXFLAGS += $(FLAGS)

CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)
VERSION_HEADER=VERSION.h

BASE_VERSION := 2.1.0

ifeq ($(origin PACKAGE_VERSION), undefined)
PACKAGE_VERSION := $(BASE_VERSION)
ifneq ("$(wildcard .git)","")
GIT_DESCRIBE := $(shell git describe --tags --always --dirty 2>/dev/null)
ifneq ("$(strip $(GIT_DESCRIBE))","")
PACKAGE_VERSION := $(GIT_DESCRIBE)
endif
endif
endif

all: $(VERSION_HEADER) ngsRelate

# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif

.PHONY: misc clean test force

-include $(OBJ:.o=.d)

$(VERSION_HEADER): force
	@printf '#ifndef NGSRELATE_VERSION_H\n#define NGSRELATE_VERSION_H\n#define NGSRELATE_VERSION "%s"\n#endif\n' "$(PACKAGE_VERSION)" > $(VERSION_HEADER).tmp
	@cmp -s $(VERSION_HEADER).tmp $(VERSION_HEADER) || mv $(VERSION_HEADER).tmp $(VERSION_HEADER)
	@rm -f $(VERSION_HEADER).tmp

force:

ifdef HTSSRC
%.o: %.cpp $(VERSION_HEADER)
	$(CXX) -c  $(CXXFLAGS) $*.cpp -I$(HTS_INCDIR) -D__WITH_BCF__ ## -DDB_EMIS -DDB_MP
	$(CXX) -MM $(CXXFLAGS) $*.cpp -I$(HTS_INCDIR) -D__WITH_BCF__ >$*.d

ngsRelate: $(VERSION_HEADER) $(OBJ)
	$(CXX) $(FLAGS)  -o ngsRelate *.o $(HTS_LIBDIR) -lz -lm -lbz2 -llzma -lpthread -lcurl -D__WITH_BCF__
else
%.o: %.cpp $(VERSION_HEADER)
	$(CXX) -c  $(CXXFLAGS)  -D__WITH_BCF__ $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -D__WITH_BCF__ $*.cpp >$*.d

ngsRelate: $(VERSION_HEADER) $(OBJ)
	$(CXX) $(FLAGS)  -o ngsRelate *.o -lz -lpthread -lhts  -lbz2 -llzma 
endif

clean:
	rm  -f *.o *.d ngsRelate $(VERSION_HEADER) $(VERSION_HEADER).tmp *~

.PHONY: all clean install install-all install-misc misc test

test:
	echo "Only subset of analyses is being tested"
	cd test;./testAll.sh ../ngsRelate
