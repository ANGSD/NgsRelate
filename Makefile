# modified from htslib makefile
FLAGS ?= -O3

CXXFLAGS += $(FLAGS)
CPPFLAGS += -D__WITH_BCF__

TARGET := ngsRelate
CXXSRC := $(wildcard *.cpp)
OBJ := $(CXXSRC:.cpp=.o)
VERSION_HEADER := VERSION.h
LDLIBS_COMMON := -lz -lpthread

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

all: $(VERSION_HEADER) $(TARGET)

# Adjust $(HTSSRC) to point to your top-level htslib directory
USE_HTSSRC := 0
ifneq ($(strip $(HTSSRC)),)
ifneq ($(strip $(HTSSRC)),systemwide)
USE_HTSSRC := 1
$(info HTSSRC defined)
HTS_ROOT := $(realpath $(strip $(HTSSRC)))
ifneq ("$(wildcard $(HTS_ROOT)/include/htslib/hts.h)","")
HTS_INCDIR := $(HTS_ROOT)/include
else
HTS_INCDIR := $(HTS_ROOT)
endif
ifneq ("$(wildcard $(HTS_ROOT)/lib/libhts.a)","")
HTS_LIBDIR := $(HTS_ROOT)/lib/libhts.a
else
HTS_LIBDIR := $(HTS_ROOT)/libhts.a
endif
else
$(info HTSSRC set to systemwide, using -lhts)
endif
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif

.PHONY: all clean misc test force

-include $(OBJ:.o=.d)

$(VERSION_HEADER): force
	@printf '#ifndef NGSRELATE_VERSION_H\n#define NGSRELATE_VERSION_H\n#define NGSRELATE_VERSION "%s"\n#endif\n' "$(PACKAGE_VERSION)" > $(VERSION_HEADER).tmp
	@cmp -s $(VERSION_HEADER).tmp $(VERSION_HEADER) || mv $(VERSION_HEADER).tmp $(VERSION_HEADER)
	@rm -f $(VERSION_HEADER).tmp

force:

ifeq ($(USE_HTSSRC),1)
%.o: %.cpp $(VERSION_HEADER)
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -I$(HTS_INCDIR) $< -o $@ ## -DDB_EMIS -DDB_MP
	$(CXX) -MM $(CXXFLAGS) $(CPPFLAGS) -I$(HTS_INCDIR) $< >$*.d

$(TARGET): $(VERSION_HEADER) $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJ) $(HTS_LIBDIR) $(LDLIBS_COMMON) -lbz2 -llzma -lm -lcurl
else
%.o: %.cpp $(VERSION_HEADER)
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
	$(CXX) -MM $(CXXFLAGS) $(CPPFLAGS) $< >$*.d

$(TARGET): $(VERSION_HEADER) $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJ) $(LDLIBS_COMMON) -lhts
endif

clean:
	rm -f *.o *.d $(TARGET) $(VERSION_HEADER) $(VERSION_HEADER).tmp *~

test:
	echo "Only subset of analyses is being tested"
	cd test;./testAll.sh ../$(TARGET)
