SHELL := /bin/bash

DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
EBROOTHTSLIB ?= ${PWD}/src/htslib/
BLEND2DSRC ?= ${PWD}/src/blend2d/
BLEND2DLIB ?= ${BLEND2DSRC}/build/libblend2d.a

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir ?= $(exec_prefix)/bin

# Flags
CXX=g++
CXXFLAGS += -std=c++17 -isystem ${EBROOTHTSLIB} -isystem ${BLEND2DSRC} -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS -fno-strict-aliasing -fpermissive
LDFLAGS += -L${EBROOTHTSLIB} -lboost_iostreams -lboost_filesystem -lboost_program_options -lboost_date_time

# Flags for static compile
ifeq (${STATIC}, 1)
	LDSTATIC = --static
	HTSSTATIC = --disable-libcurl
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz -llzma -lbz2 -ldeflate -lrt
else
	LDSTATIC =
	HTSSTATIC =
	LDFLAGS += -pthread -lhts -lz -llzma -lbz2 -lrt -Wl,-rpath,${EBROOTHTSLIB}
endif

# Flags for debugging, profiling and releases
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG
endif

# Submodules
ifeq (${EBROOTHTSLIB}, ${PWD}/src/htslib/)
	SUBMODULES += .htslib
endif
ifeq (${BLEND2DSRC}, ${PWD}/src/blend2d/)
	SUBMODULES += .blend2d
endif

# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
BLEND2DSOURCES = $(wildcard src/blend2d/blend2d/*.h) $(wildcard src/blend2d/blend2d/*/*.h)
SOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
BUILT_PROGRAMS = src/wally
TARGETS = ${SUBMODULES} ${BUILT_PROGRAMS}

all:   	$(TARGETS)

.htslib: $(HTSLIBSOURCES)
	if [ -r src/htslib/Makefile ]; then cd src/htslib && autoreconf -i && ./configure ${HTSSTATIC} --disable-plugins && $(MAKE) && $(MAKE) lib-static && cd ../../ && touch .htslib; fi

.blend2d: ${BLEND2DSOURCES}
	if [ -f ${BLEND2DSRC}/CMakeLists.txt ]; then cd ${BLEND2DSRC} && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBLEND2D_STATIC=ON -DBLEND2D_NO_JIT=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON && cmake --build build --target blend2d -j$(shell nproc) && cd ../../ && touch .blend2d; fi

src/wally: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(BLEND2DLIB) $(LDFLAGS)

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	if [ -r src/htslib/Makefile ]; then cd src/htslib && $(MAKE) clean; fi
	rm -f $(TARGETS) $(TARGETS:=.o) ${SUBMODULES}
	rm -rf ${BLEND2DSRC}/build/

distclean: clean
	rm -f ${BUILT_PROGRAMS}

.PHONY: clean distclean install all
