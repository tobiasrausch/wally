DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
EBROOTHTSLIB ?= ${PWD}/src/htslib/
OPENCVSRC ?= ${PWD}/src/opencv/
OPENCV ?= ${PWD}/src/ocv/

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir ?= $(exec_prefix)/bin

# Flags
CXX=g++
CXXFLAGS += -isystem ${EBROOTHTSLIB} -isystem ${OPENCV}/include/opencv4 -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS -fno-strict-aliasing -fpermissive
LDFLAGS += -L${EBROOTHTSLIB} -L${EBROOTHTSLIB}/lib -L ${OPENCV}/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time -lopencv_core -lopencv_imgcodecs -lopencv_imgproc -lopencv_highgui

# Flags for static compile
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz -llzma -lbz2
else
	LDFLAGS += -pthread -lhts -lz -llzma -lbz2 -Wl,-rpath,${EBROOTHTSLIB}
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
ifeq (${EBROOTHTSLIB}, ${PWD}/src/htslib/)
	SUBMODULES += .htslib
endif
ifeq (${OPENCVSRC}, ${PWD}/src/opencv/)
	SUBMODULES += .opencv
endif


# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
OPENCVSOURCES = $(wildcard src/opencv/CMakeLists.txt)
SOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
BUILT_PROGRAMS = src/wally
TARGETS = ${SUBMODULES} ${BUILT_PROGRAMS}

all:   	$(TARGETS)

.htslib: $(HTSLIBSOURCES)
	if [ -r src/htslib/Makefile ]; then cd src/htslib && autoheader && autoconf && ./configure --disable-s3 --disable-gcs --disable-libcurl --disable-plugins && $(MAKE) && $(MAKE) lib-static && cd ../../ && touch .htslib; fi

.opencv: $(OPENCVSOURCES)
	if [ -f src/opencv/CMakeLists.txt ]; then cd src/opencv/ && mkdir build && cd build/ &&  cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=${PWD}/src/ocv -D WITH_GTK=ON .. &&  make -j 4 && make install && cd ../ && rm -rf build/ && cd ../../ && touch .opencv; fi

src/wally: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	if [ -r src/htslib/Makefile ]; then cd src/htslib && $(MAKE) clean; fi
	rm -f $(TARGETS) $(TARGETS:=.o) ${SUBMODULES}
	rm -rf src/ocv/ src/opencv/build/

distclean: clean
	rm -f ${BUILT_PROGRAMS}

.PHONY: clean distclean install all
