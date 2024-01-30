SHELL := /bin/bash

DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
EBROOTHTSLIB ?= ${PWD}/src/htslib/
OPENCVSRC ?= ${PWD}/src/opencv/
OPENCV ?= ${PWD}/src/ocv

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir ?= $(exec_prefix)/bin

# Flags
CXX=g++
CXXFLAGS += -std=c++11 -isystem ${EBROOTHTSLIB} -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS -fno-strict-aliasing -fpermissive
LDFLAGS += -L${EBROOTHTSLIB} -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time

# Flags for static compile
ifeq (${STATIC}, 1)
	LDSTATIC = --static
	HTSSTATIC = --disable-libcurl
	CVSHARED = "OFF"
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz -llzma -lbz2 -ldeflate
else
	LDSTATIC =
	HTSSTATIC =
	CVSHARED = "ON"
	LDFLAGS += -pthread -lhts -lz -llzma -lbz2 -Wl,-rpath,${EBROOTHTSLIB} -Wl,-rpath,${OPENCV}/lib/ -Wl,-rpath,${OPENCV}/lib64/
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
	if [ -r src/htslib/Makefile ]; then cd src/htslib && autoreconf -i && ./configure ${HTSSTATIC} --disable-plugins && $(MAKE) && $(MAKE) lib-static && cd ../../ && touch .htslib; fi

.opencv: $(OPENCVSOURCES)
	if [ -f src/opencv/CMakeLists.txt ]; then cd src/opencv/ && mkdir build && cd build/ &&  cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=${PWD}/src/ocv -D BUILD_SHARED_LIBS=${CVSHARED} -DOPENCV_GENERATE_PKGCONFIG=ON -D BUILD_ZLIB=ON -D BUILD_PNG=ON -D BUILD_OPENJPEG=ON -D WITH_OPENEXR=OFF -D WITH_JPEG=OFF -D WITH_JASPER=OFF -D WITH_TIFF=OFF -D WITH_WEBP=OFF -D WITH_OPENCL=OFF -D WITH_GTK=${CVSHARED} -D WITH_FFMPEG=OFF -D WITH_1394=OFF -D WITH_IPP=OFF -D BUILD_TESTS=OFF -D BUILD_PERF_TESTS=OFF -D BUILD_opencv_apps=OFF .. &&  make -j 4 && make install && cd ../ && rm -rf build/ && cd ../../ && touch .opencv; fi

src/wally: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $(shell export PKG_CONFIG_PATH=${PWD}/src/ocv/lib/pkgconfig/:${PWD}/src/ocv/lib64/pkgconfig/:${PKG_CONFIG_PATH} && pkg-config --cflags opencv4) $@.cpp -o $@ $(shell export PKG_CONFIG_PATH=${PWD}/src/ocv/lib/pkgconfig/:${PWD}/src/ocv/lib64/pkgconfig/:${PKG_CONFIG_PATH} && pkg-config --libs ${LDSTATIC} opencv4) $(LDFLAGS)

.emsdk:
	cd src/ && git clone https://github.com/emscripten-core/emsdk && cd emsdk && git pull && ./emsdk install latest && ./emsdk activate latest && cd ../../ && touch .emsdk

.opencvjs: .emsdk $(OPENCVSOURCES)
	if [ -f src/opencv/CMakeLists.txt ]; then cd src/ && ./emsdk/emsdk activate latest && source "${PWD}/src/emsdk/emsdk_env.sh" && EMSCRIPTEN="${PWD}/src/emsdk/upstream/emscripten" python3 opencv/platforms/js/build_js.py opencv_js --build_wasm && cd ../ && touch .opencvjs; fi

.htslibjs: .emsdk .htslib
	if [ -r src/htslib/Makefile ]; then cp -r src/htslib/ src/htslib_js/ && cd src/htslib_js/ && make clean && make distclean && ../emsdk/emsdk activate latest && source "${PWD}/src/emsdk/emsdk_env.sh" && autoheader && autoconf && EMSCRIPTEN="${PWD}/src/emsdk/upstream/emscripten" emconfigure ./configure CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1" --disable-lzma --disable-plugins && emmake make CFLAGS="-O2 -s USE_ZLIB=1 -s USE_BZIP2=1" LDFLAGS="-O2 -s ERROR_ON_UNDEFINED_SYMBOLS=0" && cd ../../ && touch .htslibjs; fi

html/image.js: .emsdk .opencvjs .htslibjs $(SOURCES)
	mkdir -p html/ && ./src/emsdk/emsdk activate latest && source "${PWD}/src/emsdk/emsdk_env.sh" && em++ -Isrc/ocv/include/opencv4/ -Isrc/htslib_js/ src/image.cpp -o html/image.js -Lsrc/opencv_js/lib/ -Lsrc/htslib_js/ -lopencv_core -lopencv_imgproc -lhts -lz -lbz2 -O3 -s USE_BOOST_HEADERS=1 -s NO_EXIT_RUNTIME=1 -s "EXPORTED_RUNTIME_METHODS=['ccall']" -s ASSERTIONS=1  --bind

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	if [ -r src/htslib/Makefile ]; then cd src/htslib && $(MAKE) clean; fi
	rm -f $(TARGETS) $(TARGETS:=.o) ${SUBMODULES} .emsdk .opencvjs html/*.js html/*.wasm
	rm -rf src/ocv/ src/opencv/build/ src/emsdk/ src/opencv_js/

distclean: clean
	rm -f ${BUILT_PROGRAMS}

.PHONY: clean distclean install all
