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
	if [ -f src/opencv/CMakeLists.txt ]; then cd src/opencv/ && mkdir build && cd build/ &&  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${PWD}/src/ocv -DBUILD_SHARED_LIBS=${CVSHARED} -DOPENCV_GENERATE_PKGCONFIG=ON -DENABLE_PIC=FALSE -DBUILD_JAVA=OFF -DWITH_1394=OFF -DWITH_ADE=OFF -DWITH_VTK=OFF -DWITH_EIGEN=OFF -DWITH_FFMPEG=OFF -DWITH_GSTREAMER=OFF -DWITH_GTK=${CVSHARED} -DWITH_GTK_2_X=OFF -DWITH_IPP=OFF -DWITH_JASPER=OFF -DWITH_JPEG=OFF -DWITH_WEBP=OFF -DWITH_OPENEXR=OFF -DWITH_OPENGL=OFF -DWITH_OPENVX=OFF -DWITH_OPENNI=OFF -DWITH_OPENNI2=OFF -DWITH_PNG=ON -DBUILD_PNG=ON -DWITH_TBB=OFF -DWITH_TIFF=OFF -DWITH_V4L=OFF -DWITH_OPENCL=OFF -DWITH_OPENCL_SVM=OFF -DWITH_OPENCLAMDFFT=OFF -DWITH_OPENCLAMDBLAS=OFF -DWITH_GPHOTO2=OFF -DWITH_LAPACK=OFF -DWITH_ITT=OFF -DWITH_QUIRC=OFF -DCV_TRACE=OFF -DBUILD_OPENJPEG=ON -DWITH_JPEG=ON -DBUILD_ZLIB=ON -DBUILD_opencv_dnn=OFF -DBUILD_opencv_highgui=ON -DBUILD_opencv_ml=OFF -DBUILD_opencv_objdetect=OFF -DBUILD_opencv_photo=OFF -DBUILD_opencv_python3=OFF -DBUILD_opencv_shape=OFF -DBUILD_opencv_stitching=OFF -DBUILD_opencv_superres=OFF -DBUILD_opencv_videoio=OFF -DBUILD_opencv_videostab=OFF -DBUILD_opencv_xfeatures2d=OFF -DBUILD_opencv_bgsegm=OFF -DBUILD_opencv_bioinspired=OFF -DBUILD_opencv_fuzzy=OFF -DBUILD_opencv_hfs=OFF -DBUILD_opencv_img_hash=OFF -DBUILD_opencv_intensity_transform=OFF -DBUILD_opencv_line_descriptor=OFF -DBUILD_opencv_optflow=OFF -DBUILD_opencv_phase_unwrapping=OFF -DBUILD_opencv_plot=OFF -DBUILD_opencv_rapid=OFF -DBUILD_opencv_reg=OFF -DBUILD_opencv_rgbd=OFF -DBUILD_opencv_saliency=OFF -DBUILD_opencv_stereo=OFF -DBUILD_opencv_structured_light=OFF -DBUILD_opencv_surface_matching=OFF -DBUILD_opencv_tracking=OFF -DBUILD_opencv_ximgproc=OFF -DBUILD_opencv_xobjdetect=OFF -DBUILD_opencv_xphoto=OFF -DBUILD_opencv_calib3d=ON -DBUILD_opencv_core=ON -DBUILD_opencv_features2d=ON -DBUILD_opencv_flann=ON -DBUILD_opencv_imgcodecs=ON -DBUILD_opencv_video=ON -DBUILD_opencv_apps=OFF -DBUILD_DOCS=OFF -DBUILD_EXAMPLES=OFF -DBUILD_IPP_IW=OFF -DBUILD_PACKAGE=OFF -DBUILD_PERF_TESTS=OFF -DBUILD_TESTS=OFF -DBUILD_WITH_DEBUG_INFO=OFF -DWITH_PTHREADS_PF=OFF -DCV_ENABLE_INTRINSICS=OFF .. && make -j 8 && make install && cd ../ && rm -rf build/ && cd ../../ && touch .opencv; fi

src/wally: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $(shell export PKG_CONFIG_PATH=${PWD}/src/ocv/lib/pkgconfig/:${PWD}/src/ocv/lib64/pkgconfig/:${PKG_CONFIG_PATH} && pkg-config --cflags opencv4) $@.cpp -o $@ $(shell export PKG_CONFIG_PATH=${PWD}/src/ocv/lib/pkgconfig/:${PWD}/src/ocv/lib64/pkgconfig/:${PKG_CONFIG_PATH} && pkg-config --libs ${LDSTATIC} opencv4) $(LDFLAGS)

.emsdk:
	cd src/ && git clone https://github.com/emscripten-core/emsdk && cd emsdk && git pull && ./emsdk install latest && ./emsdk activate latest && cd ../../ && touch .emsdk

.opencvjs: .emsdk $(OPENCVSOURCES)
	if [ -f src/opencv/CMakeLists.txt ]; then cd src/ && ./emsdk/emsdk activate latest && source "${PWD}/src/emsdk/emsdk_env.sh" && EMSCRIPTEN="${PWD}/src/emsdk/upstream/emscripten" python3 opencv/platforms/js/build_js.py opencv_js --build_wasm && cd ../ && touch .opencvjs; fi

html/image.js: .emsdk .opencvjs $(SOURCES)
	mkdir -p html/ && ./src/emsdk/emsdk activate latest && source "${PWD}/src/emsdk/emsdk_env.sh" && em++ -Isrc/ocv/include/opencv4/ src/image.cpp -o html/image.js -Lsrc/opencv_js/lib/ -lopencv_core -lopencv_imgproc -O3 -s NO_EXIT_RUNTIME=1 -s "EXPORTED_RUNTIME_METHODS=['ccall']" -s ASSERTIONS=1  --bind

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	if [ -r src/htslib/Makefile ]; then cd src/htslib && $(MAKE) clean; fi
	rm -f $(TARGETS) $(TARGETS:=.o) ${SUBMODULES} .emsdk .opencvjs
	rm -rf src/ocv/ src/opencv/build/ src/emsdk/ src/opencv_js/ html/

distclean: clean
	rm -f ${BUILT_PROGRAMS}

.PHONY: clean distclean install all
