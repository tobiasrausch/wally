SHELL := /bin/bash

DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
EBROOTHTSLIB ?= ${PWD}/src/htslib/
BLEND2DSRC ?= ${PWD}/src/blend2d/
BLEND2DLIB ?= ${BLEND2DSRC}/build/libblend2d.a

# Optional: wasm
EMCXX ?= em++
EMCONFIGURE ?= emconfigure
EMCMAKE ?= emcmake
EMMAKE ?= emmake
BOOST_VER ?= 1.83.0
BOOST_USCORE = $(subst .,_,${BOOST_VER})
WASMDIR = ${PWD}/wasm
WASMHTS = ${WASMDIR}/htslib/libhts.a
WASMB2D = ${BLEND2DSRC}/build-wasm/libblend2d.a
WASMBOOST = ${WASMDIR}/boost/libboost_wasm.a
NPROC = $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
EMCONF = ${WASMDIR}/emscripten_config
EMCACHE = ${WASMDIR}/emcache
export EM_CONFIG = ${EMCONF}

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
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz -llzma -lbz2 -ldeflate
else
	LDSTATIC =
	HTSSTATIC =
	LDFLAGS += -pthread -lhts -lz -llzma -lbz2 -Wl,-rpath,${EBROOTHTSLIB}
endif

# librt
ifeq ($(shell uname -s), Linux)
	LDFLAGS += -lrt
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
	if [ -f ${BLEND2DSRC}/CMakeLists.txt ]; then unset CXXFLAGS CFLAGS LDFLAGS; cd ${BLEND2DSRC} && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBLEND2D_STATIC=ON -DBLEND2D_NO_JIT=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON && cmake --build build --target blend2d -j${NPROC} && cd ../../ && touch .blend2d; fi

src/wally: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(BLEND2DLIB) $(LDFLAGS)

# Wasm
wasm: src/wally client/wally.js

${EMCONF}:
	mkdir -p ${WASMDIR}
	EM_CONFIG=${EMCONF} emcc --generate-config
	sed -i "/^FROZEN_CACHE/d" ${EMCONF}
	printf "\nFROZEN_CACHE = False\nCACHE = '${EMCACHE}'\n" >> ${EMCONF}
	echo 'int main(){return 0;}' > ${WASMDIR}/.warm.c
	EMCC_CORES=1 emcc -O2 -sUSE_ZLIB=1 ${WASMDIR}/.warm.c -o ${WASMDIR}/.warm_c.js
	EMCC_CORES=1 ${EMCXX} -O2 -fexceptions -sUSE_ZLIB=1 -xc++ ${WASMDIR}/.warm.c -lworkerfs.js -lnodefs.js -o ${WASMDIR}/.warm_cpp.js
	rm -f ${WASMDIR}/.warm.c ${WASMDIR}/.warm_c.* ${WASMDIR}/.warm_cpp.*

${WASMHTS}: ${EMCONF}
	mkdir -p ${WASMDIR}
	rsync -a --exclude='*.o' --exclude='*.a' --exclude='*.so*' --exclude='.git' src/htslib/ ${WASMDIR}/htslib/
	cd ${WASMDIR}/htslib && ${EMCONFIGURE} ./configure --host=wasm32-unknown-emscripten --disable-libcurl --disable-lzma --disable-bz2 --without-libdeflate CFLAGS="-O2 -sUSE_ZLIB=1" && ${EMMAKE} make -j${NPROC} libhts.a

${WASMB2D}: ${EMCONF}
	unset CXXFLAGS CFLAGS LDFLAGS; cd ${BLEND2DSRC} && ${EMCMAKE} cmake -S . -B build-wasm -DCMAKE_BUILD_TYPE=Release -DBLEND2D_STATIC=ON -DBLEND2D_NO_JIT=ON && ${EMMAKE} cmake --build build-wasm --target blend2d -j${NPROC}

${WASMBOOST}: ${EMCONF}
	mkdir -p ${WASMDIR}
	if [ ! -d ${WASMDIR}/boost ]; then curl -sSL -o ${WASMDIR}/boost.tar.gz https://archives.boost.io/release/${BOOST_VER}/source/boost_${BOOST_USCORE}.tar.gz && tar xzf ${WASMDIR}/boost.tar.gz -C ${WASMDIR} && mv ${WASMDIR}/boost_${BOOST_USCORE} ${WASMDIR}/boost && rm -f ${WASMDIR}/boost.tar.gz; fi
	cd ${WASMDIR}/boost && mkdir -p obj && export CF="-std=c++17 -O2 -I. -sUSE_ZLIB=1 -DBOOST_ALL_NO_LIB -DBOOST_FILESYSTEM_NO_CXX20_ATOMIC_REF -DNDEBUG -fvisibility=hidden" && ${EMCXX} $$CF -c libs/system/src/error_code.cpp -o obj/system_error_code.o && for f in codecvt_error_category directory exception operations path path_traits portability unique_path utf8_codecvt_facet; do ${EMCXX} $$CF -c libs/filesystem/src/$$f.cpp -o obj/fs_$$f.o || exit 1; done && ${EMCXX} $$CF -c libs/iostreams/src/zlib.cpp -o obj/io_zlib.o && ${EMCXX} $$CF -c libs/iostreams/src/gzip.cpp -o obj/io_gzip.o && emar rcs libboost_wasm.a obj/*.o

client/wally.js: src/wasm.cpp ${EMCONF} ${WASMHTS} ${WASMB2D} ${WASMBOOST}
	${EMCXX} -O2 -std=c++17 -fexceptions -DNDEBUG -DBOOST_FILESYSTEM_NO_CXX20_ATOMIC_REF -I ${WASMDIR}/htslib -I ${BLEND2DSRC} -I ${WASMDIR}/boost src/wasm.cpp ${WASMHTS} ${WASMB2D} ${WASMBOOST} -sUSE_ZLIB=1 -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=64MB -sMAXIMUM_MEMORY=2GB -sMODULARIZE=1 -sEXPORT_NAME=createWally -sFORCE_FILESYSTEM=1 -sEXPORTED_FUNCTIONS=_wally_region,_malloc,_free -sEXPORTED_RUNTIME_METHODS=FS,ccall,cwrap -lworkerfs.js -lnodefs.js -o client/wally.js

wasm-clean:
	rm -rf ${WASMDIR}/htslib ${WASMDIR}/boost ${WASMDIR}/emcache ${EMCONF} ${BLEND2DSRC}/build-wasm client/wally.js client/wally.wasm

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	if [ -r src/htslib/Makefile ]; then cd src/htslib && $(MAKE) clean; fi
	rm -f $(TARGETS) $(TARGETS:=.o) ${SUBMODULES}
	rm -rf ${BLEND2DSRC}/build/

distclean: wasm-clean clean
	rm -rf ${BUILT_PROGRAMS}

.PHONY: clean distclean install all wasm wasm-clean
