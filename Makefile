# Use gcc-5 from MacPorts on OSX
OSXGCC5    := yes
# Use clang from MacPorts on OSX
# OSXMPCLANG := yes
#
# Can also just set CXX here:
# CXX := my-c++-compiler

# Vector report level for icc
VEC_REP := 1

# Define to build for AVX_512, the new mic (KNL) and latest generation Xeons.
# AVX_512 := 1

# Set to 1 to also do KNC mic build
# MIC := 1

# Compiler
ifdef INTEL_LICENSE_FILE
  CXX     := icc
  ifdef AVX_512
    VECHOST := -xMIC-AVX512 -vec-report=${VEC_REP}
  else
    VECHOST := -mavx -vec-report=${VEC_REP}
  endif
  VECMIC  := -mmic -vec-report=${VEC_REP}
else ifdef OSXGCC5
  CXX := c++-mp-5
  TBB_PREFIX := /opt/local
else ifdef OSXMPCLANG
  CXX := clang++-mp-3.9 -Wall -Wno-unknown-pragmas -Wno-missing-braces
  TBB_PREFIX := /opt/local
endif

### To disable vectorization set USER_CXXFLAGS := -no-simd -no-vec
# Setting only one of the above has little effect.
# Note, this also screws-up prefetching so it's a lousy deal.

# -opt-prefetch-distance=64,8
OPT      := -O3

CPPFLAGS := -I. ${USER_CPPFLAGS} ${DEFS}
CXXFLAGS := ${OPT} -openmp -std=gnu++0x ${USER_CXXFLAGS}

LDFLAGS  := ${USER_LDFLAGS}

TGTS     := t1 t2 t3 t4 mkFit/mkFit

EXES     := ${TGTS}

ifdef MIC
  EXES += $(addsuffix -mic, ${TGTS})
endif

# Was: all: ${EXES}
# t3 t4 and mkFit do now work yet after codas "cleanup"
all: t1

auto-matriplex:
	${MAKE} -C Matriplex auto && touch $@

auto:: auto-matriplex
.PHONY: auto

%.o: %.cxx *.h auto-matriplex
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VECHOST} -c -o $@ $<

%.om: %.cxx *.h auto-matriplex
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VECMIC} -c -o $@ $<


clean::
	rm -f ${EXES} *.o *.om mkFit/mkFit mkFit/mkFit-mic mkFit/*.o mkFit/*.om

distclean:: clean
	${MAKE} -C Matriplex clean
	rm -f auto-matriplex

echo:
	@echo CPPFLAGS = ${CPPFLAGS}
	@echo CXXFLAGS = ${CXXFLAGS}
	@echo LDFLAGS  = ${LDFLAGS}
	@echo TGTS     = ${TGTS}
	@echo EXES     = ${EXES}

################################################################

### mkFit test

MKFSRCS := $(wildcard mkFit/*.cc)
MKFHDRS := $(wildcard mkFit/*.h)

MKFOBJS     := $(MKFSRCS:.cc=.o)  Timing.o
MKFOBJS_MIC := $(MKFSRCS:.cc=.om) Timing.om

MKEXES   := mkFit/mkFit mkFit/mkFit-mic

### To run without validation
MK_HOST_CFLAGS := -DNO_ROOT
MK_HOST_LIBS   :=
### To run with validation (host only)
# MK_HOST_CFLAGS := -I${ROOTSYS}/include
# MK_HOST_LIBS   := -L${ROOTSYS}/lib -lCore -lRIO -lTree

auto-mkfit: mkFit/GenMPlexOps.pl
	(cd mkFit; ./GenMPlexOps.pl)
	touch auto-mkfit

distclean::
	rm -f auto-mkfit
	rm -f mkFit/*.ah

.PHONY: mkFit
mkFit: ${MKEXES}

mkFit/mkFit: auto auto-mkfit ${MKFOBJS}
	${CXX} ${CXXFLAGS} ${VECHOST} ${LDFLAGS} ${MK_HOST_LIBS} -o $@ ${MKFOBJS}

mkFit/mkFit-mic: auto auto-mkfit ${MKFOBJS_MIC}
	${CXX} ${CXXFLAGS} ${VECMIC}  ${LDFLAGS} -o $@ ${MKFOBJS_MIC}
	scp $@ mic0:

mkFit/%.o: mkFit/%.cc mkFit/*.h Matriplex/*
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VECHOST} -IMatriplex ${MK_HOST_CFLAGS} -c -o $@ $<

mkFit/%.om: mkFit/%.cc mkFit/*.h Matriplex/*
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VECMIC} -DNO_ROOT -IMatriplex -c -o $@ $<

### t1

T1_DEPS := t1 ArrayTest Timing

t1: $(addsuffix .o, ${T1_DEPS})
	${CXX} ${CXXFLAGS} ${VECHOST} ${LDFLAGS} -o $@ $^

t1-mic:  $(addsuffix .om, ${T1_DEPS})
	${CXX} ${CXXFLAGS} ${VECMIC} ${LDFLAGS} -o $@ $^
	scp $@ mic0:

run-t1:	t1 t1-mic
	./t1
	ssh mic0 ./t1-mic

### t2

T2_DEPS := t2 ArrayTest Timing

t2: $(addsuffix .o, ${T2_DEPS})
	${CXX} ${CXXFLAGS} ${VECHOST} ${LDFLAGS} -o $@ $^

t2-mic:  $(addsuffix .om, ${T2_DEPS})
	${CXX} ${CXXFLAGS} ${VECMIC} ${LDFLAGS} -o $@ $^
	scp $@ mic0:

run-t2:	t2 t2-mic
	./t2
	ssh mic0 ./t2-mic

### t3

T3_DEPS     := t3 MPlexTest Timing
T3_OBJS     := $(addsuffix .o,  ${T3_DEPS})
T3_OBJS_MIC := $(addsuffix .om, ${T3_DEPS})

MPlexTest.o MPlexTest.om: Matriplex/Matriplex.h Matriplex/MatriplexSym.h Matriplex/MatriplexVector.h

t3: auto ${T3_OBJS}
	${CXX} ${CXXFLAGS} ${VECHOST} ${LDFLAGS} -o $@ ${T3_OBJS}

t3-mic: auto ${T3_OBJS_MIC}
	${CXX} ${CXXFLAGS} ${VECMIC}  ${LDFLAGS} -o $@ ${T3_OBJS_MIC}
	scp $@ mic0:

run-t3:	t3 t3-mic
	./t3
	ssh mic0 ./t3-mic

### t4

T4_DEPS     := t4 MPlexTest
T4_OBJS     := $(addsuffix .o,  ${T4_DEPS})
T4_OBJS_MIC := $(addsuffix .om, ${T4_DEPS})

MPlexTest.o MPlexTest.om: Matriplex/Matriplex.h Matriplex/MatriplexSym.h Matriplex/MatriplexVector.h

t4: auto ${T4_OBJS}
	${CXX} ${CXXFLAGS} ${VECHOST} ${LDFLAGS} -o $@ ${T4_OBJS}

t4-mic: auto ${T4_OBJS_MIC}
	${CXX} ${CXXFLAGS} ${VECMIC}  ${LDFLAGS} -o $@ ${T4_OBJS_MIC}
	scp $@ mic0:

run-t4:	t4 t4-mic
	./t4
	ssh mic0 ./t4-mic

