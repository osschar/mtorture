# Request latest gcc!
# . /opt/rh/devtoolset-2/enable

# -mavx -msse4.2
VECHOST  := -mavx -vec-report=1
VECMIC   := -mmic -vec-report=1

### To disable vectorization set USER_CXXFLAGS := -no-simd -no-vec
# Setting only one of the above has little effect.
# Note, this also screws-up prefetching so it's a lousy deal.

OPT      := -O3

CPPFLAGS := -I. ${USER_CPPFLAGS}
CXXFLAGS := ${OPT} -openmp -std=gnu++0x ${USER_CXXFLAGS}

LDFLAGS  := ${USER_LDFLAGS}

all: t1 t1-mic t2 t2-mic

%.o: %.cxx *.h
	icc ${CPPFLAGS} ${CXXFLAGS} ${VECHOST} -c -o $@ $<

%.om: %.cxx *.h
	icc ${CPPFLAGS} ${CXXFLAGS} ${VECMIC} -c -o $@ $<


clean:
	rm -f t1 t1-mic *.o *.om

echo:
	@echo CPPFLAGS = ${CPPFLAGS}
	@echo CXXFLAGS = ${CXXFLAGS}
	@echo LDFLAGS  = ${LDFLAGS}


### t1

T1_DEPS := t1 ArrayTest Timing

t1: $(addsuffix .o, ${T1_DEPS})
	icc ${CXXFLAGS} ${VECHOST} ${LDFLAGS} -o $@ $^

t1-mic:  $(addsuffix .om, ${T1_DEPS})
	icc ${CXXFLAGS} ${VECMIC} ${LDFLAGS} -o $@ $^
	scp $@ mic0:

run-t1:	t1 t1-mic
	./t1
	ssh mic0 ./t1-mic

### t2

T2_DEPS := t2 ArrayTest Timing

t2: $(addsuffix .o, ${T2_DEPS})
	icc ${CXXFLAGS} ${VECHOST} ${LDFLAGS} -o $@ $^

t2-mic:  $(addsuffix .om, ${T2_DEPS})
	icc ${CXXFLAGS} ${VECMIC} ${LDFLAGS} -o $@ $^
	scp $@ mic0:

run-t2:	t2 t2-mic
	./t2
	ssh mic0 ./t2-mic

