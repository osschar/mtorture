# Request latest gcc!
# . /opt/rh/devtoolset-2/enable

# -mavx -msse4.2
VECHOST  := -mavx -vec-report=1
VECMIC   := -mmic -vec-report=1

VECOPTS  :=
#VECOPTS  := -no-simd -no-vec
# Just setting one of the above has little effect

CPPFLAGS := -I.
CXXFLAGS := -O3 -openmp -std=gnu++0x ${VECOPTS}

LDFLAGS  :=

all: t1

%.o: %.cxx *.h
	icc ${CPPFLAGS} ${CXXFLAGS} ${VECHOST} -c -o $@ $<


%.om: %.cxx *.h
	icc ${CPPFLAGS} ${CXXFLAGS} ${VECMIC} -c -o $@ $<

t1: t1.o ArrayTest.o Timing.o
	icc ${CXXFLAGS} ${VECHOST} ${LDFLAGS} -o $@ $^

t1-mic: t1.om ArrayTest.om Timing.om
	icc ${CXXFLAGS} ${VECMIC} ${LDFLAGS} -o $@ $^
	scp $@ mic0:

run-t1:	t1 t1-mic
	./t1
	ssh mic0 ./t1-mic

clean:
	rm -f t1 t1-mic *.o *.om
