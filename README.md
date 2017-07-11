# mtorture -- Machine Torture


## Determining clock speed

### Linux

```
lscpu
# CPU flags, look for sse4.2, avx, avx2, avx512
cat /proc/cpuinfo  | grep flags | tail -n1
su -
echo performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
cat /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_cur_freq
cat /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq
# Enter frequency into Timing.cxx (note, this is in kHz, multiply by 1000)
# Run the tests ...
echo ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
```

### OSX

```
sysctl machdep.cpu.features
sysctl hw.cpufrequency_min hw.cpufrequency hw.cpufrequency_max
# Enter the max value into Timing.cxx (note, this is in Hz, copy as is)
```


## External documentation of interest

- Lists of instruction latencies and throughput for Intel, AMD and VIA CPUs:
http://www.agner.org/optimize/instruction_tables.pdf

- Intel intrinsics guide:
https://software.intel.com/sites/landingpage/IntrinsicsGuide/


## Running tests

### Makefile

Variables:
- `OPT` -- set optimiztion level, default is `-O3`
- `USER_CPPFLAGS` -- set extra `cpp` flags, especially `-D` for compilation
- `USER_CXXFLAGS` -- whatever you want to pass to compiler
- `USER_LDFLAGS`  -- whatever you want to pass to linker

There is a special target *echo* that will print out the relevant varables, e.g:
```
  make CXX:=my-c++ OPT:=-O2 USER_CPPFLAGS:="-DTEST_FUNC=sum3_cube -DFOO=bar" echo
```


### Defines used in tests

- `TEST_FUNC` -- the default function that gets called on the test object

To use this from a script, you'd have to do something like this:
```
  make USER_CPPFLAGS:="-DTEST_FUNC=sum3_cube" -W t1.cxx t1 t1-mic
```

### Global variables and environment variables

Format: `global type and name | environment name | default value`, see also `common.h`.

- `double g_test_duration | TEST_DURATION | 1.0  ` -- duration of each test in seconds
- `double g_pre_test_frac | PRE_TEST_FRAC | 0.01 ` -- fraction of the above to use for loop calibration

- `int    g_n_vec_min     | N_VEC_MIN     | 8    ` -- starting vector size for the test
- `int    g_n_vec_max     | N_VEC_MAX     | 64*1024*1024 ` -- ending vector size
