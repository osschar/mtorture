# mtorture -- Machine Torture


## Running tests


### Makefile

Variables:
- `OPT` -- set optimiztion level, default is `-O3`
- `USER_CPPFLAGS` -- set extra `cpp` flags, especially `-D` for compilation
- `USER_CXXFLAGS` -- whatever you want to pass to compiler
- `USER_LDFLAGS`  -- whatever you want to pass to linker

There is a special target *echo* that will print out the relevant varables, e.g:
```
  make OPT:=-O2 USER_CPPFLAGS:="-DTEST_FUNC=sum3_cube -DFOO=bar" echo
```


### Defines used in tests

- `TEST_FUNC` -- the function that gets called on the test object

To use this from a script, you'd have to do something like this:
```
  make USER_CPPFLAGS:="-DTEST_FUNC=sum3_cube" -W t1.cxx t1 t1-mic
```

### Global variables and environment variables

Format: `global type and name | environment name | default value`, see also `common.h`.

- `double g_test_duration | TEST_DURATION | 1.0  ` -- duration of each test in seconds
- `double g_pre_test_frac | PRE_TEST_FRAC | 0.01 ` -- fraction of the above to
  use for loop calibration

- `int    g_n_vec_min     | N_VEC_MIN     | 8    ` -- starting vector size for the test
- `int    g_n_vec_max     | N_VEC_MAX     | 64*1024*1024 ` -- ending vector size
