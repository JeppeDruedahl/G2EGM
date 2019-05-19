mex cfuncs\mex_acon_grid.cpp COMPFLAGS="$COMPFLAGS /Ox /openmp /arch:AVX2"
mex cfuncs\mex_E.cpp COMPFLAGS="$COMPFLAGS /Ox /openmp /arch:AVX2"
mex cfuncs\mex_E_vec.cpp COMPFLAGS="$COMPFLAGS /Ox /openmp /arch:AVX2"
mex cfuncs\mexUpperEnvelopeToCommon.cpp COMPFLAGS="$COMPFLAGS /Ox /openmp /arch:AVX2"
mex cfuncs\mexVFI_NLopt.cpp COMPFLAGS="$COMPFLAGS /Ox /openmp /arch:AVX2" cfuncs/nlopt-2.4.2-dll64/libnlopt-0.lib
mex cfuncs\mexVFI.cpp COMPFLAGS="$COMPFLAGS /Ox /openmp /arch:AVX2"

run accuracy.m
run main.m
run accuracy_3d.m
run main_3d.m
run illustrate_upperenvelope.m