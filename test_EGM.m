%% setup

clear;
close all;
clc;

mex -setup:'C:\Program Files\MATLAB\R2016a\bin\win64\mexopts\mingw64_g++.xml'
mex cfuncs\mex_acon_grid.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -fopenmp -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
mex cfuncs\mex_E.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -fopenmp -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
mex cfuncs\mex_E_vec.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -fopenmp -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
mex cfuncs\mexUpperEnvelopeToCommon.cpp CXXFLAGS="$CXXFLAGS -fopenmp -std=c++11 -O3 -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
mex cfuncs\mexVFI_NLopt.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -fopenmp -march=native -ffast-math -fabi-version=0 -mavx2" cfuncs\libnlopt-0.lib C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;

%% run

par   = SetupPar_3d();
par.T = 20;

par.print = 1;
par.Nm    = 50;
par.T     = 2;

%profile on
sol = solve.egm(par);
%profile viewer