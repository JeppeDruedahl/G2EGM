# G2EGM

Code for [A General Endogenous Grid Method for Multi-Dimensional Models with Non-Convexities and Constraints](https://linkinghub.elsevier.com/retrieve/pii/S0165188916301920), 2017, *Journal of Economic Dynamics and Control*, 74.

# READ ME 

*This version of the code differ from the one used in the paper due to changes in MATLABs MEX setup.*

The code is written and tested for 64-bit Windows machines using MATLAB 2018b.
It requires MATLAB and a C++ compiler.

Tested using Microsoft Visual Studio 2017 Community Edition (free)
https://visualstudio.microsoft.com/downloads/

For NLopt see cfuncs\ReadMe\NLopt.txt
For vectorclass see cfuncs\ReadMe_vectorclass.txt.

# FOLDER STRUCTURE 

All MATLAB files are in the main folder.
The subfolder \cfuncs contains all the C-files.
The subfolder \data is used for storing solution data.
The subfolder \figures are used for storing printed figures.

# MAIN 

fun.m: contains various MATLAB functions.
egmfun.m contains MATLAB functions unique for EGM.
solve.m: contains MATLAB functions for solving.

## 2-dimensional model
SetupPar.m: function for setting the parameter struct.
accuracy.m: computes the speed and accuracy measures used in the paper for.
main.m: solves the model for selected specifications and prints figures. 
funfig.m: containt MATLAB functions for printing.

## 3-dimensional model
SetupPar_3d.m: function for setting the parameter struct.
accuracy_3d.m: computes the speed and accuracy measures used in the paper for.
main_3d.m: solves the model for selected specifications and prints figures. 
funfig_3d.m: containt MATLAB functions for printing.

## other
test_EGM.m: simple test file for running the algorithm.
illustrate_upperenvelope.m: create upper envelope illustrations.

## C++-functions
mexVFI.cpp: solution algorithm with VFI.
mexVFI_NLopt.cpp: solution algorithm with VFI and NLopt.
mexUpperEnvelopeToCommon.cpp: G2EGM upper envelope and interpolation to common grid.
mex_acon_grid.cpp: construction of acon grid for 3-dimensional model.
mesh_E.cpp: construction of w, wa and wb.
mesh_E_vec.cpp: construction of w, wa and wb using some vector operations.
base_funcs.c: various c++-functions.
mesh_interp.cpp: interpolation function.
HighResTimer_class.hpp: timing function.