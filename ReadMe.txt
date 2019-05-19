#############
## READ ME ##
#############

The code is written and tested for 64-bit Windows machines using MATLAB 2016a.
It requires MATLAB and a C++ compiler.

Best behavior is with GCC 4.9.2.
MATLAB installation: Apps -> Get more Apps -> MinGW-w64 C/C++ COmpiler from TDM-GCC.
For parallization also install separately: http://tdm-gcc.tdragon.net/ (with OpenMP)

For NLopt see cfuncs\ReadMe\NLopt.txt
For vectorclass see cfuncs\ReadMe_vectorclass.txt.


######################
## FOLDER STRUCTURE ##
######################

All MATLAB files are in the main folder.
The subfolder \cfuncs contains all the C-files.
The subfolder \data is used for storing solution data.
The subfolder \figures are used for storing printed figures.


##########
## MAIN ##
##########

fun.m: contains various MATLAB functions.
egmfun.m contains MATLAB functions unique for EGM.
solve.m: contains MATLAB functions for solving.

# 2-dimensional model
SetupPar.m: function for setting the parameter struct.
accuracy.m: computes the speed and accuracy measures used in the paper for.
main.m: solves the model for selected specifications and prints figures. 
funfig.m: containt MATLAB functions for printing.

# 3-dimensional model
SetupPar_3d.m: function for setting the parameter struct.
accuracy_3d.m: computes the speed and accuracy measures used in the paper for.
main_3d.m: solves the model for selected specifications and prints figures. 
funfig_3d.m: containt MATLAB functions for printing.

# other
test_EGM.m: simple test file for running the algorithm.
illustrate_upperenvelope.m: create upper envelope illustrations.

# c++-functions
mexVFI.cpp: solution algorithm with VFI.
mexVFI_NLopt.cpp: solution algorithm with VFI and NLopt.
mexUpperEnvelopeToCommon.cpp: G2EGM upper envelope and interpolation to common grid.
mex_acon_grid.cpp: construction of acon grid for 3-dimensional model.
mesh_E.cpp: construction of w, wa and wb.
mesh_E_vec.cpp: construction of w, wa and wb using some vector operations.
base_funcs.c: various c++-functions.
mesh_interp.cpp: interpolation function.
HighResTimer_class.hpp: timing function.