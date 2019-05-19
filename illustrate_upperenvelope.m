%% cleaning

clc;
clear;
close all;

mex -setup:'C:\Program Files\MATLAB\R2016a\bin\win64\mexopts\mingw64_g++.xml'
mex cfuncs\mex_acon_grid.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -fopenmp -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
mex cfuncs\mex_E.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -fopenmp -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
mex cfuncs\mex_E_vec.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -fopenmp -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
mex cfuncs\mexUpperEnvelopeToCommon.cpp CXXFLAGS="$CXXFLAGS -fopenmp -std=c++11 -O3 -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
mex cfuncs\mexVFI_NLopt.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -fopenmp -march=native -ffast-math -fabi-version=0 -mavx2" cfuncs\libnlopt-0.lib C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;

%% settings

FIGVISIBLE = 'off';

cases      = {'nonsmooth'};
methodlist = {'egm'};

%% run all

for i = 1:numel(cases);
    
    % a. setup    
    sols       = struct;
    casenow    = cases{i};
    filename   = ['data\' casenow '.mat']; 
    simfilename = ['data\' casenow '_sim.mat']; 
    fprintf([casenow '\n']);
    
    % b. parametrization
    par               = SetupPar();
    par.Nm            = 200;
    par.save_segments = 1; 

    if strcmp(casenow,'nonsmooth')
        par.var_eta = 0;
        par.sigma   = 0;
    elseif strcmp(casenow,'smooth')
        par.var_eta = 0.1;        
        par.sigma   = 0.1;        
    else
        error('unknown case');
    end
    
    % c. solve (or load)
    sols = solve.all(sols,methodlist,par);        
    save(filename,'sols','-v7.3');
    
    % d. figures: value and policy functions
    figfolder = ['figures\' casenow];
    if isdir(figfolder) == 0
        mkdir(figfolder)
    end  
    figfolder = ['figures\' casenow '\pdf\'];
    if isdir(figfolder) == 0
        mkdir(figfolder)
    end      

    kvec = 5;
    I = kvec < par.T;
    for k = kvec(I)
        
        fprintf('figs, T-%d\n',k)
        figs = struct();
        
        % segment: irregular grids
        figs = funfig.segments(figs,par,sols.egm,sols.egm,k,'clean',{'acon'},{});
        figs = funfig.segments(figs,par,sols.egm,sols.egm,k,'final',{'acon'},{});
        figs = funfig.segments(figs,par,sols.egm,sols.egm,k,'final_opt',{'acon'},{});
        fun.printfig(figs.grid.acon.clean,FIGVISIBLE,casenow); 
        fun.printfig(figs.grid.acon.final,FIGVISIBLE,casenow); 
        fun.printfig(figs.grid.acon.final_opt,FIGVISIBLE,casenow); 

    end % kvec
        
end