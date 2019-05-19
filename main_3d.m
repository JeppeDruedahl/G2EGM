%% setup

clc;
clear;
close all;

%% all solution figures

methods = {'egm','vfi'};

FIGVISIBLE  = 'off';
FIG         = 1;
EULER       = 0;

for imethod = 1:numel(methods)
    
    method = methods{imethod};
    
    % a. load data
    casefig = 'nonsmooth_3d';
    if strcmp(method,'egm')
        load('data\egm_Nm150_nonsmooth_3d.mat')
    else
        load('data\vfi_Nm150_nonsmooth_3d.mat')
    end
    
    % b. figures
    if FIG == 1
    figfolder = ['figures\' casefig];
    if isdir(figfolder) == 0
        mkdir(figfolder)
    end
    figfolder = ['figures\' casefig '\pdf\'];
    if isdir(figfolder) == 0
        mkdir(figfolder)
    end
   
    par = fun.solprep(par);
    for t = [1]
        
        sol(2,t).a = par.grid_m_nd + par.rk*par.grid_k_nd.*sol(2,t).l - sol(2,t).c - sol(2,t).d;
        sol(2,t).b = par.grid_n_nd + sol(2,t).d + fun.f_pens(sol(2,t).d,par);
        sol(2,t).q = (1-par.delta)*par.grid_k_nd + sol(2,t).l;

        figs = struct();
        par = fun.solprep(par);
        figs = funfig_3d.plot(figs,par,sol,par.T-1,floor(par.Nk/1.5),floor(par.Nn/3),method,{'c','d','l','v','a','b','q'});    
        fprintf('k = %g, n = %g\n',par.grid_k(floor(par.Nk/1.5)),par.grid_n(floor(par.Nn/3)));

        fun.printfig(figs.c.(method).mn,FIGVISIBLE,casefig)
        fun.printfig(figs.d.(method).mn,FIGVISIBLE,casefig)
        fun.printfig(figs.l.(method).mn,FIGVISIBLE,casefig)
        fun.printfig(figs.a.(method).mn,FIGVISIBLE,casefig)
        fun.printfig(figs.b.(method).mn,FIGVISIBLE,casefig)
        fun.printfig(figs.q.(method).mn,FIGVISIBLE,casefig)
        fun.printfig(figs.v.(method).mn,FIGVISIBLE,casefig)
 
        fun.printfig(figs.c.(method).mk,FIGVISIBLE,casefig)
        fun.printfig(figs.d.(method).mk,FIGVISIBLE,casefig)
        fun.printfig(figs.l.(method).mk,FIGVISIBLE,casefig)
        fun.printfig(figs.a.(method).mk,FIGVISIBLE,casefig)
        fun.printfig(figs.b.(method).mk,FIGVISIBLE,casefig)
        fun.printfig(figs.q.(method).mk,FIGVISIBLE,casefig)
        fun.printfig(figs.v.(method).mk,FIGVISIBLE,casefig)
   
    end
    end
    
    % c. euler errors
    if EULER == 1
        fun.simulate_euler(sol,100,par.T,par);
    end
    
end