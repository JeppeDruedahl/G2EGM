%% setup

clear;
close all;
clc;

%% settings

ONLY_SOL_EGM    = 0;
LOAD_SOL        = 0;

FIG        = 1;
FIGVISIBLE = 'on';
EXTRA_FIG  = 1;

cases      = {'nonsmooth','smooth'};
methodlist = {'egm','vfi'};

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
    par = SetupPar();
    par.print = 1;
    par.max_threads = 7;

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
    if ONLY_SOL_EGM == 1
        load(filename);        
        sols = solve.all(sols,{'egm'},par);         
    elseif LOAD_SOL == 1
        load(filename);        
    else
        sols = solve.all(sols,methodlist,par);
    end
    if LOAD_SOL == 0
        save(filename,'sols','-v7.3');
    end
    
    % d. figures: value and policy functions
    if FIG == 1
    figfolder = ['figures\' casenow];
    if isdir(figfolder) == 0
        mkdir(figfolder)
    end  
    figfolder = ['figures\' casenow '\pdf\'];
    if isdir(figfolder) == 0
        mkdir(figfolder)
    end      

    kvec = [5 19];
    I = kvec < par.T;
    for k = kvec(I)
        
        fprintf('figs, T-%d\n',k)
        figs = struct();
        
        % i. overview    
        t1 = tic;
        for j = 1:numel(methodlist)            
            method = methodlist{j};
            if strcmp(method,'egm')
                figs = funfig.overview(figs,par,sols.(method),k,method);
                fun.printfig(figs.overview.(method),FIGVISIBLE,casenow);       
                fun.printfig(figs.discrete.(method),FIGVISIBLE,casenow);    
            end
        end
        fprintf('  overview: %g secs\n',round(toc(t1)*10)/10);
        
    end
    
    kvec = [19];
    I = kvec < par.T;
    for k = kvec(I)
        
        % ii. comparisons
        t1 = tic;
        
        figs = funfig.compare(figs,par,sols.vfi,sols.egm,k,{'c','d','v'});      
        fun.printfig(figs.c.egm,FIGVISIBLE,casenow)
        fun.printfig(figs.d.egm,FIGVISIBLE,casenow)
        fun.printfig(figs.v.egm,FIGVISIBLE,casenow)       
       
            if EXTRA_FIG == 1        
            % extra
            figs = funfig.compare(figs,par,sols.vfi,sols.egm,k,{'a','b'});       
            fun.printfig(figs.a.egm,FIGVISIBLE,casenow)
            fun.printfig(figs.b.egm,FIGVISIBLE,casenow)         
            end % extrafig
            
        fprintf('  comparison: %g secs\n',round(toc(t1)*10)/10);

        close all;
        clearvars figs;
        
    end % kvec
end % FIG    
    
end