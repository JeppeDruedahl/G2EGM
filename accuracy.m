%% setup

clc;
clear;
close all;

casenow = 'accuracy';
figfolder = ['figures\' casenow];
if isdir(figfolder) == 0
    mkdir(figfolder)
end
figfolder = ['figures\' casenow '\pdf\'];
if isdir(figfolder) == 0
    mkdir(figfolder)
end

MEX_GCC  = 1;
LOAD_MAX = 1;
LOAD_VFI = 0;
LOAD_EGM = 0;

MAX_THREADS = 1;
T = 20;

cases = {'nonsmooth','smooth'};
Nms   = [600,500,400,300,200,100];
simN  = 100; % number of indivudals simulated to calculate Euler errors.

N_guess_vfi_max = 400;
Nm_max          = 800; 

%% 1. maximum solution

if MEX_GCC == 1
    mex -setup:'C:\Program Files\MATLAB\R2016a\bin\win64\mexopts\mingw64_g++.xml'
    mex cfuncs\mex_acon_grid.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
    mex cfuncs\mex_E.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
    mex cfuncs\mex_E_vec.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
    mex cfuncs\mexUpperEnvelopeToCommon.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -march=native -ffast-math -fabi-version=0 -mavx2" C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
    mex cfuncs\mexVFI_NLopt.cpp CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -march=native -ffast-math -fabi-version=0 -mavx2" cfuncs\libnlopt-0.lib C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a;
else
    mex cfuncs\mexVFI.cpp COMPFLAGS="/openmp $COMPFLAGS";
    mex cfuncs\mex_E_vec.cpp COMPFLAGS="/openmp $COMPFLAGS";
    mex cfuncs\mex_E.cpp COMPFLAGS="/openmp $COMPFLAGS";
end

true_v = cell(numel(cases));
I      = cell(numel(cases));
for icase = 1:numel(cases);
        
    par_max             = SetupPar();
    par_max.Nm          = Nm_max;
    par_max.N_guess_vfi = N_guess_vfi_max;
    par_max.do_NLopt    = 0;
    par_max.max_threads = 16;
    par_max.print       = 1;
    
    name = cases{icase};
    if strcmp(name,'nonsmooth') == 1
        % non-smooth
    else
        % smooth
        par_max.var_eta = 0.1;
        par_max.sigma   = 0.1;
    end
    
    if LOAD_MAX == 0
        t1 = tic;
        sol_max = solve.vfi(par_max);
        par_max.time = toc(t1);    
        save(['data\true_max_' name '.mat'],'sol_max','par_max','-v7.3');    
    else                 
        load(['data\true_max_' name '.mat']);        
    end  
    
    par_max = fun.solprep(par_max);
    Im      = par_max.grid_m < par_max.fig_max_m;
    In      = par_max.grid_n < par_max.fig_max_n;
    [m, n]  = ndgrid(par_max.grid_m(Im),par_max.grid_n(In));

    interp = griddedInterpolant(par_max.grid_m_nd,par_max.grid_n_nd,sol_max(2,1).v);
    true_v{icase} = fun.trans_inv(interp(m,n),par_max);

    interp = griddedInterpolant(par_max.grid_m_nd,par_max.grid_n_nd,sol_max(2,1).d);
    true_d = interp(m,n);

    a = par_max.grid_m_nd - sol_max(2,1).c - sol_max(2,1).d;
    interp = griddedInterpolant(par_max.grid_m_nd,par_max.grid_n_nd,a);
    true_a = interp(m,n);

    I{icase} = true_d > 0 | true_a > 1e-8;

end


%% 2. VFI loop

errors_vfi = NaN(numel(Nms),numel(cases));
time_vfi   = NaN(numel(Nms),numel(cases));
euler_vfi  = NaN(numel(Nms),numel(cases));

for icase = 1:numel(cases)
for i = 1:numel(Nms)

    % a. settings
    par = SetupPar();
    par.Nm            = Nms(i);
    par.T             = T;
    par.max_threads   = MAX_THREADS;

    name = cases{icase};
    if strcmp(name,'nonsmooth') == 1
        % non-smooth
    else
        % smooth
        par.var_eta = 0.1;
        par.sigma   = 0.1;
    end    
    fprintf(['\n\n' name ': VFI, Nm = %d\n'],par.Nm);
    
    % b. solve or load
    if LOAD_VFI == 0
        t1 = tic;
        sol = solve.vfi(par);
        time_vfi(i,icase) = toc(t1);
        par.time = time_vfi(i,icase);
        save(['data\vfi_Nm' num2str(par.Nm) '_' name '.mat'],'sol','par','-v7.3');        
    else        
        load(['data\vfi_Nm' num2str(par.Nm) '_' name '.mat']);
        time_vfi(i,icase) = par.time;        
    end

    % c. errors
    par       = fun.solprep(par);
    interp    = griddedInterpolant({par.grid_m,par.grid_n},sol(2,1).v);
    now_v     = fun.trans_inv(interp(m,n),par);
    errorvec  = fun.vec(abs((now_v - true_v{icase})./true_v{icase}));
    errors_vfi(i,icase) = mean(fun.vec(errorvec(I{icase})));
    
    % d. Euler errors 
    sim                = fun.simulate_euler(sol,simN,par.T,par);
    euler_vfi(i,icase) = nanmean(-log10( abs(sim.euler_work(:)./sim.c(:)) + 1.0e-16));

end
end

%% 3. EGM loop
                               
errors_egm = NaN(numel(Nms),numel(cases));
time_egm   = NaN(numel(Nms),numel(cases));
euler_egm  = NaN(numel(Nms),numel(cases));

for icase = 1:numel(cases)
for i = 1:numel(Nms)
    
    % a. setup
    par             = SetupPar();    
    par.Nm          = Nms(i);
    par.T           = T;
    par.max_threads = MAX_THREADS;
    
    name = cases{icase};
    if strcmp(name,'nonsmooth') == 1
        % non-smooth
    else
        % smooth
        par.var_eta = 0.1;
        par.sigma   = 0.1;
    end
    fprintf(['\n\n' name ': Nm = %d\n'],par.Nm);
    
    % b. sols
    if LOAD_EGM == 0        
        t1                = tic;
        [sol, par]        = solve.egm(par);
        time_egm(i,icase) = toc(t1);
        par.time = time_egm(i,icase);
        save(['data\egm_Nm' num2str(par.Nm) '_' name '.mat'],'sol','par','-v7.3');
    else    
        load(['data\egm_Nm' num2str(par.Nm) '_' name '.mat'])
        time_egm(i,icase) = par.time;
    end
    
    % c. errors
    par           = fun.solprep(par);
    interp        = griddedInterpolant({par.grid_m,par.grid_n},sol(2,1).v);
    now_v         = fun.trans_inv(interp(m,n),par);
    errorvec      = fun.vec(abs((now_v - true_v{icase})./true_v{icase}));
    errors_egm(i,icase) = mean(fun.vec(errorvec(I{icase})));
    
    % c. Euler errors 
    sim                = fun.simulate_euler(sol,simN,par.T,par);
    euler_egm(i,icase) = nanmean(-log10( abs(sim.euler_work(:)./sim.c(:)) + 1.0e-16));

end
end

%% 4. Figures

for icase = 1:numel(cases);
    
    leg = cell(2,1);
    leg{1} = 'VFI';
    leg{2} = 'G$^2$EGM'; 
    
    % a. accuracy
    name = [cases{icase} '_accuracy'];      
    fig = figure('Name',name);
    fig.Color = [1 1 1];
    fig.Visible = 'off';
    ax = axes;    
    ax.FontSize = par.fontsize_small;
    ax.TickLabelInterpreter = 'latex';
    box(ax,'on');    
    semilogy(ax,Nms,errors_vfi(:,icase),'.--','Color','red','linewidth',2.5,'MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','black')
    set(gca,'FontSize',par.fontsize_small)
    hold on;
    semilogy(ax,Nms,errors_egm(:,icase),'.-','Color','black','linewidth',2.5,'MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','black')    
    fun.mylegend(par,leg,'NorthWest');
    xlim([100 600])
    ax.XTick = 100:50:600;
    ylim([5*10^(-7) 5*10^(-4)])
    ax.YTick = [10^(-6) 10^(-5) 10^(-4)];
    ax.FontSize = par.fontsize_small;
    grid on
    xlabel('nodes, $\#_m$','Interpreter','latex','FontSize', par.fontsize_small);
    ylabel('mean abs. rel. error','Interpreter','latex','FontSize', par.fontsize_small);
    set(gca,'FontSize',par.fontsize_small)
    fun.printfig(fig,'off',casenow);

    % b. time
    name = [cases{icase} '_time'];
    fig = figure('Name',name);
    fig.Color = [1 1 1];
    fig.Visible = 'off';
    ax = axes;    
    ax.FontSize = par.fontsize_small;
    ax.TickLabelInterpreter = 'latex';
    box(ax,'on');    
    plot(Nms,time_vfi(:,icase)./60,'.--','Color','red','linewidth',2.5,'MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','black')
    hold on
    set(gca,'FontSize',par.fontsize_small)
    plot(Nms,time_egm(:,icase)./60,'.-','Color','black','linewidth',2.5,'MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','black')
    fun.mylegend(par,leg,'NorthWest');
    hold off
    xlim([100 600])
    ax.XTick = 100:50:600;    
    ylim([0 18])
    ax.YTick = 0:3:18;
    ax.FontSize = par.fontsize_small;
    grid on
    xlabel('nodes, $\#_m$','Interpreter','latex','FontSize', par.fontsize_small);
    ylabel('minutes','Interpreter','latex','FontSize', par.fontsize_small);
    ax.FontSize = par.fontsize_small;
    set(gca,'FontSize',par.fontsize_small)
    fun.printfig(fig,'off',casenow);
        
    % c. time_rel
    name = [cases{icase} '_time_rel'];
    [figs.name, ax] = fun.myfigure(par,name);
    plot(Nms,time_vfi(:,icase)./time_egm(:,icase),'.--','Color','red','linewidth',2.5,'MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','black')
    xlim([100 600])
    ax.XTick = 100:50:600;
    ylim([0 35])
    ax.YTick = 0:5:35;        
    grid on
    xlabel('nodes, $\#_m$','Interpreter','latex','FontSize', par.fontsize_small);
    ylabel('time relative to G$^2$EGM','Interpreter','latex','FontSize', par.fontsize_small);
    ax.FontSize = par.fontsize_small;
    set(gca,'FontSize',par.fontsize_small)
    fun.printfig(figs.name,'off',casenow);    
        
end

%% 8. Euler errors

for icase = 1:numel(cases)
    
    leg = cell(2,1);
    leg{1} = 'VFI';
    leg{2} = 'G$^2$EGM';

    name = [cases{icase} '_euler_accuracy'];      
    fig = figure('Name',name);
    fig.Color = [1 1 1];
    fig.Visible = 'off';
    ax = axes;    
    ax.FontSize = par.fontsize_small;
    ax.TickLabelInterpreter = 'latex';
    box(ax,'on');  
    plot(Nms,-euler_vfi(:,icase),'.--red','linewidth',2.5,'MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','black')
    set(gca,'FontSize',par.fontsize_small)
    hold on
    plot(Nms,-euler_egm(:,icase),'.-black','linewidth',2.5,'MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','black')
    hold off
    fun.mylegend(par,leg,'NorthEast');
    xlim([100 600])
    ax.XTick = 100:50:600;    
    ylim([-6.5 -2.5])
    ax.YTick = -6.5:0.5:-2.5;    
    grid on
    xlabel('nodes, $\#_m$','Interpreter','latex','FontSize', par.fontsize_small);
    ylabel('Accuracy, $\log_{10}$ rel. Euler error','Interpreter','latex','FontSize', par.fontsize_small);
    ax.FontSize = par.fontsize_small;
    fun.printfig(fig,'off',casenow);
        
end