%% setup
clc;
clear;
close all;

casenow = 'accuracy_3d';
figfolder = ['figures\' casenow];
if isdir(figfolder) == 0
    mkdir(figfolder)
end
figfolder = ['figures\' casenow '\pdf\'];
if isdir(figfolder) == 0
    mkdir(figfolder)
end

LOAD_VFI = 0;
LOAD_EGM = 0;
MAX_THREADS = 1;
T = 20;

cases   = {'nonsmooth_3d','smooth_3d'};
Nms     = [75, 100, 125, 150];
simN    = 100; % number of indivudals simulated to calculate Euler errors.

%% 2. VFI loop

time_vfi   = NaN(numel(Nms),numel(cases));
euler_vfi  = NaN(numel(Nms),numel(cases));

for icase = 1:numel(cases)

    if strcmp('nonsmooth_3d',cases{icase}) ~= 1
        continue;
    end    
           
for i = 1:numel(Nms)
    
    % a. settings
    par = SetupPar_3d();      
    par.Nm          = Nms(i);
    par.max_threads = MAX_THREADS;
    par.T           = T;
    
    name = cases{icase};
    if strcmp(name,'nonsmooth_3d') == 1
        % non-smooth
    else
        % smooth
        par.var_eta = 0.1;
        par.sigma   = 0.1;
    end
    
    fprintf(['\n\n' name ': VFI, Nm = %d\n'],par.Nm);
    
    % b. solve or load   
    if LOAD_VFI == 0;
        t1 = tic;
        sol = solve.vfi(par);
        time_vfi(i,icase) = toc(t1);
        par.time = time_vfi(i,icase);
        save(['data\vfi_Nm' num2str(par.Nm) '_' name '.mat'],'sol','par','-v7.3');       
    else        
        load(['data\vfi_Nm' num2str(par.Nm) '_' name '.mat']);
        time_vfi(i,icase) = par.time;        
    end

    % c. Euler errors 
    sim                      = fun.simulate_euler(sol,simN,par.T,par);
    euler_vfi(i,icase) = nanmean(-log10( abs(sim.euler_work(:)./sim.c(:)) + 1.0e-16));

end
end

%% 3. EGM loop
                               
time_egm   = NaN(numel(Nms),numel(cases));
euler_egm  = NaN(numel(Nms),numel(cases));

for icase = 1:numel(cases)
    
    if strcmp('nonsmooth_3d',cases{icase}) ~= 1
        continue;
    end
     
for i = 1:numel(Nms)
    
    % a. setup
    par             = SetupPar_3d();    
    par.Nm          = Nms(i);
    par.max_threads = MAX_THREADS;
    par.T           = T;

    name = cases{icase};
    if strcmp(name,'nonsmooth_3d') == 1
        % non-smooth
    else
        % smooth
        par.var_eta = 0.1;
        par.sigma   = 0.1;
    end
    fprintf(['\n\n' name ': EGM, Nm = %d\n'],par.Nm);
    
    % b. sols
    if LOAD_EGM == 0        
        t1                   = tic;
        [sol, par]           = solve.egm(par);
        time_egm(i,icase)    = toc(t1);
        par.time = time_egm(i,icase);
        save(['data\egm_Nm' num2str(par.Nm) '_' name '.mat'],'sol','par','-v7.3');
    else    
        load(['data\egm_Nm' num2str(par.Nm) '_' name '.mat'])
        time_egm(i,icase) = par.time;
    end
    
    % c. Euler errors 
    sim                      = fun.simulate_euler(sol,simN,par.T,par);
    euler_egm(i,icase) = nanmean(-log10( abs(sim.euler_work(:)./sim.c(:)) + 1.0e-16));

end
end

%% 4. Figures

for icase = 1:1;
    
    leg = cell(2,1);
    leg{1} = 'VFI';
    leg{2} = 'G$^2$EGM'; 
    
    % a. time
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
    ax.FontSize = par.fontsize_small;
    grid on
    xlim([70 155])
    ax.XTick = 75:25:150;  
    ylim([0 500])
    ax.YTick = 0:50:500;
    xlabel('nodes, $\#_m$','Interpreter','latex','FontSize', par.fontsize_small);
    ylabel('minutes','Interpreter','latex','FontSize', par.fontsize_small);
    ax.FontSize = par.fontsize_small;
    set(gca,'FontSize',par.fontsize_small)
    fun.printfig(fig,'off',casenow);
        
    % b. time_rel
    name = [cases{icase} '_time_rel'];
    [figs.name, ax] = fun.myfigure(par,name);
    plot(Nms,time_vfi(:,icase)./time_egm(:,icase),'.-','Color','black','linewidth',2.5,'MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','black')
    hold on;
    grid on
    ylim([0 50])
    xlim([70 155])
    ax.XTick = 75:25:150;
    ax.YTick = 0:10:100;
    xlabel('nodes, $\#_m$','Interpreter','latex','FontSize', par.fontsize_small);
    ylabel('time relative to G$^2$EGM','Interpreter','latex','FontSize', par.fontsize_small);
    ax.FontSize = par.fontsize_small;
    set(gca,'FontSize',par.fontsize_small)
    fun.printfig(figs.name,'off',casenow);    
    
end

%% 8. Euler errors

for icase = 1:1
    
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
    grid on
    xlim([70 155])
    ax.XTick = 75:25:150;
    ylim([-4.5 -3.0])
    ax.YTick = -4.5:0.25:-3.0;
    xlabel('nodes, $\#_m$','Interpreter','latex','FontSize', par.fontsize_small);
    ylabel('Accuracy, $\log_{10}$ rel. Euler error','Interpreter','latex','FontSize', par.fontsize_small);
    ax.FontSize = par.fontsize_small;
    fun.printfig(fig,'off',casenow);
        
end