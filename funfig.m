classdef funfig
methods (Static)
    
function figs = overview(figs,par,sol,k,ext)
    % create overview figures
    
    % data
    par = fun.solprep(par);
    
    a = par.grid_m_nd - sol(2,par.T-k).c - sol(2,par.T-k).d;
    d = max(sol(2,par.T-k).d , 0);
    m = par.grid_m_nd;
    n = par.grid_n_nd;
    z = sol(2,par.T-k).z;
    
    I = a < 1e-7;
    a(I) = 0;
    
    I = m < par.fig_max_m & n < par.fig_max_n;
    Icon = a == 0 & d == 0 & I == 1;
    Iucon = a > 0 & d > 0 & I == 1;
    Iacon = a == 0 & d > 0 & I == 1;
    Idcon = a > 0 & d == 0 & I == 1;
    
    % sections
    name = ['T' num2str(k,'%02.0f') '_overview_' ext ];
    [figs.overview.(ext), ~] = fun.myfigure(par,name);     
        ileg = 0;
        if sum(fun.vec(Icon)) > 0
            ileg = ileg+1;
            leg{ileg} = 'con';
            scatter(m(Icon),n(Icon),4,par.color{3},'filled');
        end
        if sum(fun.vec(Iacon)) > 0
            ileg = ileg+1;
            leg{ileg} = 'acon';
            scatter(m(Iacon),n(Iacon),4,par.color{1},'filled');
        end
        if sum(fun.vec(Idcon)) > 0
            ileg = ileg+1;
            leg{ileg} = 'dcon';
            scatter(m(Idcon),n(Idcon),4,par.color{2},'filled');
        end
        if sum(fun.vec(Iucon)) > 0
            ileg = ileg+1;
            leg{ileg} = 'ucon';
            scatter(m(Iucon),n(Iucon),4,'black','filled');  
        end
    fun.mylegend(par,leg,'NorthWest');
    
    % discrete choice
    if par.sigma == 0
        
        retire = z(:) == 0;
        work = z(:) == 1;

        name = ['T' num2str(k,'%02.0f') '_discrete_' ext ];
        [figs.discrete.(ext), ~] = fun.myfigure(par,name);        
            scatter(m(retire),n(retire),4,par.color{1},'filled');
            scatter(m(work),n(work),5,par.color{2},'filled');
        fun.mylegend(par,{'retire','work'},'best');
        
    else
        
        name = ['T' num2str(k,'%02.0f') '_discrete_' ext ];
        [figs.discrete.(ext), ~] = fun.myfigure(par,name);            
            scatter(m(:),n(:),4,z(:),'filled'); % colopmap off due to limits
            colormap(jet);
            colorbar;
            
    end
    
end

function figs = compare(figs,par,sol_compare,sol_egm,k,varlist)
    % create comparision figures
    
    par = fun.solprep(par);    
    
    sol_egm(2,par.T-k).a = par.grid_m_nd - sol_egm(2,par.T-k).c - sol_egm(2,par.T-k).d;
    sol_egm(2,par.T-k).b = par.grid_n_nd + sol_egm(2,par.T-k).d + fun.f_pens(sol_egm(2,par.T-k).d,par);  

    sol_compare(2,par.T-k).a = par.grid_m_nd - sol_compare(2,par.T-k).c - sol_compare(2,par.T-k).d;
    sol_compare(2,par.T-k).b = par.grid_n_nd + sol_compare(2,par.T-k).d + fun.f_pens(sol_compare(2,par.T-k).d,par);  
    
    % loop over variables
    for i = 1:numel(varlist)
        
        var = varlist{i};
        
        % egm - data
        m_egm = par.grid_m_nd(:);
        n_egm = par.grid_n_nd(:);
        I_egm = m_egm < par.fig_max_m & n_egm < par.fig_max_n;
        egm    = sol_egm(2,par.T-k).(var)(:);
        
        % vfi - data
        m_vfi  = par.grid_m_nd(:);
        n_vfi  = par.grid_n_nd(:);
        I_vfi  = m_vfi < par.fig_max_m & n_vfi < par.fig_max_n;
        vfi    = sol_compare(2,par.T-k).(var)(:);
        
        % 1. vf  - points            
        % name = ['T' num2str(k,'%02.0f') '_' var '_vfi'];    
        % [figs.(var).vfi, ~] = fun.myfigure(par,name);
        %    scatter3(m_vfi(I_vfi),n_vfi(I_vfi),vfi(I_vfi),2,vfi(I_vfi),'filled')
        %    view([-37.5 30]);
        %    grid on;
            
        % 2. egm - points            
        name = ['T' num2str(k,'%02.0f') '_' var '_egm'];    
        [figs.(var).egm, ~] = fun.myfigure(par,name);
            scatter3(m_egm(I_egm),n_egm(I_egm),egm(I_egm),2,egm(I_egm),'filled')
            view([-37.5 30]);
            grid on;
            
    end
    
end

function figs = segments(figs,par,sol_compare,sol,k,stage,segs,varlist)
    % creates segment figures
    
    par = fun.solprep(par);
        
    % loop over segments
    for i = 1:numel(segs)
        

        seg = segs{i};
        if strcmp(stage,'clean') == 1
            stagename = ['a_' stage];            
        elseif strcmp(stage,'upper') == 1
            stagename = ['b_' stage];
        elseif strcmp(stage,'final')
            stagename = ['c_' stage];  
        elseif strcmp(stage,'final_opt')
            stagename = ['d_' stage];
        else
            error('unknown stage')
        end
        if  strcmp(stage,'final') || strcmp(stage,'final_opt')
            solnow = sol(2,par.T-k).(seg);
            solnow.m = par.grid_m_nd;
            solnow.n = par.grid_n_nd;
        else
            solnow = sol(2,par.T-k).(seg).(stage);
        end
        solnow.a = solnow.m - solnow.c - solnow.d;
        solnow.b = solnow.n + solnow.d + fun.f_pens(solnow.d,par);
            
        % mn based:
        I = solnow.m < par.fig_max_m & solnow.n < par.fig_max_n & isnan(solnow.v) == 0 & isinf(-solnow.v) == 0;
        
        % i. (m,n)
        if strcmp(stage,'final_opt')
            name = ['T' num2str(k,'%02.0f') '_grid_'  seg '_' stagename];
            figs.grid.(seg).(stage) = fun.myfigure(par,name);
            scatter(solnow.m(I),solnow.n(I),5,.7*ones(1,3),'filled');
            a = par.grid_m_nd - sol(2,par.T-k).c - sol(2,par.T-k).d;
            d = sol(2,par.T-k).d;
            m = par.grid_m_nd;
            n = par.grid_n_nd;
            z = sol(2,par.T-k).z;
            if strcmp(seg,'con')
            I_seg = a < 1e-8 & d == 0;
            elseif strcmp(seg,'ucon')
            I_seg = a >= 1e-8 & d > 0;
            elseif strcmp(seg,'acon')
            I_seg = a < 1e-8 & d > 0;
            elseif strcmp(seg,'dcon')
            I_seg = a >= 1e-8 & d == 0;
            end
            scatter(m(I_seg),n(I_seg),5,'black','filled');
            grid on;
        else
            name = ['T' num2str(k,'%02.0f') '_grid_'  seg '_' stagename];
            figs.grid.(seg).(stage) = fun.myfigure(par,name);
            scatter(solnow.m(I),solnow.n(I),5,'black','filled');
            grid on;
            
            name = ['T' num2str(k,'%02.0f') '_grid_'  seg '_' stagename];
            figs.gridab.(seg).(stage) = fun.myfigure(par,name);
            scatter(solnow.a(I),solnow.b(I),5,'black','filled');
            grid on;
            xlabel('a','Interpreter','latex','FontSize', par.fontsize_small);
            ylabel('b','Interpreter','latex','FontSize', par.fontsize_small);            
            
        end    
    end
           
end

%% over and out
end
end