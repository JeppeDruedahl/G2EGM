classdef fun   
% contains: functions used when solving, plotting and simulating the model.
methods (Static)

% 1. figure functions
function [fig, ax] = myfigure(par,name,varargin)
    
    if nargin == 2
        xlabelnow = 'm';
        ylabelnow = 'n';
    else
        xlabelnow = varargin{1};
        ylabelnow = varargin{2};
    end
    
    fig = figure('Name',name);
    fig.Color = [1 1 1];
    fig.Visible = 'off';
    
    ax = axes;    
    ax.FontSize = par.fontsize_small;
    ax.XTick = 0:1:par.fig_max_m;
    ax.YTick = 0:1:par.fig_max_n;
    ax.TickLabelInterpreter = 'latex';
    hold(ax,'on');
    box(ax,'on');
    
    xlim([0 par.fig_max_m])
    ylim([0 par.fig_max_n])
    xlabel(xlabelnow,'Interpreter','latex','FontSize', par.fontsize_small);
    ylabel(ylabelnow,'Interpreter','latex','FontSize', par.fontsize_small);
    
end
function [] = mylegend(par,leg,placement)
    
    h = legend(leg,'Location',placement,'FontSize',par.fontsize_small-2);
    h.Box = 'on';
    h.Interpreter = 'latex';
    
end
function [] = printfig(figin,visible,casenow)
    
    fig = figure(figin);
    fig.Visible = visible;
    fig.PaperUnits = 'centimeters';   
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0 0 16 12];
    fig.PaperSize = [16 12];
    
    if strcmp(casenow,'') == 0
        
         filename = ['figures\' casenow '\pdf\' get(fig,'name') ''];    
         print('-dpdf',['' filename '.pdf']);

        filename = ['figures\' casenow '\' get(fig,'name') ''];    
        print(fig,'-dpng',[filename '.png'],'-opengl','-r300');
        
    else
     
        filename = ['figures\' get(fig,'name') ''];    
        print('-dpdf',['' filename '.pdf']);

        filename = ['figures\' get(fig,'name') ''];    
        print(fig,'-dpng',[filename '.png'],'-opengl','-r300');
    end
    
end

% 2. generic functions
function y = vec(x)
    y = x(:);
end
function x = nonlinspace(lo,hi,n,phi)
    % recursively constructs an unequally spaced grid.
    % phi > 1 -> more mass at the lower end of the grid.
    % lo can be a vector (x then becomes a matrix).
    
    x      = NaN(n,length(lo));
    x(1,:) = lo;
    for i = 2:n
        x(i,:) = x(i-1,:) + (hi-x(i-1,:))./((n-i+1)^phi);
    end
    
end
function [x,w] = GaussHermiteNodes(n)
    % creates Gauss-Hermite nodes x and associated weights w.
    
    if n == 1
        xw = [0 sqrt(pi)];
    elseif n == 2
        xw = ...
            [-7.071067811865476e-1   8.86226925452758e-1,
            7.071067811865476e-1    8.86226925452758e-1];
    elseif n == 4
        xw = ...
            [-1.650680123885785e0    8.13128354472452e-2,
            -5.246476232752904e-1   8.04914090005513e-1,
            5.246476232752904e-1    8.04914090005513e-1,
            1.650680123885785e0     8.13128354472452e-2];
    elseif n == 8
        xw = ...
            [-2.930637420257244e0    1.996040722113676e-4,
            -1.981656756695843e0    1.707798300741347e-2,
            -1.15719371244678e0     2.078023258148919e-1,
            -3.811869902073221e-1   6.611470125582414e-1,
            3.811869902073221e-1    6.611470125582414e-1,
            1.15719371244678e0      2.078023258148919e-1,
            1.981656756695843e0     1.707798300741347e-2,
            2.930637420257244e0     1.996040722113676e-4];
    elseif n == 10
        xw = ...
            [-3.436159118 0.7640432855e-5,
            -2.532731674 0.1343645746e-2,
            -1.756683649 0.3387439445e-1,
            -1.036610829 0.2401386110,
            -0.3429013272 0.6108626337,
            0.3429013272 0.6108626337,
            1.036610829 0.2401386110,
            1.756683649 0.3387439445e-1,
            2.532731674 0.1343645746e-2,
            3.436159118 0.7640432855e-5];
    elseif n == 16
        xw = ...
            [-4.688738939305818e0    2.654807474011183e-10,
            -3.869447904860123e0    2.320980844865211e-7,
            -3.176999161979956e0    2.711860092537881e-5,
            -2.546202157847481e0    9.32284008624181e-4,
            -1.951787990916254e0    1.288031153550997e-2,
            -1.380258539198881e0    8.38100413989858e-2,
            -8.22951449144656e-1    2.806474585285337e-1,
            -2.734810461381524e-1   5.079294790166138e-1,
            2.734810461381524e-1    5.079294790166138e-1,
            8.22951449144656e-1     2.806474585285337e-1,
            1.380258539198881e0     8.38100413989858e-2,
            1.951787990916254e0     1.288031153550997e-2,
            2.546202157847481e0     9.32284008624181e-4,
            3.176999161979956e0     2.711860092537881e-5,
            3.869447904860123e0     2.320980844865211e-7,
            4.688738939305818e0     2.654807474011183e-10];
    elseif n == 32
        xw = ...
            [-7.125813909830728e0    7.310676427384165e-23,
            -6.409498149269661e0    9.23173653651829e-19,
            -5.812225949515914e0    1.197344017092849e-15,
            -5.275550986515881e0    4.215010211326448e-13,
            -4.777164503502596e0    5.933291463396639e-11,
            -4.305547953351199e0    4.098832164770897e-9,
            -3.853755485471445e0    1.574167792545594e-7,
            -3.417167492818571e0    3.650585129562376e-6,
            -2.992490825002374e0    5.416584061819983e-5,
            -2.577249537732317e0    5.36268365527972e-4,
            -2.169499183606112e0    3.654890326654428e-3,
            -1.767654109463201e0    1.755342883157343e-2,
            -1.370376410952872e0    6.045813095591262e-2,
            -9.76500463589683e-1    1.512697340766425e-1,
            -5.849787654359325e-1   2.774581423025299e-1,
            -1.948407415693993e-1   3.752383525928024e-1,
            1.948407415693993e-1    3.752383525928024e-1,
            5.849787654359325e-1    2.774581423025299e-1,
            9.76500463589683e-1     1.512697340766425e-1,
            1.370376410952872e0     6.045813095591262e-2,
            1.767654109463201e0     1.755342883157343e-2,
            2.169499183606112e0     3.654890326654428e-3,
            2.577249537732317e0     5.36268365527972e-4,
            2.992490825002374e0     5.416584061819983e-5,
            3.417167492818571e0     3.650585129562376e-6,
            3.853755485471445e0     1.574167792545594e-7,
            4.305547953351199e0     4.098832164770897e-9,
            4.777164503502596e0     5.933291463396639e-11,
            5.275550986515881e0     4.215010211326448e-13,
            5.812225949515914e0     1.197344017092849e-15,
            6.409498149269661e0     9.23173653651829e-19,
            7.125813909830728e0     7.310676427384165e-23];
    else
        error('unknown number of GaussHermite nodes');
    end
    x = xw(:,1);
    w = xw(:,2);
    
end
function [LogSum, Prob] = logsum(v1,v2,par)
    % calculates the log-sum and choice-probabilities.
    
    % 1. setup
    V           = [v1,v2];
    DIM         = size(v1,1);
    
    % 2. maximum over the discrete choices
    [mxm,id]    = max(V,[],2);
    
    % 3. logsum and probabilities
    if abs(par.sigma) > 1.0e-10
        
        % a. numerically robust log-sum
        LogSum  = mxm + par.sigma*log(sum(exp((V-mxm*ones(1,2))./par.sigma),2));
        
        % b. numerically robust probability
        Prob = exp((V-LogSum*ones(1,2))./par.sigma);
        %Prob = exp(V/par.sigma)./repmat(sum(exp(V/par.sigma),2),1,2);
            
    else % no smoothing -> max-operator
    
        LogSum  = mxm;
        
        Prob    = zeros(DIM,2);
        I       = cumsum(ones(DIM,1)) + (id-1)*DIM; % calculate linear index
        Prob    = zeros(DIM,2);
        Prob(I) = 1;
    
    end

end

% 3. basic function    
function u = u(c,work,l,par)
    % utility function
    if par.ndim == 3
        u = c.^(1-par.rho)/(1-par.rho) - par.alpha*(work==1) - par.varphi*l.^(1+par.gamma)/(1+par.gamma);
    else
        u = c.^(1-par.rho)/(1-par.rho) - par.alpha*(work==1);
    end
end
function marg_u = marg_u(c,work,par)
    % marginal utility function for consumption
    marg_u = c.^(-par.rho);
end
function marg_u_l = marg_u_l(l,work,par)
    % marginal utility function for labor supply
    marg_u_l = par.varphi*l.^(par.gamma);
end
function inv_marg_u = inv_marg_u(u,work,par)
    % inverse utility function
    inv_marg_u = u.^(-1/par.rho);
end
function f_pens = f_pens(p,par)
    % pension deposit function
    f_pens = par.chi*log(1+p);    
end
function trans = trans(v,par)
    % transformation function
    trans = -1.0./v;
end
function trans = trans_inv(v,par)
    % inverse transformation function
    trans = -1.0./v;
end

% 4. preperation
function par = solprep(par)
      
    par.eps   = 1e-6;       % very small number
          
    % 1. retirement
    
        % pre-decision states
        par.grid_m_ret = fun.nonlinspace(par.eps,par.m_max_ret,par.Nm_ret,par.phi_m);
        par.Nmcon_ret  = par.Nm_ret - par.Na_ret;
        
        % post-decision states
        par.grid_a_ret = fun.nonlinspace(0,par.a_max_ret,par.Na_ret,par.phi_m);
        
    % 2. working: state space (m,n,k)    
    par.grid_m  = fun.nonlinspace(par.eps,par.m_max,par.Nm,par.phi_m);
    
    par.Nn      = par.Nm;
    par.n_max   = par.m_max + par.n_add;
    par.grid_n  = fun.nonlinspace(0,par.n_max,par.Nn,par.phi_n);
       
    if par.ndim == 3
        par.Nk = par.Nm;
        par.grid_k  = fun.nonlinspace(1e-2,par.k_max,par.Nk,par.phi_k);
        % nd
        [par.grid_m_nd, par.grid_n_nd, par.grid_k_nd] = ...
            ndgrid(par.grid_m, par.grid_n, par.grid_k);        
    else
        par.Nk     = 1;
        par.grid_k = [];
        % nd
        [par.grid_m_nd, par.grid_n_nd] = ...
            ndgrid(par.grid_m, par.grid_n);
        
    end
                      
    % 2. working: w interpolant (and wa and wb and wq)
    par.Na_pd      = floor(par.pd_fac*par.Nm);
    par.a_max      = par.m_max + par.a_add;
    par.grid_a_pd  = fun.nonlinspace(0,par.a_max,par.Na_pd,par.phi_m);
    
    par.Nb_pd      = floor(par.pd_fac*par.Nn);
    par.b_max      = par.n_max + par.b_add;
    par.grid_b_pd  = fun.nonlinspace(0,par.b_max,par.Nb_pd,par.phi_n);
    
    if par.ndim == 3
        par.Nq_pd      = floor(par.pd_fac*par.Nk);
        par.q_max      = par.k_max + par.q_add;
        par.grid_q_pd  = fun.nonlinspace(0,par.q_max,par.Nq_pd,par.phi_k);
        
        % nd vectors
        [par.grid_a_pd_nd, par.grid_b_pd_nd, par.grid_q_pd_nd] = ...
            ndgrid(par.grid_a_pd, par.grid_b_pd, par.grid_q_pd);
        
    else
        
        par.Nq_pd      = 1;
        par.q_max      = [];
        par.grid_q_pd  = [];
        
        % nd vectors
        [par.grid_a_pd_nd, par.grid_b_pd_nd] = ...
            ndgrid(par.grid_a_pd, par.grid_b_pd);
        par.grid_q_pd_nd = [];
        
    end
    
    % 3. working: egm (seperate grids for each segment)
    
        % a. ucon
        % same as pd
        
        % b. dcon
        % same as pd
        par.d_dcon = zeros(size(par.grid_a_pd_nd));
            
        % c. acon
        par.Nc_acon     = floor(par.Na_pd*par.acon_fac);
        par.Nb_acon     = floor(par.Nb_pd*par.acon_fac);       
        par.grid_b_acon = fun.nonlinspace(0,par.b_max,par.Nb_acon,par.phi_n);
        
        if par.ndim == 3
            par.Nq_acon     = floor(par.Nq_pd*par.acon_fac);
            par.grid_q_acon = fun.nonlinspace(par.eps,par.q_max,par.Nq_acon,par.phi_k);
            [par.b_acon, par.q_acon] = ndgrid(par.grid_b_acon, par.grid_q_acon);
        else
            par.Nq_acon     = [];
            par.grid_q_acon = [];           
            par.b_acon = par.grid_b_acon;
        end
        par.a_acon = zeros(size(par.b_acon));
            
        % d. con
        par.Nc_con     = floor(par.Na_pd*par.con_fac); 
        par.Nb_con     = floor(par.Nb_pd*par.con_fac); 
        par.grid_c_con = fun.nonlinspace(par.eps,par.m_max,par.Nc_con,par.phi_m);
        par.grid_b_con = fun.nonlinspace(0,par.b_max,par.Nb_con,par.phi_n);
        
        if par.ndim == 3
            par.Nq_con     = floor(par.Nq_pd*par.con_fac);     
            par.grid_q_con = fun.nonlinspace(par.eps,par.q_max,par.Nq_con,par.phi_k);
            [par.c_con, par.b_con, par.q_con] = ...
                ndgrid(par.grid_c_con, par.grid_b_con, par.grid_q_con);
        else
            par.Nq_con     = [];     
            par.grid_q_con = [];
            [par.c_con, par.b_con] = ...
                ndgrid(par.grid_c_con, par.grid_b_con);          
        end
        par.a_con = zeros(size(par.c_con));
        par.d_con = zeros(size(par.c_con));
        
    % 4. shocks
    if par.Neta == 1 || par.var_eta == 0
        par.eta   = 1;
        par.w_eta = 1;
        par.Neta  = 1;
    else
        [n, w]    = fun.GaussHermiteNodes(par.Neta);
        nodes     = n*sqrt(2.0);
        par.eta   = exp(sqrt(par.var_eta)*nodes - .5*par.var_eta);
        par.w_eta = w*pi^(-1/2);
    end
       
end
function [sol, interp] = last_period(sol,work,par)
    % allocate memory and solve last period of working.
    
    t1 = tic;
    
    % a. allocate memory
    if par.ndim == 3
        vars = {'c','d','l','v'};
    else
        vars = {'c','d','v'};        
    end
    
        for i = 1:numel(vars);
        for t = 1:par.T
            sol(2,t).(vars{i}) = NaN(par.Nm,par.Nn,par.Nk);
        end
        end
        
    if par.do_derivatives == 1
        vm_next_work = NaN(par.Nm,par.Nn,par.Nk);
        vn_next_work = NaN(par.Nm,par.Nn,par.Nk);
        if par.ndim == 3
            vk_next_work = NaN(par.Nm,par.Nn,par.Nk);
        end
    end
    
    % b. solve last period
    t = par.T;
    for i_n = 1:par.Nn        
    for i_k = 1:par.Nk  
        
        % i. states
        m  = par.grid_m;
        n  = par.grid_n(i_n);
                
        % ii. no labor supply, consume everything and save nothing
        sol(2,t).d(:,i_n,i_k) = zeros(par.Nm,1);        
        if par.ndim == 3
            k                      = par.grid_k(i_k);
            sol(2,t).l(:,i_n,i_k)  = zeros(par.Nm,1);
            sol(2,t).c(:,i_n,i_k)  = m + n + par.rk_retire*k;
        else
            sol(2,t).c(:,i_n,i_k)  = m + n;
        end
               
        % iii. value function
        if par.ndim == 3
            sol(2,t).v(:,i_n,i_k) = fun.u(sol(2,t).c(:,i_n,i_k),work,sol(2,t).l(:,i_n,i_k),par);
        else
            sol(2,t).v(:,i_n,i_k) = fun.u(sol(2,t).c(:,i_n,i_k),work,0,par);
        end
        sol(2,t).v(:,i_n,i_k) = fun.trans(sol(2,t).v(:,i_n,i_k),par);
        
        % iv. value function derivatives
        if par.do_derivatives == 1
            vm_next_work(:,i_n,i_k) = fun.trans(fun.marg_u(sol(2,t).c(:,i_n,i_k),work,par),par);
            vn_next_work(:,i_n,i_k) = vm_next_work(:,i_n,i_k);
            if par.ndim == 3
                vk_next_work(:,i_n,i_k) = fun.trans(fun.marg_u(sol(2,t).c(:,i_n,i_k),work,par)*par.rk_retire,par);                        
            end 
        end
    end    
    end
    
    % c. interpolant for value function derivatives
    if par.do_derivatives == 1
       if par.ndim == 3
            interp.vm_next_work = griddedInterpolant({par.grid_m,par.grid_n,par.grid_k},vm_next_work);
            interp.vn_next_work = griddedInterpolant({par.grid_m,par.grid_n,par.grid_k},vn_next_work);      
            interp.vk_next_work = griddedInterpolant({par.grid_m,par.grid_n,par.grid_k},vk_next_work);
       else
            interp.vm_next_work = griddedInterpolant({par.grid_m,par.grid_n},vm_next_work);
            interp.vn_next_work = griddedInterpolant({par.grid_m,par.grid_n},vn_next_work);           
       end
    else
        interp.vm_next_work.Values = [];
        interp.vn_next_work.Values = [];
        interp.vk_next_work.Values = [];
    end
    
    if par.print == 1
        fprintf(' last period solved in %g secs\n',round(toc(t1)*10)/10);
    end
    
end
            
% 5. discrete choice
function [sol] = discrete_choice(par,sol,t)
    
    interp.v_retire  = griddedInterpolant(sol(1,t).m,sol(1,t).v);
    
    if par.ndim == 3
        m_retire  = par.grid_m_nd+par.grid_n_nd+par.rk_retire*par.grid_k_nd;
    else
        m_retire   = par.grid_m_nd+par.grid_n_nd;        
    end
    v_retire  = interp.v_retire(m_retire(:));
    v_retire  = reshape(v_retire,par.Nm,par.Nn,par.Nk);

    if par.sigma < 10e-6
        sol(2,t).z = zeros(par.Nm,par.Nn,par.Nk);
        I = sol(2,t).v > v_retire;
        sol(2,t).z(I) = 1;
    else
        [N1, N2] = size(v_retire);
        [~, prob] = fun.logsum(sol(2,t).v(:),v_retire(:),par);
        sol(2,t).z = reshape(prob(:,1),N1,N2);
    end            

end

% 6. interpolation objects
function interp = create_interp(t,par,sol,interp)
    % creates all the interpolation objects.
    
    % 1. value function
    if t == par.T-1
        if par.ndim ==3
            gridvectors_work = {par.grid_m,par.grid_n,par.grid_k};
        else
            gridvectors_work = {par.grid_m,par.grid_n};
        end
        interp.v_next_work = griddedInterpolant(gridvectors_work,sol(2,t+1).v);
    else
        interp.v_next_work.Values = sol(2,t+1).v;
    end
    interp.v_next_retire  = griddedInterpolant(sol(1,t+1).m,sol(1,t+1).v);
    
    % 2. value function derivatives for retired households
    if par.do_derivatives == 1
       
        interp.vn_next_retire        = interp.v_next_retire; 
        interp.vn_next_retire.Values = sol(1,t+1).vn;

        interp.vm_next_retire        = interp.v_next_retire; 
        interp.vm_next_retire.Values = sol(1,t+1).vm;        
        
        if par.ndim == 3
            interp.vk_next_retire        = interp.v_next_retire; 
            interp.vk_next_retire.Values = sol(1,t+1).vk; 
        end
        
    end
    
    % 3. post-decision value function (w) and derivates (wa and wb and wq)
    if par.ndim == 3
        if par.Neta ~= 4
            [interp.w_values, wa, wb, wq] = mex_E(par,sol(2,t+1).v,interp.vm_next_work.Values,interp.vn_next_work.Values,interp.vk_next_work.Values,sol(1,t+1));
        else
            [interp.w_values, wa, wb, wq] = mex_E_vec(par,sol(2,t+1).v,interp.vm_next_work.Values,interp.vn_next_work.Values,interp.vk_next_work.Values,sol(1,t+1));
        end      
    else
        if par.Neta ~= 4
            [interp.w_values, wa, wb] = mex_E(par,sol(2,t+1).v,interp.vm_next_work.Values,interp.vn_next_work.Values,sol(1,t+1));
        else
            [interp.w_values, wa, wb] = mex_E_vec(par,sol(2,t+1).v,interp.vm_next_work.Values,interp.vn_next_work.Values,sol(1,t+1));
        end
    end

        % main interpolant
        if t == par.T-1
            if par.ndim == 3
                interp.w = griddedInterpolant({par.grid_a_pd,par.grid_b_pd,par.grid_q_pd},...
                    interp.w_values);
            else
                interp.w = griddedInterpolant({par.grid_a_pd,par.grid_b_pd},...
                    interp.w_values);                
            end
        else
            interp.w.Values = interp.w_values;
        end
        
        % derivatives
        if par.do_derivatives == 1
            
            if t == par.T-1 
                interp.wa = interp.w;
                interp.wb = interp.w;
                if par.ndim == 3
                    interp.wq = interp.w;
                end
            end
            
            interp.wa.Values = wa;            
            interp.wb.Values = wb;            
            if par.ndim == 3
                interp.wq.Values = wq;
            end
        end        
end

% 7. upperenvelope
function [sol,par] = upperenvelopetocommon(par,seg,interp)
    
    if par.ndim == 3
        [seg.Na, seg.Nb, seg.Nq] = size(seg.v); 
    else
        [seg.Na, seg.Nb] = size(seg.v); 
        seg.Nq = 1;
        seg.k  = [];
        seg.q  = [];
    end
        
    % 1. indicator for valid (and interesting choice or not)
    valid = imag(seg.c) == 0 & imag(seg.d) == 0 & isnan(seg.v) == 0;
    valid = valid == 1 & seg.c >= -0.50 & seg.d >= -0.50;
    valid = valid == 1 & seg.m > -0.1 & seg.n > -0.1 ;
    valid = valid == 1 & seg.m < par.m_max + 1 & seg.n < par.n_max + 1; 
    
    if par.ndim == 3
        valid = valid == 1 & imag(seg.l) == 0 & seg.l >= -0.50;
        valid = valid == 1 & seg.k > -2 & seg.k < par.k_max + 2;
        if par.lmax > 0           
            valid = valid == 1 & seg.l <= par.lmax;
        end
    end
            
    if isfield(seg,'valt') == 1
        for j = 1:numel(seg.valt)
            valid = valid == 1 & seg.v > seg.valt{j};
        end
    end
    sol.percent_cleaned = sum(valid(:)==0)/numel(valid(:))*100;  
   
    % 2. upper envelope
    if sum(valid(:)) < 100
            
        if par.ndim == 3
            vars = {'c','d','l','v'};
            sizenow = [par.Nm,par.Nn,par.Nk];
        else
            vars = {'c','d','v'};
            sizenow = [par.Nm,par.Nn];
        end
        for i = 1:numel(vars)
            sol.(vars{i})  = NaN(sizenow);
        end
        sol.time_upper = 0;
        
    else   
        
        t1 = tic;
        seg.valid = logical(valid);  
        if par.ndim == 3           
            [sol.c, sol.d, sol.l, sol.v, sol.holes] = ...
                mexUpperEnvelopeToCommon(par,interp,seg);
        else
            [sol.c, sol.d, sol.v, sol.holes] = ...
                mexUpperEnvelopeToCommon(par,interp,seg);            
        end
        sol.time_upper = toc(t1);
                
    end  
    
    % print
    if par.print == 1
        names = {'   ucon','   dcon','   acon','    con','   lcon','  ldcon','  lacon','fullcon'};
        fprintf(['   ' names{seg.num} ': %4.1f secs (upper: %3.1f) (nodes %dk, %5.3f%% cleaned)\n'],...
                    toc(seg.time),...
                    sol.time_upper,...
                    floor(numel(seg.v)/1000),...
                    sol.percent_cleaned); 
    end
    
end

% 8. simulate
function sim = simulate_euler(sol,N,T,par)
     
    % a. setup
    Nini = N;
    T    = T-1;
    N    = N^par.ndim; 
    par  = fun.solprep(par);
     
    % b. grid
    min_m = 0.50;
    min_n = 0.01;    
    min_k = 1.00;
     
    m_max = 5.00;
    n_max = 5.00;
    k_max = 30.0;
     
    if par.ndim == 3
        m_grid = linspace(min_m,m_max,Nini);
        n_grid = linspace(min_n,n_max,Nini);
        k_grid = linspace(min_k,k_max,Nini);    
        [m_grid,n_grid,k_grid] = ndgrid(m_grid,n_grid,k_grid);
        m = m_grid(:);
        n = n_grid(:); 
        k = k_grid(:); 
    else
        m_grid = linspace(min_m,m_max,Nini);
        n_grid = linspace(min_n,n_max,Nini);
        [m_grid,n_grid] = ndgrid(m_grid,n_grid);
        m = m_grid(:);
        n = n_grid(:);         
    end
    z_lag = ones(size(m)); 
     
    % c. shocks
    rng(9210);
    uniform = rand(N,T); % uniform draws in the case with taste shocks
     
    % d. loop over time
    sim.euler      = NaN(N,T-1);
    sim.euler_work = NaN(N,T-1);
     
    for t = 1:T-1
         
        % i. interpolants: current period
        if par.ndim == 3
            gridvectors_work = {par.grid_m,par.grid_n,par.grid_k};
        else
            gridvectors_work = {par.grid_m,par.grid_n};
        end         
        interp_v_work = griddedInterpolant(gridvectors_work,sol(2,t).v);
        interp_c_work = griddedInterpolant(gridvectors_work,sol(2,t).c);
        interp_d_work = griddedInterpolant(gridvectors_work,sol(2,t).d); 
        if par.ndim == 3
            interp_l_work = griddedInterpolant(gridvectors_work,sol(2,t).l);
        end
        interp_v_retire  = griddedInterpolant(sol(1,t).m,sol(1,t).v);
        interp_c_retire  = griddedInterpolant(sol(1,t).m,sol(1,t).c);
         
        % ii. interpolants: next-period
        interp_v_work_next = griddedInterpolant(gridvectors_work,sol(2,t+1).v);
        interp_c_work_next = griddedInterpolant(gridvectors_work,sol(2,t+1).c);
        if par.ndim == 3
            interp_l_work_next = griddedInterpolant(gridvectors_work,sol(2,t+1).l);
        end
        interp_v_retire_next  = griddedInterpolant(sol(1,t+1).m,sol(1,t+1).v);
        interp_c_retire_next  = griddedInterpolant(sol(1,t+1).m,sol(1,t+1).c);
         
        % iii. optimal retirement choice       
        if par.ndim == 3
            v_work       = fun.trans_inv(interp_v_work(m,n,k),par);
            m_retire     = m+n+par.rk_retire*k;
        else
            v_work       = fun.trans_inv(interp_v_work(m,n),par);            
            m_retire     = m+n;            
        end
        v_retire         = fun.trans_inv(interp_v_retire(m_retire),par);
        [~, Prob]        = fun.logsum(v_retire,v_work,par);
        Prob(z_lag==0,1) = 1;
        Prob(z_lag==0,2) = 0;
         
        z = sum(repmat(uniform(:,t),1,2) > cumsum(Prob,2),2);
        sim.z(:,t) = z;
         
        % iv. optimal continuous choices 
        if par.ndim == 3
            sim.l(:,t) = (z==1).*interp_l_work(m,n,k);
            cmax_work  = m + par.rk*k.*sim.l(:,t);         
            sim.c(:,t) = (z==1).*min(interp_c_work(m,n,k),cmax_work) + ...
                         (z==0).*min(interp_c_retire(m_retire),m_retire);
            sim.d(:,t) = (z==1).*interp_d_work(m,n,k);
        else
            sim.c(:,t) = (z==1).*min(interp_c_work(m,n),m) + ...
                         (z==0).*min(interp_c_retire(m_retire),m_retire);
            sim.d(:,t) = (z==1).*interp_d_work(m,n);            
        end
         
        % v. states after retirement choice
        sim.n(:,t)  = n.*(z==1) + 0.*(z==0);
        sim.m(:,t)  = m.*(z==1) + (m_retire).*(z==0);
        if par.ndim == 3
            sim.k(:,t)  = k.*(z==1) + 0.*(z==0);
        else
        end
         
        % vi. post-decision states
        sim.b(:,t)  = sim.n(:,t) + sim.d(:,t) + fun.f_pens(sim.d(:,t),par); 
        if par.ndim == 3        
            sim.a(:,t)  = sim.m(:,t) + par.rk*sim.k(:,t).*sim.l(:,t) - sim.c(:,t) - sim.d(:,t);
            sim.q(:,t)  = (1-par.delta)*sim.k(:,t) + sim.l(:,t);
        else
             sim.a(:,t)  = sim.m(:,t) - sim.c(:,t) - sim.d(:,t);
        end
         
        % vii. loop over shocks
        E = 0.0;
        for i_eta=1:par.Neta
             
            % o. state variables
            n_next    = par.Rb*sim.b(:,t);            
            if par.ndim == 3
                k_next = par.eta(i_eta)*sim.q(:,t);
                m_next = par.Ra*sim.a(:,t) + par.yret.*(z==0);
            else
                m_next = par.Ra*sim.a(:,t) + par.yret.*(z==0) + par.eta(i_eta).*(z==1);
            end
            % oo. optimal discrete choice
            if par.ndim == 3
                v_work   = fun.trans_inv(interp_v_work_next(m_next,n_next,k_next),par);
            else
                v_work   = fun.trans_inv(interp_v_work_next(m_next,n_next),par);
            end
            v_retire     = fun.trans_inv(interp_v_retire_next(m_next),par);
            [~, Prob]    = fun.logsum(v_retire,v_work,par);
            Prob(z==0,1) = 1;
            Prob(z==0,2) = 0;
 
            z_next = sum(repmat(uniform(:,t+1),1,2) > cumsum(Prob,2),2);
             
            % ooo. optimal continuous choices
            if par.ndim == 3
                m_next_retire = m_next + n_next + par.rk_retire*k_next;
                l_next          = (z_next==1).*interp_l_work_next(m_next,n_next,k_next);
                cmax_work_next  = m_next + par.rk*k_next.*l_next;            
                c_next          = (z_next==1).*min(interp_c_work_next(m_next,n_next,k_next),cmax_work_next) + ...
                                  (z_next==0).*min(interp_c_retire_next(m_next_retire),m_next_retire);
            else
                m_next_retire   = m_next + n_next;
                c_next          = (z_next==1).*min(interp_c_work_next(m_next,n_next),m_next) + ...
                                  (z_next==0).*min(interp_c_retire_next(m_next_retire),m_next_retire);    
            end
            
            % oooo. weighted sum
            E = E + par.w_eta(i_eta)* par.beta*par.Ra*fun.marg_u(c_next,z_next,par);         
         
        end
         
        % Euler error
        I = sim.a(:,t) >= 0.001;
        sim.euler(I,t) = sim.c(I,t) - fun.inv_marg_u(E(I),z,par);
        I = I==1 & z==1;
        sim.euler_work(I,t) = sim.euler(I,t);
         
        mean_log10_euler = nanmean(log10( abs(sim.euler_work(:,t)./sim.c(:,t)) + 1.0e-16));
        fprintf(' t = %d: mean_log10_euler = %5.3f (%d)\n',t,mean_log10_euler,sum(I));
 
    end
end

%% over and out
end
end