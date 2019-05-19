function [par] = SetupPar()
    
    par.ndim = 2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. demograhpics and preferences %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    par.T  = 20; % life-span
    
    % preferences
    par.beta    = 0.98;  % discount factor
    par.rho     = 2.00;  % risk aversion
    par.alpha   = 0.25;  % disutility of labor
    par.sigma   = 0.00;  % variance of taste shock 
    
    % 3d model
    par.gamma     = [];  % labor supply elasticity
    par.varphi    = [];  % scaling of disutility of labor supply
    par.rk        = [];  % return rate on human capital
    par.delta     = [];  % depreciation of human capital
    par.rk_retire = [];  % percentage one-time lump sump income at retirement
    par.lmax      = 0;   % maximum l, only active if > 0 
    
    
    %%%%%%%%%%%%%
    % 2. income %
    %%%%%%%%%%%%%
    
    par.yret    = 0.50; % retirement income
    par.var_eta = 0.00; % variance of income shock when working
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. borrowing and saving %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    par.Ra  = 1.02; % gross return on liquid assets
    par.Rb  = 1.04; % gross return on illiquid pension assets
    par.chi = 0.10; % scaling factor in pension deposit function
    
    
    %%%%%%%%%%%%%%%%%%
    % 4. state space %
    %%%%%%%%%%%%%%%%%%
    
    % solution algorithm settings    
    par.max_threads    = 1;     % maximum number of threads in mex-files

        % VFI
        par.N_guess_vfi    = 75; % guesses in brute force VFI
        par.do_NLopt       = 1;  % use NLopt in VFI
        
        % EGM
        par.save_segments  = 0;     % saves additional data
        par.egm_extrap_add = 2;     % indexes added to bounding box
        par.egm_extrap_w   = -0.25; % minimum weight for extrapolation

        par.delta_con = 0.001; % change of choice in tests
        par.do_dcon_x = 1;     % dcon check of deviation from constraint
        par.do_acon_x = 0;     % acon check of deviation from constraint
        par.do_con_x  = [];    % con check of deviation from constraints
        
    % general
    par.phi_m = 1.1;  % non-linearity of m and a grids
    par.phi_n = 1.25; % non-linearity of n and b grids
    par.phi_k = [];   % non-linearity of k and q grids
    
    % retirement

        % pre-decision states 
        par.Nm_ret    = 500;
        par.m_max_ret = 50;

        % post-decision states        
        par.Na_ret    = par.Nm_ret-floor(par.Nm_ret*0.20);        
        par.a_max_ret = 25;
        
    % working
    
        par.Nm     = 600;
        par.m_max  = 10;
        par.n_add  = 2;
        par.k_max  = [];
        
    % interpolant (w, wa, wb)
    
        par.pd_fac = 2;
        par.a_add  = -2;
        par.b_add  = 2;
        par.q_add  = [];
                
    % egm

        par.acon_fac = 0.25;
        par.con_fac  = 0.50;
        
    % shocks

        par.Neta = 4;           

        
    %%%%%%%%%%%%%%%%%
    % 5. simulation %
    %%%%%%%%%%%%%%%%%

    par.simN     = 10000;
    par.simT     = par.T;
    par.sim_seed = 9210;

    
    %%%%%%%%%%%%%%%%
    % 6. graphical %
    %%%%%%%%%%%%%%%%

    par.fontsize_small = 18;
    par.fontsize_big   = 22;

    par.color      = cell(5);
    par.color{1}   = [3.0/255.0 103.0/255.0 166.0/255.0];
    par.color{2}   = [242.0/255.0 62.0/255.0 46.0/255.0];
    par.color{3}   = [3.0/255.0 166.0/255.0 166.0/255.0];
    par.color{4}   = [242.0/255.0 131.0/255.0 68.0/255.0];
    par.color{5}   = [242.0/255.0 100.0/255.0 48.0/255.0];
    par.color_grey = [65.0/255.0 68.0/255.0 81.0/255.0];

    par.fig_max_m  = 5;
    par.fig_max_n  = 5;
    par.fig_min_k  = [];
    par.fig_max_k  = [];
    par.print      = 0;
    
end