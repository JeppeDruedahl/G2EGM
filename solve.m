classdef solve
% contains: functions for solving the model
methods (Static)
    
% 0. solve the model using various methods
function sols = all(sols,methodlist,par)
    for i = 1:numel(methodlist)
        method = methodlist{i};   
        sols.(method) = solve.(method)(par);
    end
end

% 1. retirement: EGM
function sol = ret(par)
    
    t1 = tic;
    
    work_next = 0;
    work      = 0;  
    
    % 1. allocate memory
    sol(2,par.T).v = [];
    if par.ndim == 3
        vars = {'m','v','c','a','vm','vn','vk'};
    else
        vars = {'m','v','c','a','vm','vn'};        
    end
    for i = 1:numel(vars);
    for t = 1:par.T
        sol(1,t).(vars{i}) = NaN(par.Nm_ret,1);
    end
    end
    
    % 2. last period
    t = par.T;
    sol(1,t).m  = par.grid_m_ret;
    sol(1,t).c  = sol(1,t).m;
    sol(1,t).a  = zeros(par.Nm_ret,1);
    sol(1,t).v  = fun.trans(fun.u(sol(1,t).c,work,0,par),par);
    sol(1,t).vm = fun.trans(fun.marg_u(sol(1,t).c,work,par),par);
    sol(1,t).vn = sol(1,t).vm;
    if par.ndim == 3
        sol(1,t).vk = fun.trans(par.rk_retire*fun.marg_u(sol(1,t).c,work,par));
    end
    
    % 3. backwards induction
    for t = par.T-1:-1:1
        
        % a. optimal c choice
        m_next   = par.Ra*par.grid_a_ret + par.yret;
        c_next   = interp1(sol(1,t+1).m,sol(1,t+1).c,m_next,'linear','extrap');
        vm_next  = fun.marg_u(c_next,work_next,par);        
        sol(1,t).c(par.Nmcon_ret+1:end) = fun.inv_marg_u(par.beta*par.Ra*vm_next,work,par);
        
        % b. constraint
        sol(1,t).c(1:par.Nmcon_ret) = fun.nonlinspace(1e-6,sol(1,t).c(par.Nmcon_ret+1)*0.999,par.Nmcon_ret,par.phi_m);

        % c. end-of-period assets and value-of-choice
        sol(1,t).a = [par.grid_a_ret(1)*ones(par.Nmcon_ret,1);par.grid_a_ret];        
        v_next     = fun.trans_inv(interp1(sol(1,t+1).m,sol(1,t+1).v,m_next,'linear','extrap'));                       
        sol(1,t).v = fun.u(sol(1,t).c,work,0,par) + par.beta*[v_next(1)*ones(par.Nmcon_ret,1);v_next];
        sol(1,t).v = fun.trans(sol(1,t).v,par);
                
        % d. endogenous grid
        sol(1,t).m = sol(1,t).a + sol(1,t).c;

        % e. marginal v
        sol(1,t).vm = fun.trans(fun.marg_u(sol(1,t).c,work,par),par);
        sol(1,t).vn = sol(1,t).vm;
        if par.ndim == 3
            sol(1,t).vk = fun.trans(par.rk_retire*fun.marg_u(sol(1,t).c,work,par));
        end;
        
    end % t
    
    if par.print == 1
        fprintf(' retirement solved in %3.1f secs\n',toc(t1));
    end
    
end

% 1. working: value function iteration
function sol = vfi(par)
    
    if par.print == 1
        t0 = tic;        
        fprintf('vfi:\n');       
    end
    work = 1;
    
    % settings
    par = fun.solprep(par);
    par.do_derivatives = 0;

    % 1. solve retirement by EGM
    sol = solve.ret(par); 
    
    % 2. last period
    [sol, interp] = fun.last_period(sol,work,par);

    % 3. backward iteration
    for t = par.T-1:-1:1
               
        if par.print == 1
            t1 = tic;
            fprintf(' t = %d\n',t);
        end
        
        % a. interpolation objects        
        t2 = tic;
        interp = fun.create_interp(t,par,sol,interp);
        if par.print == 1
            fprintf('  w-interpolant: %3.1f secs\n',toc(t2));             
        end
        
        % b. loop over beginning-of-period states
        if par.do_NLopt == 0
            [sol(2,t).c,sol(2,t).d,sol(2,t).v] = ...
                mexVFI(par,interp);  
        elseif par.ndim == 3
            [sol(2,t).c,sol(2,t).d,sol(2,t).l,sol(2,t).v] = ...
                mexVFI_NLopt(par,interp);               
        else
            [sol(2,t).c,sol(2,t).d,sol(2,t).v] = ...
                mexVFI_NLopt(par,interp);                           
        end
        sol(2,t).v = fun.trans(sol(2,t).v,par);
        
        % c. discrete choice
        [sol] = fun.discrete_choice(par,sol,t);
        
        % d. timing
        if par.print == 1
            fprintf('  time: %4.1f secs\n',toc(t1));             
        end
        
    end % t
    
    if par.print == 1
        fprintf(' solution time: %4.2f mins (%4.0f secs)\n',...
            toc(t0)/60,...
            toc(t0)); 
    end
    
end

% 2. working: endogenous grid method
function [sol, par] = egm(par)
    
    t0 = tic;    
    if par.print == 1
        fprintf('egm:\n');     
    end
    work = 1;    
    time_tot_tot_upper = 0;
        
    % settings
    par = fun.solprep(par);
    par.do_derivatives = 1;

    % 1. solve retirement by EGM
    sol = solve.ret(par); 
    
    % 2. allocate memory
    [sol, interp] = fun.last_period(sol,work,par);
           
    % 3. backward iterations
    for t = par.T-1:-1:1

        t1 = tic;        
        if par.print == 1
            fprintf(' t = %d\n',t);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % a. interpolation objects %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        t2 = tic;
        interp = fun.create_interp(t,par,sol,interp);
        if par.print == 1
            fprintf('  w-interpolant: %3.1f secs\n',toc(t2));             
            fprintf('  segments:\n');
        end     
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % b. unconstrained (ucon) %%     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ucon.time = tic;
        ucon.num = 1;
        
        % i. choices
        ucon.c = egmfun.foc_c(interp.wa.Values,work,par);
        ucon.d = egmfun.foc_d(interp.wa.Values,interp.wb.Values,par);        
        if par.ndim == 3
            ucon.l = egmfun.foc_l(par.grid_q_pd_nd,interp.wa.Values,interp.wq.Values,par);
        else
            ucon.l  = [];
        end
        
        % ii. states and value
        [ucon.v, ucon.k, ucon.m, ucon.n] = egmfun.inv_kmn_and_v(...
            par.grid_a_pd_nd,par.grid_b_pd_nd,par.grid_q_pd_nd,...
            ucon.c,ucon.d,ucon.l,...
            interp.w.Values, work, par);
        
        % iii. upperenvelope and intepr to common
        [sol(2,t).ucon, par] = fun.upperenvelopetocommon(par,ucon,interp);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % c. constrained: d = 0 (dcon) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        dcon.time = tic;
        dcon.num  = 2;
        
        % i. decisions                
        dcon.c = ucon.c;
        dcon.d = par.d_dcon;
        dcon.l = ucon.l;
        
        % ii. states and value
        [dcon.v, dcon.k, dcon.m, dcon.n] = egmfun.inv_kmn_and_v(...
            par.grid_a_pd_nd,par.grid_b_pd_nd,par.grid_q_pd_nd,...
            dcon.c,dcon.d,dcon.l,...
            interp.w.Values,work,par);
                           
        % iii. value of deviating a bit from the constraint
        if par.do_dcon_x == 1         
            dcon.valt{1} = egmfun.deviate_d_con(...
                par.grid_a_pd_nd,par.grid_q_pd_nd,dcon,work,par,interp);
        end
        
        % v. upperenvelope and interp to common
        [sol(2,t).dcon, par] = fun.upperenvelopetocommon(par,dcon,interp);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % d. constrained: a = 0 (acon) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        acon.time = tic; 
        acon.num  = 3;
        
        if par.ndim == 3
            
            % i. basic setup
            wb_acon                  = interp.wb(par.a_acon,par.b_acon,par.q_acon);
            [c_min, c_max]           = egmfun.c_min_max_acon(wb_acon,work,par);
            [acon.c, b_acon, q_acon] = mex_acon_grid(c_min,c_max,par);
            a_acon                   = zeros(size(acon.c));

            % ii. interpolants        
            wb_acon = interp.wb(a_acon,b_acon, q_acon);
            wq_acon = interp.wq(a_acon,b_acon,q_acon);
            w_acon  = interp.w(a_acon,b_acon,q_acon);
            
            % iii. choices
            acon.d  = egmfun.foc_d_ccon(acon.c,wb_acon,work,par);             
            acon.l  = egmfun.foc_l_ccon(acon.c,q_acon,wq_acon,work,par);
        
        else
            
            % i. setup
            wb_acon = interp.wb(par.a_acon,par.b_acon);
            c_min = fun.inv_marg_u((par.chi+1)*par.beta*wb_acon,work,par);
            c_max = fun.inv_marg_u(par.beta*wb_acon,work,par)';
            
            % ii. choices
            acon.c  = fun.nonlinspace(c_min,c_max,par.Nc_acon,par.phi_m)';
            wb_acon = wb_acon*ones(1,par.Nc_acon);
            acon.d  = par.chi./(fun.marg_u(acon.c,work,par)./(par.beta*wb_acon)-1)-1;
            acon.l  = [];
            
            % iii. post-decision states and value function
            a_acon = zeros(size(acon.c));
            b_acon = par.b_acon*ones(1,par.Nc_acon);
            w_acon = interp.w(par.a_acon,par.b_acon)*ones(1,par.Nc_acon);            
            q_acon = [];
            
        end
        
        % iv. states and value
        [acon.v, acon.k, acon.m, acon.n] = egmfun.inv_kmn_and_v(...
            a_acon,b_acon,q_acon,...
            acon.c,acon.d,acon.l,...
            w_acon,work,par);
                    
        % v. value of deviating a bit from the constraint
        if par.do_acon_x == 1
            acon.valt{1} = egmfun.deviate_a_con(b_acon,q_acon,acon,work,par,interp);
        end   
        
        % vi. upperenvelope and interp to common
        [sol(2,t).acon, par] = fun.upperenvelopetocommon(par,acon,interp);

            if par.save_segments == 1
                sol(2,t).acon.clean = acon;
            end
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % e. constrained: a = 0, d = 0 (con) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
        con.time = tic;
        con.num  = 4;
        
        if par.ndim == 3
                        
            % i. interpolations       
            wq_con  = interp.wq(par.a_con,par.b_con,par.q_con);
            w_con   = interp.w(par.a_con,par.b_con,par.q_con);

            % ii. decisions
            con.c = par.c_con;
            con.l = egmfun.foc_l_ccon(con.c,par.q_con,wq_con,work,par);            
            con.d = par.d_con;

            % iii. value pre-decision states
            [con.v, con.k, con.m, con.n] = egmfun.inv_kmn_and_v(...
                par.a_con,par.b_con,par.q_con,...
                con.c,con.d,con.l,...
                w_con,work,par);

            % iv. value of deviating a bit from the constraints
            if par.do_con_x == 1
                con.valt{1} = egmfun.deviate_d_con(par.a_con,par.q_con,con,work,par,interp);
                con.valt{2} = egmfun.deviate_a_con(par.b_con,par.q_con,con,work,par,interp);            
            end

            % v. upperenvelope and interp to common
            [sol(2,t).con, par] = fun.upperenvelopetocommon(par,con,interp);
            
        else
                     
            % i. choices
            sol(2,t).con.c = par.grid_m_nd;
            sol(2,t).con.d = zeros(par.Nm,par.Nn);
        
            % ii. post-decision states and value
            sol(2,t).con.v = fun.u(sol(2,t).con.c,work,[],par) + par.beta*...
                interp.w(zeros(par.Nm,par.Nn),par.grid_n_nd);            
            sol(2,t).con.time_upper = 0;            
            if par.print == 1
                fprintf('   con: %4.1f secs\n',toc(con.time));  
            end
        
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % f. constrained: l = lmax (lcon) %     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if par.lmax > 0

            lcon.time = tic;
            lcon.num  = 5;
            
            % i. choices
            lcon.c = ucon.c;
            lcon.d = ucon.d;
            lcon.l = par.lmax*ones(size(par.grid_a_pd_nd));
            
            % ii. states and value
            [lcon.v, lcon.k, lcon.m, lcon.n] = egmfun.inv_kmn_and_v(...
                par.grid_a_pd_nd,par.grid_b_pd_nd,par.grid_q_pd_nd,...
                lcon.c,lcon.d,lcon.l,...
                interp.w.Values, work, par);

            % iii. value of deviating a bit from the constraint
            if par.do_lcon_x == 1
                lcon.valt{1} = egmfun.deviate_l_con(...
                    par.grid_a_pd_nd,par.grid_b_pd_nd,lcon,work,par,interp);
            end
            
            % iv. upperenvelope and intepr to common
            [sol(2,t).lcon, par] = fun.upperenvelopetocommon(par,lcon,interp);
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % g. constrained: d = 0, l = lmax (ldcon) %     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if par.lmax > 0

            ldcon.time = tic;
            ldcon.num  = 6;
            
            % i. choices
            ldcon.c = ucon.c;
            ldcon.d = zeros(size(par.grid_a_pd_nd));
            ldcon.l = par.lmax*ones(size(par.grid_a_pd_nd));
            
            % ii. states and value
            [ldcon.v, ldcon.k, ldcon.m, ldcon.n] = egmfun.inv_kmn_and_v(...
                par.grid_a_pd_nd,par.grid_b_pd_nd,par.grid_q_pd_nd,...
                ldcon.c,ldcon.d,ldcon.l,...
                interp.w.Values, work, par);

            % iii. value of deviating a bit from the constraints
            if par.do_ldcon_x == 1
                ldcon.valt{1} = egmfun.deviate_d_con(...
                    par.grid_a_pd_nd,par.grid_q_pd_nd,ldcon,work,par,interp);
                ldcon.valt{2} = egmfun.deviate_l_con(...
                    par.grid_a_pd_nd,par.grid_b_pd_nd,ldcon,work,par,interp);
            end
            
            % iv. upperenvelope and intepr to common
            [sol(2,t).ldcon, par] = fun.upperenvelopetocommon(par,ldcon,interp);
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % h. constrained: a = 0, l = lmax (lacon) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if par.lmax > 0
        
            lacon.time = tic; 
            lacon.num  = 7;

            % i. choices
            lacon.c = acon.c;
            lacon.d = acon.d;             
            lacon.l = par.lmax*ones(size(lacon.c));

            % ii. states and value
            [lacon.v, lacon.k, lacon.m, lacon.n] = egmfun.inv_kmn_and_v(...
                a_acon,b_acon,q_acon,...
                lacon.c,lacon.d,lacon.l,...
                w_acon,work,par);

            % iii. value of deviating a bit from the constraint
            if par.do_lacon_x == 1
                lacon.valt{1} = egmfun.deviate_a_con(...
                    b_acon,q_acon,lacon,work,par,interp);
                lacon.valt{2} = egmfun.deviate_l_con(...
                    a_acon,b_acon,lacon,work,par,interp);
            end
            
            % iv. upperenvelope and interp to common
            [sol(2,t).lacon, par] = fun.upperenvelopetocommon(par,lacon,interp);
            
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % i. constrained: a = 0, d = 0, l = lmax (fullcon) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if par.lmax > 0 
                    
            fullcon.time = tic; 
            fullcon.num  = 8;
                    
            % i. choices
            sol(2,t).fullcon.c = par.grid_m_nd + par.rk*par.grid_k_nd*par.lmax;
            sol(2,t).fullcon.d = zeros(par.Nm,par.Nn,par.Nk);
            sol(2,t).fullcon.l = par.lmax*ones(par.Nm,par.Nn,par.Nk);
            
            a_fullcon = zeros(par.Nm,par.Nn,par.Nk);
            b_fullcon = par.grid_n_nd;
            q_fullcon = (1-par.delta)*par.grid_k_nd + sol(2,t).fullcon.l;
            
            % ii. post-decision states and value
            sol(2,t).fullcon.v = fun.u(sol(2,t).fullcon.c,work,sol(2,t).fullcon.l,par) +  par.beta*...
                interp.w(a_fullcon,b_fullcon,q_fullcon);
            
            sol(2,t).fullcon.time_upper = 0;
            if par.print == 1
                fprintf('   fullcon: %4.1f secs\n',toc(fullcon.time));
            end
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%
        % j. upper envelope %
        %%%%%%%%%%%%%%%%%%%%%
        
        if par.print == 1
            fprintf('  after segments:\n');        
            t2 = tic;
        end
        
        % i. stack segments
        if par.ndim == 3
            
            vars = {'v','c','d','l'};
            if par.lmax > 0
                segs = {'ucon','dcon','acon','con','lcon','ldcon','lacon','fullcon'};
            else
                segs = {'ucon','dcon','acon','con'};
            end
            numsegs = numel(segs);
            vs = struct();
            for i = 1:numel(vars)
                vs.(vars{i}) = NaN(par.Nm,par.Nn,par.Nk,numsegs);
                for j = 1:numsegs
                    vs.(vars{i})(:,:,:,j) = sol(2,t).(segs{j}).(vars{i});
                end
            end

            % ii. find max
            [~,I] = max(vs.v,[],4);

            % iii. index for maximum
            [I1, I2, I3] = ndgrid(1:par.Nm,1:par.Nn,1:par.Nk);
            index = sub2ind(size(vs.v),I1(:),I2(:),I3(:),I(:));

            % iv. over-arching optimal choices
            for i = 1:numel(vars)
                sol(2,t).(vars{i}) = reshape(vs.(vars{i})(index),[par.Nm,par.Nn,par.Nk]);
            end
            sol(2,t).v = fun.trans(sol(2,t).v,par);
        
        else
            
            vars = {'v','c','d'};
            vs = struct();
            for i = 1:numel(vars)
                vs.(vars{i}) = NaN(par.Nm,par.Nn,4);
                vs.(vars{i})(:,:,1) = sol(2,t).ucon.(vars{i});
                vs.(vars{i})(:,:,2) = sol(2,t).con.(vars{i});
                vs.(vars{i})(:,:,3) = sol(2,t).acon.(vars{i});
                vs.(vars{i})(:,:,4) = sol(2,t).dcon.(vars{i});
            end

            % ii. find max
            [~,I] = max(vs.v,[],3);

            % iii. index for maximum
            [I1, I2] = ndgrid(1:par.Nm,1:par.Nn);
            index = sub2ind(size(vs.v),I1(:),I2(:),I(:));

            % iv. over-arching optimal choices
            for i = 1:numel(vars)
                sol(2,t).(vars{i}) = reshape(vs.(vars{i})(index),[par.Nm,par.Nn]);
            end
            sol(2,t).v = fun.trans(sol(2,t).v,par);

        
        end
        if par.print == 1
            fprintf('   max over common: %3.1f secs\n',toc(t2));
        end
        

        %%%%%%%%%%%%%%%%%%
        % k. derivatives %
        %%%%%%%%%%%%%%%%%%
        
        t2 = tic;        
        
        interp.vm_next_work.Values = fun.trans(fun.marg_u(sol(2,t).c,work,par),par);       
                
        if par.ndim == 3
            a = par.grid_m_nd + par.rk*par.grid_k_nd.*sol(2,t).l - sol(2,t).c - sol(2,t).d;
            b = par.grid_n_nd + sol(2,t).d + fun.f_pens(sol(2,t).d,par);
            q = (1-par.delta)*par.grid_k_nd + sol(2,t).l;
            interp.vn_next_work.Values = fun.trans(par.beta*interp.wb(a,b,q),par);
            interp.vk_next_work.Values = fun.trans(fun.marg_u(sol(2,t).c,work,par)*par.rk.*sol(2,t).l + par.beta*(1-par.delta)*interp.wq(a,b,q),par);
        else
            a = par.grid_m_nd - sol(2,t).c - sol(2,t).d;
            b = par.grid_n_nd + sol(2,t).d + fun.f_pens(sol(2,t).d,par);            
            interp.vn_next_work.Values = fun.trans(par.beta*interp.wb(a,b),par);
        end
        
        if par.print == 1
            fprintf('   derivatives: %3.1f secs\n',toc(t2));         
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%
        % l. discrete choice %
        %%%%%%%%%%%%%%%%%%%%%%
        
        t2 = tic;
        [sol] = fun.discrete_choice(par,sol,t);
        if par.print == 1       
            fprintf('   discrete choice: %3.1f secs\n',toc(t2));         
        end
        
        
        %%%%%%%%%%
        % m. end %
        %%%%%%%%%%
        
        if par.print == 1
            fprintf('  total: %4.1f secs\n',toc(t1));
            time_tot_upper = sol(2,t).ucon.time_upper...
                             + sol(2,t).dcon.time_upper...
                             + sol(2,t).acon.time_upper...
                             + sol(2,t).con.time_upper;
            fprintf('   upper: %3.1f secs\n',time_tot_upper);   
            time_tot_tot_upper = time_tot_tot_upper + time_tot_upper;
            mem = memory;         
            fprintf('  Memory: %2.1f GB\n',mem.MemUsedMATLAB/10^9);
        end
        
            % delete segment if not to be saved
            if par.save_segments == 0
                sol(2,t).ucon = [];
                sol(2,t).con = [];
                sol(2,t).acon = [];
                sol(2,t).dcon = [];
                if par.ndim == 3
                    sol(2,t).lcon = [];
                    sol(2,t).lacon = [];
                    sol(2,t).ldcon = [];
                    sol(2,t).fullcon = [];                    
                end
            end
    end
    
    if par.print == 1
        fprintf(' solution time: %4.2f mins (%4.0f secs)\n',...
            toc(t0)/60,...
            toc(t0));    
        fprintf('  total upper: %3.1f%%\n',...
            time_tot_tot_upper/toc(t0)*100); 
    end

end

%% over and out
end
end