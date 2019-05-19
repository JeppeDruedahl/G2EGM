classdef egmfun
% contains: functions for solving the model
methods (Static)
   
    function k = inv_k(q,l,par)
        k = (q-l)/(1-par.delta);
    end
    
    function m = inv_m(a,c,d,par)
        m = a + c + d;
    end
    
    function m = inv_m_3d(a,k,l,c,d,par)
        m = a + c + d - par.rk*k.*l ;
    end
    
    function n = inv_n(b,d,par)
        n = b-d-fun.f_pens(d,par);
    end
    
    function [v,k,m,n] = inv_kmn_and_v(a,b,q,c,d,l,w,work,par)
        v = fun.u(c,work,l,par) + par.beta*w;
        if par.ndim == 3           
            k = egmfun.inv_k(q,l,par);
            m = egmfun.inv_m_3d(a,k,l,c,d,par);
            n = egmfun.inv_n(b,d,par);
        else
            k = [];
            m = egmfun.inv_m(a,c,d,par);
            n = egmfun.inv_n(b,d,par);            
        end
    end
    
    function c = foc_c(wa,work,par)
        c = fun.inv_marg_u(par.beta*wa,work,par);
    end
    
    function d = foc_d(wa,wb,par)
        d = (par.chi.*wb)./(wa-wb) - 1;
    end
    
    function d = foc_d_ccon(c,wb,work,par)
        d = par.chi./(fun.marg_u(c,work,par)./(par.beta*wb)-1)-1;
    end
    
    function l = foc_l(q,wa,wq,par)
        lfac = (par.beta*par.rk)/(1-par.delta)*wa;
        l    = (lfac.*q + par.beta*wq)./(lfac+par.varphi);
    end
    
    function l = foc_l_ccon(c,q,wq,work,par)
        lfac = fun.marg_u(c,work,par)*par.rk/(1-par.delta);
        l    = (lfac.*q + par.beta*wq)./(lfac+par.varphi);
    end
    
    function [c_min, c_max] = c_min_max_acon(wb,work,par)
        c_min = fun.inv_marg_u((par.chi+1)*par.beta*wb,work,par);
        c_max = fun.inv_marg_u(par.beta*wb,work,par);
    end
    
    function v_x = deviate_d_con(a,q,seg,work,par,interp)
        
        % choices
        d_x = par.delta_con*seg.c; % > 0
        c_x = (1-par.delta_con)*seg.c; % < c            
            
        % post-decision states            
        b_x = seg.n + d_x + fun.f_pens(d_x,par);
        
            I = imag(b_x) ~= 0;
            b_x(I) = NaN;
            
        if par.ndim == 3
            w_x = interp.w(a,b_x,q);
        else
            w_x = interp.w(a,b_x);            
        end

        % value
        v_x = fun.u(c_x,work,seg.l,par) + par.beta*w_x;
        
            I = imag(v_x) ~= 0;
            v_x(I) = NaN;   
            
    end
    
    function v_x = deviate_a_con(b,q,seg,work,par,interp)
        
        % choices
        c_x = (1-par.delta_con)*seg.c;            

        % post-decision states
        a_x = par.delta_con*seg.c;
        if par.ndim == 3
            w_x = interp.w(a_x,b,q);
        else
            w_x = interp.w(a,b);
        end
            
        % value
        v_x = fun.u(c_x,work,seg.l,par) + par.beta*w_x;
            
    end
    
    function v_x = deviate_l_con(a,b,seg,work,par,interp)
        
        % choices
        l_x = (1-par.delta_con)*seg.l;            
        c_x = seg.c - par.rk*seg.k.*par.delta_con.*seg.l;
        
        % post-decision states
        q_x = (1-par.delta)*seg.k + l_x;
        w_x = interp.w(a,b,q_x);
            
        % value
        v_x = fun.u(c_x,work,l_x,par) + par.beta*w_x;
            
    end
    
end
end