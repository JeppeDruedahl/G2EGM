classdef funfig_3d
methods (Static)
    
function figs = plot(figs,par,sol,t,ik,in,method,varlist)

    par = fun.solprep(par);    

    % loop over variables
    for i = 1:numel(varlist)
        
        % mn
        var = varlist{i};
        
        m = fun.vec(par.grid_m_nd(:,:,ik));
        n = fun.vec(par.grid_n_nd(:,:,ik));
        I = m < par.fig_max_m & n < par.fig_max_n;
        y = fun.vec(sol(2,par.T-t).(var)(:,:,ik));
                  
        name = ['T' num2str(t,'%02.0f') '_' var '_' method '_mn_3d'];    
        [figs.(var).(method).mn, ~] = fun.myfigure(par,name);
            scatter3(m(I),n(I),y(I),4,y(I),'filled')
            if strcmp(var,'q') || strcmp(var,'l') 
                view([-37.5 50]);
                if strcmp(var,'l') 
                    zlim([0 par.lmax])
                else
                    zlim((1-par.delta)*par.grid_k(ik)+[0 par.lmax])
                end
            else
                view([-37.5 30]);
            end
            grid on;
            
        % mk
        m = fun.vec(par.grid_m_nd(:,in,:));
        k = fun.vec(par.grid_k_nd(:,in,:));
        I = m < par.fig_max_m & k < par.fig_max_k & k > par.fig_min_k;
        y = fun.vec(sol(2,par.T-t).(var)(:,in,:));
                  
        name = ['T' num2str(t,'%02.0f') '_' var '_' method '_mk_3d'];    
        [figs.(var).(method).mk, ax] = fun.myfigure(par,name,'m','k');
            scatter3(m(I),k(I),y(I),4,y(I),'filled')
            view([-37.5 30]);
            grid on;
            ylim([0 par.fig_max_k])
            ax.YTick = 0:5:par.fig_max_k;
            if strcmp(var,'l') 
                zlim([0 par.lmax])
            end           
    end
end

%% over and out
end
end