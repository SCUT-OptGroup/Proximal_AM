%% ***************************************************************
%% Filename: PAM
%% ***************************************************************
%%
% Proximal alternating minimization method for calculating a Wasserstein Barycenter 
%% ***************************************************************

function [opt,iter,ttime,obj,pinf,lambda,par,value,lptime,KKT_res,ssubiter] = PAM(c0,A,b,ns,bdim,OPTIONS_PAM)

%% *************** parameters for PAM ***************************

if isfield(OPTIONS_PAM,'tol');          tol      = OPTIONS_PAM.tol;       end
if isfield(OPTIONS_PAM,'printyes');     printyes   = OPTIONS_PAM.printyes;   end
if isfield(OPTIONS_PAM,'maxiter');      maxiter    = OPTIONS_PAM.maxiter;    end

if (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n ************** PAM  for calculating a Wasserstein Barycenter  ****************');
    fprintf('\n ****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter       obj          alpha          infeas          pinf       subiter      time');
end
       
OPTIONS_ADMM.printyes = 0;
OPTIONS_ADMM.maxiter = 5000;

%% ************** Initialization part for BCD *******************

N = length(bdim);

par.ns = ns;

par.bdim = bdim;

par.normb = norm(b);

rho = 1e-5;           % the proximal parameter for x

rhoN = rho*N/2;

%% ******************* Initialization part ***********************
e_vec = ones(ns,1);

Z = zeros(ns,sum(bdim));

y = c0{1,1}.w';

lambda = zeros(ns,N);

X = c0{1,1}.supp;

W = Fmap(X,A); 

meanW = mean(mean(W));

par.alpha = 100;%0.05*meanW;           % the proximal parameter for (Z,y)

par.beta = (ns/(N*max(bdim)))*meanW;

par.sn = sum(bdim);

par.normA = norm(A,'fro');

OPTIONS_ADMM.tol = min(12/mean(bdim),1);%min(1/(par.normb),0.1); 

obj_list = zeros(maxiter,1);

ssubiter = 0;

flag = 0;

printflag = 1;

tstart = clock;

%% ************************ Main Loop ***************************

for iter = 1:maxiter
    
    %% ************** solve the subproblem of (Z,y) ****************
    
    [Zk,yk,lambda,beta,pinf,Amaplambda,res2,infeas,subiter] = sPADMM(Z,y,lambda,par,OPTIONS_ADMM,W,b);
    
    ssubiter = ssubiter + subiter;
    
    SZk = sum(Zk,2);
    
    %****************** solve the subproblem of X *****************************

    Xk = (A*Zk'+rhoN*X)./(SZk+rhoN*e_vec)';

    Wk = Fmap(Xk,A);
    
    obj = (1/N)*sum(sum(Zk.*Wk));
    
    obj_list(iter) = obj;

    %*************************** stopping condition ***************************
    
    diff_Z = Zk - Z; diff_y = yk - y;
    
    infeas = infeas/(1+par.normb);
        
    OPTIONS_ADMM.tol = max(1e-5,0.8*OPTIONS_ADMM.tol);
    
    ttime = etime(clock,tstart);
    
    if (printyes)
        
        fprintf('\n %3.0d     %5.4e      %3.2e       %3.2e        %3.2e       %3.0d       %3.2f',iter,obj,par.alpha,infeas,pinf,subiter,ttime);
        
    end
    
    if iter>=2 && (abs(obj-obj_list(iter-1))/max(abs(obj),1)<1e-4)
        
        flag = flag + 1;
        
    else
        
        flag = 0;
        
    end
    
    if infeas<tol
        
        diff_X = Xk - X;
        
        temp_res1 = Zk - Wk/N - Amaplambda;
        
        temp_res1 = Simplex_matcol(temp_res1,b);
        
        res1 = norm(Zk - temp_res1,'fro')/(1+par.normA);
        
        res2 = res2/(1+par.normb);
        
        res3 = norm(rho*diff_X,'fro')/(1+par.normA);
        
        KKT_res = max([infeas,res1,res2,res3]);
        
        if (printyes)
            
            if printflag == 1
                
                fprintf('\n  iter       obj          alpha          KKT_res          pinf       subiter      time');
                
                printflag = 0;
                
            end
            
            fprintf('\n %3.0d     %5.4e      %3.2e       %3.2e        %3.2e       %3.0d       %3.2f',iter,obj,par.alpha,KKT_res,pinf,subiter,ttime);
            
        end
        
        if  (iter>=5&&KKT_res<tol) || (iter>=30&&flag>=10)
            
            Xopt = Xk; yopt = yk;
            
            opt{1}.supp = Xopt;
            
            opt{1}.w = yopt';
            
            opt{1}.Z = Zk;
            
            break;
            
        end
        
    end
    
    %% ************* update the parameter alpha_k and beta *******************
    
     if (0.5*par.alpha*norm(diff_Z,'fro')^2+norm(diff_y)^2>1.0e-5*obj)

         par.alpha = max(1.0e-4,par.alpha/2);  
         
     end
            
    par.beta = beta;
    
    X = Xk;  Z = Zk;  W = Wk;  y = yk;    
end
%% ***************************************************************

if (iter == maxiter)
    
    Xopt = Xk; yopt = yk;
    
    opt{1}.supp = Xopt;
    
    opt{1}.w = yopt';
    
    opt{1}.Z = Zk;
    
end
 
c.supp = Xopt;

c.w = yopt';

tstart = clock;

[Zopt,value] = centroid_meancost(bdim', A, b', c);

lptime = etime(clock,tstart);

opt{1}.Z = Zopt;

end
