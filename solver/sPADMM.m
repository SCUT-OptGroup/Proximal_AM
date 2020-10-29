%% ***************************************************************
%% ******************* Filename: sPADMM ************************
%% 
%% This code is the semi-proximal ADMM with a large step-size for solving 
%
%  min_{Z,y} (1/N) sum_{s=1}^N (<Z^s,F^s(x^k)> + h(Z^s) + (alpha/2)||Z^s - Z^{s,k}||_F^2) 
%             + deta_{Deta}(y) + (alpha/2)||y-y^k||^2
%
%  s.t.  Z^s e - y =0, s = 1,2,...,N
%
%% alpha: is the proximal parameter 
%% Deta:  is the simplex set
%% ***************************************************************


function [Zopt,yopt,lambda_opt,betak,pinf,Amaplambda,PAM_measure2,PAM_inf,iter,ttime] = sPADMM(Z,y,lambda,par,OPTIONS_ADMM,W,b)

%% *************** parameters for ADMM ***************************

if isfield(OPTIONS_ADMM,'tol');          tol        = OPTIONS_ADMM.tol;        end
if isfield(OPTIONS_ADMM,'printyes');     printyes   = OPTIONS_ADMM.printyes;   end
if isfield(OPTIONS_ADMM,'maxiter');      maxiter    = OPTIONS_ADMM.maxiter;    end

prim_win = 0;

dual_win = 0;

if (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n ************** The ADMM for solving the subproblem of Zypart  ****************');
    fprintf('\n ****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter     primal_inf       relgap          beta       max(sigma)  time');
end

bdim = par.bdim;

N = length(bdim);

ns = par.ns;

beta = par.beta;

alpha = par.alpha;

normb = par.normb;

normA = par.normA;

sn = par.sn;

tau = 1.618;          % the convergent factor for the ADMM !!
%% ************** Initialization part for ADMM ******************

Z0 = Z;  y0 = y;  

alpha_Z0 = alpha*Z0; alpha_y0 = alpha*y0;

WdN = W/N;  

norm_WdN = norm(WdN,'fro');  

Z0_WdN = alpha_Z0 - WdN;

slambda = sum(lambda,2); 

y_repN = repmat(y,1,N);

tstart = clock;

%% ************************ Main Loop ***************************

for iter = 1:maxiter

    sigma_temp = beta*bdim;
    
    ybeta_rep = beta*y_repN;
    
    yblam_mat = ybeta_rep - lambda;

%********************** solve the subproblem of Zpart *********************

    temp_Zk = BCD_Zpart(Z,yblam_mat,bdim,sigma_temp,alpha,beta,Z0_WdN);

    Zk = Simplex_matcol(temp_Zk,b);

   %% ****************** solve the subproblem of ypart *****************
    
    beta_SZk = beta*sum(Zk,2);
   
    betaN_alpha =1/(beta*N+alpha);
                
    temp_y = betaN_alpha*(beta_SZk + slambda + alpha_y0);
    
    dtemp_y = sort(temp_y,'descend');
    
    proj_temp_y = simplex_y(dtemp_y,1);

    yk = max(temp_y-repmat(proj_temp_y,ns,1),0);
    
    %***********************  update the Lagrange multiplier ****************
    
    y_repN = repmat(yk,1,N);
    
    tau_beta = tau*beta;
    
    diff_lambda = tau_beta*(Atmap(N,bdim,Zk) - y_repN);
    
    lambdak = lambda + diff_lambda;
    
    %*************************** stopping condition ***************************
    
    diff_Z = Zk - Z;    diff_y = yk - y; 
    
    Z = Zk;              y = yk;  
     
    lambda = lambdak;    slambda = sum(lambda,2);
    
    ndiff_lambda = norm(diff_lambda,'fro')/tau_beta;
    
    primal_inf = ndiff_lambda/(1+normb);
    
    beta_update_iter = beta_fun(iter);
    
    if primal_inf<tol || rem(iter,beta_update_iter)==0%mod(iter,10)==0
    
    Amaplambda = Amap(sn,bdim,lambda);

    temp_dinf1 = (1-alpha)*Zk - Amaplambda + Z0_WdN; 
    
    temp_dinf1 = Simplex_matcol(temp_dinf1,b);
    
    dinf1 = Zk - temp_dinf1; 
    
    temp_e1 = slambda + (1-alpha)*yk + alpha_y0;
    
    dtemp_e1 = sort(temp_e1,'descend');
    
    proj_temp_e1 = simplex_y(dtemp_e1,1);
    
    temp_e1 = max(temp_e1-repmat(proj_temp_e1,ns,1),0);
    
    dinf2 = yk - temp_e1;
    
%     dinf1 = dinfpart(diff_Z,diff_y,diff_lambda,tau,sigma_temp,beta,alpha,N,bdim);
%     
%     dinf2 = ((1-tau)/tau)*sum(diff_lambda,2);
    
    norm_dinf1 = norm(dinf1,'fro');
    
    norm_dinf2 = norm(dinf2);

    dual_inf = max(norm_dinf1/(1+normA),norm_dinf2/(1+normb));%sqrt(norm_dinf1^2 +norm_dinf2^2)/(1+normA);
    
    KKT_res = sqrt(ndiff_lambda^2 + norm_dinf1^2 + norm_dinf2^2)/(1+norm_WdN);
    
    ratio = primal_inf/dual_inf;
    
    ttime = etime(clock,tstart);
    
    if (printyes) && (mod(iter,10)==0)
        
        fprintf('\n %3.0d        %3.2e        %3.2e       %3.2e     %3.2e   %3.2f',iter,primal_inf,KKT_res,beta,max(sigma_temp),ttime);
        
    end
    
    norm_Zk0 = norm(Zk-Z0,'fro');
    
    norm_yk0 = norm(yk-y0);
    
    Xi = dinfpart(diff_Z,diff_y,diff_lambda,tau,sigma_temp,beta,alpha,N,bdim);
    
    xi = ((1-tau)/tau)*sum(diff_lambda,2);
    
%     if max(primal_inf,dual_inf)<tol
%         
%         exit = 1;
%         
%     end
    
    if  max(primal_inf,dual_inf)<tol || (iter>100&&(norm(Xi,'fro')<0.25*alpha*norm_Zk0&&norm(xi)<0.25*alpha*norm_yk0))
        
        Zopt = Z;  yopt = y;
        
        lambda_opt = lambda;  betak = beta;
        
        pinf = max(sum(diff_lambda.*diff_lambda,1).^(1/2)/tau_beta)/(1+normb);
        
        temp_e1 = slambda + yk;
        
        dtemp_e1 = sort(temp_e1,'descend');
        
        proj_temp_e1 = simplex_y(dtemp_e1,1);
        
        temp_e1 = max(temp_e1-repmat(proj_temp_e1,ns,1),0);
        
        inf2 = yk - temp_e1;
        
        PAM_measure2 = norm(inf2);
        
        PAM_inf = ndiff_lambda;
        
        break;
        
    end
    
    %% ******************* to update the value of beta *******************

    betascale = 1.2;
    
    if (ratio <1)
        
        prim_win = prim_win+1;       % the primal feasibility is better 
        
    elseif (ratio>1)
        
        dual_win = dual_win+1;
        
    end
  
    if (rem(iter,beta_update_iter)==0)
    
        betamax = 1e4;  betamin = 1e-4;
        
        if (iter <= 5000)
            
            if (prim_win > dual_win)%max(1,1.2*dual_win)
                
                prim_win = 0;
                
                beta = max(betamin,beta/betascale);
                
            elseif (dual_win > prim_win)%max(1,1.2*prim_win)
                
                dual_win = 0;
                
                beta = min(betamax,beta*betascale);
                
            end            
        end        
    end
    
    end
  
end

if (iter == maxiter)
    
     Zopt = Z;  yopt = y; lambda_opt = lambda; betak = beta;
     
     pinf = max(sum(diff_lambda.*diff_lambda,1).^(1/2)/tau_beta)/(1+normb);
     
     temp_e1 = slambda + yk;
     
     dtemp_e1 = sort(temp_e1,'descend');
     
     proj_temp_e1 = simplex_y(dtemp_e1,1);
     
     temp_e1 = max(temp_e1-repmat(proj_temp_e1,ns,1),0);
     
     inf2 = yk - temp_e1;
     
     PAM_measure2 = norm(inf2);
     
     PAM_inf = ndiff_lambda;
end

end

%%
%% ***********************************************************************
  function sigma_update_iter = beta_fun(iter)
         
       if (iter < 30)
            sigma_update_iter = 30;
        elseif (iter < 250)
            sigma_update_iter = 25;
        elseif (iter < 500)
            sigma_update_iter = 50;
        else
            sigma_update_iter = 100;
        end
 
  end


%% ************* End of the main program *************************