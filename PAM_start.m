function [opt,iter] =PAM_start(db,c0,ns)
%% db
%   db{1}.supp   --- a matrix whose column is a support vector
%   db{1}.w      --- weight of corresponding support vector
%   db{1}.stride --- dimensions of corresponding weight
%% c0 is a start point of PAM
%   c0{1}.supp
%   c0{1}.w
%% ns --- the number of support vectors in Wasserstein barycenter 
%%

OPTIONS_PAM.tol = 5e-4;      %accuracy tolerance of PAM
OPTIONS_PAM.printyes = 1;
OPTIONS_PAM.maxiter = 100;    %maximum number of PAM 

bdim = db{1,1}.stride';

b = db{1,1}.w';

A = db{1,1}.supp;

[opt,iter] = PAM(c0,A,b,ns,bdim,OPTIONS_PAM);