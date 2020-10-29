function [opt] = example_2()

%% the data is from USPS Handwritten Digits
%  the number of samples: 1100
%%

clear;

load('usps1.mat');

%% Set Initialization

load('e_2c0.mat');

ns = 80;

%% Compute Wasserstein Barycenter (PAM)

[opt] = PAM_start(usps1,c0,ns);

