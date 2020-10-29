function [opt] = example_4()

%% the data is from BBC News dataset
%  the number of samples: 510
%%

clear;

load('bbc_1.mat');

%% Set Initialization

load('e_4c0.mat');

ns = 25;

%% Compute Wasserstein Barycenter (PAM)

[opt] = PAM_start(bbc_1,c0,ns);

