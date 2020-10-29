function [opt] = example_3()

%% the data is from MNIST Handwritten Digits
%  the number of samples: 1010
%%

clear;

load('mnist_test3.mat');

%% Set Initialization

load('e_3c0.mat');

ns = 160;

%% Compute Wasserstein Barycenter (PAM)

[opt] = PAM_start(mnist_test3,c0,ns);