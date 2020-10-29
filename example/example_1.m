function [opt] = example_1()

%% Load data
clear;
max_sample_size = 1000;

randstate = 1;
randn('state',double(randstate));
rand('state',double(randstate));

s_modalities = 2;
d_modalities = [3 3];
filename = 'mountaindat.d2';
db_tmp = loaddata(max_sample_size, s_modalities, d_modalities, filename);

db=cell(1,1); db{1}=db_tmp{1}; % only use the color quantization data

%% Set Initialization

options.init_method = 'kmeans'; % {'kmeans', 'mvnrnd'}
options.support_size = 60;
options.max_support_size = options.support_size;

c0 = cell(1,1);

c0{1} = centroid_rand(db{1}.stride, db{1}.supp, db{1}.w, options);

ns = 60;

%% Compute Wasserstein Barycenter (PAM)

[opt] = PAM_start(db,c0,ns);

