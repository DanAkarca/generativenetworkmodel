%% Running the initial generative models
%{
Akarca D, et al. A generative network model of neurodevelopment.
Code written by Dr Danyal Akarca, MRC Cognition and Brain Sciences Unit.
University of Cambridge: https://www.neuroscience.cam.ac.uk/directory/profile.php?da434
Email: danyal.akarca@mrc-cbu.cam.ac.uk
2nd December 2020.
%}
%% About this script
%{
This script runs generative models, as defined by Vertes et al. 2012 & Betzel et al. 2016,
across thirteen generative rules. Note, this script was adapted to be run
on our internal cluster.

Due to NHS data restrictions, we have not shared observed connectome data. 
Instead, we have provided simualted connectomes using Maslow-Sneppen rewiring
with 100 rewires per connection, to replicate basic statistical properties.
%}
%% %%% Step 1: Load the relevant data %%%
% clear workspace and command window
clear; clc;
% set latex formatting 
set(0,'DefaultTextInterpreter','Latex',...
    'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',1,...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
set(groot,'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');  
% add the brain connectivity toolbox <<<<<<<<<<<<< SET
brainconnectivity = '/Users/da04/Desktop/PhD/Research/Toolboxes/BCT/2019_03_03_BCT';
addpath(brainconnectivity);
% set the current directory to the data <<<<<<<<<<<<< SET
datapath = '/Users/da04/Desktop/Generative Models GitHub/Example data/';
cd(datapath);
% load rewired data
load('example_binarised_connectomes');
% load coordinates
load('dk_coordinates.mat');
% load example euclidean distances
load('example_euclidean.mat');
%% 2. Set up the generative model
% initialise number of subjects
nsub         = 270;
% Define target connectomes and compute seed
Atgt_set    = example_binarised_connectomes;
proportion  = 0.2; % set proportion of edges in common across subjects (this was set to 1 for our work i.e. 100% edges common to all subjects)
connections = squeeze(mean(Atgt_set,1));
index       = find(connections==proportion);
A           = zeros(size(connections));
A(index)    = 1; A = upper(A);
D = squareform(pdist(coordinates));
modeltype = {...
    'sptl',...          
    'neighbors',...     
    'matching',...      
    'clu-avg',...       
    'clu-min',...       
    'clu-max',...       
    'clu-diff',...      
    'clu-prod',...      
    'deg-avg',...       
    'deg-min',...       
    'deg-max',...       
    'deg-diff',...      
    'deg-prod'};        
modelvar  = [{'powerlaw'},{'powerlaw'}];
% define parameters
etalimits = [-7 7]; 
gamlimits = [-7 7];
nruns     = 64; % we ran 10000
[p,q]     = meshgrid(linspace(etalimits(1),etalimits(2),sqrt(nruns)),...
                   linspace(gamlimits(1),gamlimits(2),sqrt(nruns)));
params    = unique([p(:) q(:)],'rows'); 
nparams   = size(params,1);
%% Step 3: Run the generative models
%{
Commented out code was used to run each script on our internal cluster.
In this case, we ran a single subject per script. If you would like to do
the same, please edit the below loops over subjects accordingly.

path = mfilename('fullpath');
sub1 = path(end); sub1_num = str2double(sub1);
sub2 = path(end-1); sub2_num = str2double(sub2);
sub3 = path(end-2); sub3_num = str2double(sub3);
if isnan(sub3_num);
    if isnan(sub2_num);
        w = sub1_num;
    else
        w = strcat(sub2,sub1);
        w = str2double(w);
    end
else
    w = strcat(sub3,sub2,sub1);
    w = str2double(w);
end
%}
% initialise 
generativedata  = {};
Asynthall      = {};
Eall           = {};
Kall           = {};
% loop over subjects
for w = 1:2
    tic
    Atgt = squeeze(Atgt_set(w,:,:));
    m    = nnz(Atgt)/2;
    n    = length(Atgt);
    x    = cell(4,1);
    x{1} = sum(Atgt,2);
    x{2} = clustering_coef_bu(Atgt);
    x{3} = betweenness_bin(Atgt)';
    x{4} = D(triu(Atgt,1) > 0);
    disp(sprintf('Running generative models for subject %g...',w));
    % loop over models 
    for v = 1:13
        B    = generative_model(A,D,m,modeltype{v},modelvar,params);
        nB = size(B,2);
        K  = zeros(nB,4);
        for iB = 1:nB
            b = zeros(n);
            b(B(:,iB)) = 1;
            b = b + b';
            Asynthall{w}(v,iB,:,:) = b; % keep synthetic networks: note, this is computationally expensive
            y = cell(4,1);
            y{1} = sum(b,2);
            y{2} = clustering_coef_bu(b);
            y{3} = betweenness_bin(b)';
            y{4} = D(triu(b,1) > 0);
            for j = 1:4
                K(iB,j) = fcn_ks(x{j},y{j});
            end
            Kall{w}(v,:,:) = K; % keep ks statistics
            Eall{w}(v,:)   = max(K,[],2);
        end
    end
    t = toc
    disp(sprintf('Subject %g complete (%g seconds)',w,round(t,2)));
end
generativedata = {Asynthall Eall Kall}; % single cell to save all the data: subject by parameter by node by node space.
%% Define KS function
function kstat = fcn_ks(x1,x2)
binEdges    =  [-inf ; sort([x1;x2]) ; inf];
binCounts1  =  histc (x1 , binEdges, 1);
binCounts2  =  histc (x2 , binEdges, 1);
sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);
sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);
deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
kstat = max(deltaCDF);
end