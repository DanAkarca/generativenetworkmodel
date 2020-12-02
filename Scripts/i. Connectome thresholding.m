%% Connectome thresholding
%{
Akarca D, et al. A generative network model of neurodevelopment.
Code written by Dr Danyal Akarca, MRC Cognition and Brain Sciences Unit.
University of Cambridge: https://www.neuroscience.cam.ac.uk/directory/profile.php?da434
Email: danyal.akarca@mrc-cbu.cam.ac.uk
2nd December 2020.
%}
%% About this script
%{
This script replicates the basic preprocessing steps for thresholding
connectomes. After you have set the appropriate paths (below) all code
should run without having to make further edits.

Due to NHS data restrictions, we have not shared observed connectome data. 
Instead, we have provided simualted connectomes using Maslow-Sneppen rewiring
with 100 rewires per connection, to replicate basic statistical properties.
%}
%% 0. set paths and load required data

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
load('example_rewired_connectomes');
%% 1. threshold connectomes
% number of subjects
nsub = 270;
% set threshold
thr  = 27;
% initialise
density               = [];
binarised_connectomes = [];
% loop through subjects
W_thr = [];
for sub = 1:nsub;
    W     = squeeze(example_rewired_connectomes(sub,:,:));          % take unthresholded connectome
    W_thr = threshold_absolute(W,thr);                              % absolute threshold
    density(sub) = density_und(W_thr);                              % calculate density
    B = zeros(size(W_thr));
    B(find(W_thr)) = 1;                                             % binarise
    binarised_connectomes(sub,:,:) = B;
end
% display
disp(sprintf('At a threshold of 27 streamlines, a mean density of %g%% is produced across the sample.',100*mean(density)));
%% 2. visualise thresholding
% set subject to visualise threhsolding
subVis = 1;
% visualise
figure; set(gcf,'Position',[100 100 2000 500]);
subplot(1,3,1);
histogram(density); title('Example density across the sample'); 
xlabel('Density'); ylabel('Frequency'); set(gca,'TickLength',[0 0]);
subplot(1,3,2);
imagesc(squeeze(example_rewired_connectomes(subVis,:,:))); title(sprintf('Subject %g: Unthresholded network',subVis));
xlabel('Node i'); ylabel('Node j'); set(gca,'TickLength',[0 0]);
subplot(1,3,3);
imagesc(squeeze(binarised_connectomes(subVis,:,:))); title(sprintf('Subject %g: Thresholded and binarised networks',subVis));
xlabel('Node i'); ylabel('Node j'); set(gca,'TickLength',[0 0]);