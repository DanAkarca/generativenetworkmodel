%% Computing the seed network
%{
Akarca D, et al. A generative network model of neurodevelopment.
Code written by Dr Danyal Akarca, MRC Cognition and Brain Sciences Unit.
University of Cambridge: https://www.neuroscience.cam.ac.uk/directory/profile.php?da434
Email: danyal.akarca@mrc-cbu.cam.ac.uk

2nd December 2020.
%}
%% About this script
%{
This script replicates how the seed network were computed for running of the
generative models across our sample. After you have set the appropriate paths 
(below) all code should run without having to make further edits.

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
load('example_binarised_connectomes');
load('dk_coordinates');
%% 1. Determine the seed network
% note, as these are rewired simulated connectomes, there are less edges in
% common that you would expect from observed connectomes.
Atgt_set    = example_binarised_connectomes;
proportion  = 0.2; % set proportion of edges in common across subjects (this was set to 1 for our work i.e. 100% edges common to all subjects)
connections = squeeze(mean(Atgt_set,1));
index       = find(connections==proportion);
A           = zeros(size(connections));
A(index)    = 1; A = upper(A);
%% 2. Compute characteristics of the seed for visualisation
degree      = degrees_und(A);
%% 3. Visualise
figure; 
set(gcf,'position',[0 100 2000 500]);
subplot(1,3,1); imagesc(A); title('Seed network adjacency matrix'); xlabel('Node i'); ylabel('Node j');
subplot(1,3,2); imagesc(connections); title('Proportion of edges conserved across subjects'); 
xlabel('Node i'); ylabel('Node j'); 
c = colorbar; c.Label.String = 'Proportion of shared connections'; 
c.Label.Interpreter = 'latex'; c.TickLabelInterpreter = 'latex';
subplot(1,3,3); h = plot(graph(A),...
    'LineWidth',2,...
    'NodeColor','k',...
    'EdgeAlpha',0.75,...
    'XData',coordinates(:,1),...
    'YData',coordinates(:,2),...
    'ZData',coordinates(:,3));
labelnode(h,1:68,''); 
title('Example seed network visualisation'); xlabel('x'); ylabel('y'); zlabel('z');
for node = 1:68
    highlight(h,node,'MarkerSize',3*degree(node)+2);
end
hold on;
B = squeeze(example_binarised_connectomes(120,:,:));
C = B-A; 
h = plot(graph(C),...
    'LineWidth',2,...
    'NodeColor','k',...
    'EdgeColor','r',...
    'EdgeAlpha',0.2,...
    'XData',coordinates(:,1),...
    'YData',coordinates(:,2),...
    'ZData',coordinates(:,3));
labelnode(h,1:68,'');
k = (1.5*degrees_und(B)+2);
for node = 1:68
    highlight(h,node,'MarkerSize',k(node));
end
%% 4. Number of connections
nz = [];
for sub = 1:270
    nz(sub) = nnz(squeeze(example_binarised_connectomes(sub,:,:)))/2;
end
figure; 
histogram(nz,'facecolor','b','facealpha',0.6); xlim([0 350]); 
title('Connection distribution across the example sample'); ylabel('Frequency'); xlabel('Number of connections');