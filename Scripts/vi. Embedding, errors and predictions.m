%% Embedding, errors and predictions.
%{
Akarca D, et al. A generative network model of neurodevelopment.
Code written by Dr Danyal Akarca, MRC Cognition and Brain Sciences Unit.
University of Cambridge: https://www.neuroscience.cam.ac.uk/directory/profile.php?da434
Email: danyal.akarca@mrc-cbu.cam.ac.uk
2nd December 2020.
%}
%% About this script
%{
This script replicates analysis concerning predictions of measures from inside and outside of the energy equation.
This replicates the code (not findings - as different data) for Figures 3 and Supplementary Figure 3.

Due to NHS data restrictions, we have not shared observed connectome data.
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
% load rewired observed and simulated connectomes
load('example_binarised_connectomes');
load('example_optimised_simulations');
load('example_euclidean');
load('generative50000_parameters');
%% 1. form observed and simulated measures
% set number of subjects
nsub = 270;
% initialise
observed_local_statistics   = [];
simulated_local_statistics  = [];
observed_global_statistics  = [];
simulated_global_statistics = [];
% loop over subjects
for sub = 1:nsub;
        % take observed network
    observed_network  = ...
        squeeze(example_binarised_connectomes(sub,:,:));
    % take simulated network
    simulated_network = ...
        squeeze(example_optimised_simulations(sub,:,:));
    % degree
    observed_local_statistics(sub,:,1)     = ...
        degrees_und(observed_network);
    simulated_local_statistics(sub,:,1)    = ...
        degrees_und(simulated_network);
    % betweeness
    observed_local_statistics(sub,:,2)     = ...
        betweenness_bin(observed_network);
    simulated_local_statistics(sub,:,2)    = ...
        betweenness_bin(simulated_network);
    % clustering
    observed_local_statistics(sub,:,3)     = ...
        clustering_coef_bu(observed_network);
    simulated_local_statistics(sub,:,3)    = ...
        clustering_coef_bu(simulated_network);
    % edge length
    j = triu(observed_network,1)>0;
    r = triu(simulated_network,1)>0;
    for node = 1:68
        h = j(node,:); % edges from node
        s = r(node,:);
        observed_local_statistics(sub,node,4)  = ...
            sum(euclidean(node,h)); % sum total distance from the nodes
        simulated_local_statistics(sub,node,4) = ...
            sum(euclidean(node,s));
    end
    % local efficiency
    observed_local_statistics(sub,:,5)     = ...
        efficiency_bin(observed_network,1);
    simulated_local_statistics(sub,:,5)    = ...
        efficiency_bin(simulated_network,1);
    % rich club
    a{sub}                                   = ...
        rich_club_bu(observed_network);
    b{sub}                                   = ...
        rich_club_bu(simulated_network);
    % assortativity
    observed_global_statistics(sub,2)      = ...
        assortativity_bin(observed_network,0);
    simulated_global_statistics(sub,2)     = ...
        assortativity_bin(simulated_network,0);
    % transitivity
    observed_global_statistics(sub,3)      = ...
        transitivity_bu(observed_network);
    simulated_global_statistics(sub,3)     = ...
        transitivity_bu(simulated_network);
    % maximum modularity statistic
    [~,observed_global_statistics(sub,4)]  = ...
        modularity_und(observed_network);
    [~,simulated_global_statistics(sub,4)] = ...
        modularity_und(simulated_network);
end
% number of rich club nodes
for sub = 1:nsub;
    [~,observed_global_statistics(sub,1)]  = ...
        max(a{sub});
    [~,simulated_global_statistics(sub,1)] = ...
        max(b{sub});
end
%% 2. plot cumulative distribution functions and compute averaged correlations for local measures [Figure 3 and Supplementary Figure 2a]
% cumulative distribution functions
rng('default') 
% local measures
localmeasures = string(...
    {'degree','betweeness','clustering','edge length','local efficiency'});
% select local measure to plot
measure =  1;
% take measure
observed_measures  = squeeze(observed_local_statistics(:,:,measure));
simulated_measures = squeeze(simulated_local_statistics(:,:,measure));
% visualise cdf
figure;
histogram(observed_measures(:),...
    'normalization','cdf',...
    'displaystyle','stairs',...
    'edgecolor','r',...
    'edgealpha',0.5,...
    'linewidth',3);
hold on;
histogram(simulated_measures(:),...
    'normalization','cdf',...
    'displaystyle','stairs',...
    'edgecolor','b',...
    'edgealpha',0.5,...
    'linewidth',3);
xlim([0 max([observed_measures(:);simulated_measures(:)])]);
xlabel(localmeasures(measure));
ylabel(sprintf('F(%s)',localmeasures(measure)));
set(gca,'TickLength',[0 0]);
title('');
set(gca,'color','none');
legend(string({'Observed','Simulated'}),'location','northeastoutside');
% compute correlations
mean_observed  = mean(observed_measures);
mean_simulated = mean(simulated_measures);
[r p]          = corr(mean_observed',mean_simulated');
% compute line of best fit
a = polyfit(mean_observed,mean_simulated,1); 
f = polyval(a,mean_observed);
% visualise
figure;
scatter(mean_observed,mean_simulated,40,'filled','k');
hold on;
a = plot(mean_observed,f,'k'); a.LineWidth = 2;
title(sprintf('R = %g. P = %g',r,p));
xlabel(sprintf('Observed %s',localmeasures(measure)));
ylabel(sprintf('Simulated %s',localmeasures(measure)));
set(gca,'TickLength',[0 0]);
%% 3. Form an error matrix, across measures, for later visualisation [Supplementary Figure 3a-d]
% Compute an error, which shows how much the simulated network exagerates measure compared to the observed
ErrorData = [];
for measure =  1:4;
    % take measure
    mean_observed        = mean(squeeze(observed_local_statistics(:,:,measure)));
    mean_simulated       = mean(squeeze(simulated_local_statistics(:,:,measure)));
    ErrorData(:,measure) = (mean_simulated - mean_observed)';
end
%% 4. Rank nodes by cumulative error rank [Supplementary Figure 3e]
absError  = abs(ErrorData);
absError  = normalize(absError,'range');
absError  = mean(absError')';
[~,~,rnk] = unique(absError);
% visualise
figure; 
scatter(rnk,absError,40,'filled'); xlabel('Error rank'); ylabel('Absolute generative error');
%% 5. plot cumulative distribution functions and compute averaged correlations for global measures [Supplementary Figure 2b-e]
% cumulative distribution functions
rng('default') 
% local measures
globalmeasures = string(...
    {'rich clubs','assortativity','transitivity','modularity'});
% select global measure to plot *** 
measure = 1;
% take measure
observed_measures  = squeeze(observed_global_statistics(:,measure));
simulated_measures = squeeze(simulated_global_statistics(:,measure));
% compute cdf
observed  = cdfplot(observed_measures(:)); ...
    xo = observed.XData; yo = observed.YData;
simulated = cdfplot(simulated_measures(:)); ...
    xs = simulated.XData; ys = simulated.YData;
% visualise cdf
figure;
plot(xo,yo,'r');
hold on;
plot(xs,ys,'b');
xlabel(localmeasures(measure));
ylabel(sprintf('F(%s)',localmeasures(measure)));
set(gca,'TickLength',[0 0]);
title('');
set(gca,'color','none');
legend(string({'Observed','Simulated'}));
% compute correlations
[r p]          = corr(observed_measures,simulated_measures);
% compute line of best fit
a = polyfit(observed_measures,simulated_measures,1); 
f = polyval(a,observed_measures);
% take data
data = [observed_measures simulated_measures];
% change size based on repititions
[xyUnique, ignore, ixs] = unique(data,'rows');
counts = zeros(size(xyUnique,1),1);
for ix = 1:size(counts,1);
    counts(ix) = sum(ixs == ix);
end
% visualise 
figure;
scatter(...
    xyUnique(:,1), ... 
    xyUnique(:,2), ... 
    counts*100, ... 
    'k', ...            
    'filled');       
hold on;
a = plot(observed_measures,f,'k'); a.LineWidth = 6;
title(sprintf('R = %g. P = %g',r,p));
xlabel(sprintf('Observed %s',globalmeasures(measure)));
ylabel(sprintf('Simulated %s',globalmeasures(measure)));
set(gca,'TickLength',[0 0]);
%% 6. Present parameter-based predictions [part of Figure 4 and Supplementary Table 3]
% nsub
nsub = 270;
% mean across top n simulations
n = 500;
% initialise
params_set = [];
% calculate optimal mean n parameters
for sub = 1:nsub; % nsub of study
    p                 = squeeze(generative50000_parameters(sub,1:n,:));
    p                 = mean(p,1); % mean across parameters
    params_set(sub,:) = p;                              
end
% correlations with eta and gamma with local and global connectome measures
[r p] = corr(params_set,[squeeze(mean(observed_local_statistics,2)) observed_global_statistics]);
% visualise
figure;
set(gcf, 'Position',  [100, 100, 1400, 800]);
subplot(2,1,1); imagesc(r); caxis([-1 1]); colorbar;
sgtitle(sprintf('%g top parameter combinations',n));
title('Correlation coefficient');
yticks([1:2]); yticklabels({'$$ \eta $$','$$ \gamma $$'});
set(gca,'TickLength',[0 0]);
subplot(2,1,2); imagesc(p); caxis([0 1]); colorbar;
title('Correlation significance');
yticks([1:2]); yticklabels({'$$ \eta $$','$$ \gamma $$'});
set(gca,'TickLength',[0 0]);
% note, in out work we also compared cortical morphology, which is unable
% to be shared.