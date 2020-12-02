%% Gene expression PLS analysis
%{
Akarca D, et al. A generative network model of neurodevelopment.
Code written by Dr Danyal Akarca, MRC Cognition and Brain Sciences Unit.
University of Cambridge: https://www.neuroscience.cam.ac.uk/directory/profile.php?da434
Email: danyal.akarca@mrc-cbu.cam.ac.uk
2nd December 2020.
%}
%% About this script
%{
This script replicates PLS analysis using data from the Allen Brain Atlas.
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
% set the current directory to the data <<<<<<<<<<<<< SET
datapath = '/Users/da04/Desktop/Generative Models GitHub/Example data/';
% node-wise parameterised costs over time
load('matching_D_average');
% individual node-wise matching averaged over time
load('matching_Fksum_average');
% RNA gene data
load('RNAgenematrix');
load('RNAgenenames');
%% 1. set up PLS hyperparameters
% set the predictor variable 
X     = RNAgenematrix;
% set number of subjects
nsub = 270;
% set number of nodes
nnode = 34;
% set the number of components
ncomp = 6;
% set the number of permutations per gene
nperm = 1000;
%% 2. PLS analysis for variability in parameterised costs
% initialise observed variables
Xloading      = zeros(nsub,length(X),ncomp);
Xscore        = zeros(nsub,nnode,ncomp);
var           = zeros(nsub,2,ncomp);
% initialise permutated variables
Xscore_perm   = zeros(nsub,nnode,ncomp,nperm);
var_perm      = zeros(nsub,2,ncomp,nperm);
% initialise corrected p values for gene loadings
pcorr         = zeros(nsub,length(X),ncomp);
% display estimated time for completion
disp('Estimated time for completion of approximately 30 minutes (~5s per subject)');
% loop through subjects
for sub = 1:270
    % clock the time it takes per subject
    tic;
    % response variable for the subject is its average parameterised costs
    Y = matching_D_average(sub,:)';
    % take only the left hemisphere only
    Y = Y(35:end);
    % starting for this subject
    disp(sprintf('Subject %g of 270: Running PLSR...',sub));
    % run the observed PLS for a specific subject
    [Xloading(sub,:,:),~,Xscore(sub,:,:),~,~,var(sub,:,:)] = plsregress(X,Y,ncomp);
     disp(sprintf('Subject %g of 270: Permuting %g PLSRs...',sub,nperm));
    % run a permutation test for each subject
    for permutation = 1:nperm
        % randomise the response variable
        Y = Y(randperm(length(Y)));
        % run the permuted PLS
        [Xl,~,Xs,~,~,V] = plsregress(X,Y,ncomp);
        % keep permutations
        Xloading_perm(permutation,:,:) = Xl;
        Xscore_perm(sub,:,:,permutation) = Xs;
        var_perm(sub,:,:,permutation) = V;
    end
    % calculate the p vector for this subject for all components
    for comp = 1:ncomp
        null = squeeze(squeeze(Xloading_perm(:,:,comp,:))); % select the null matrix for this component  
        for gene = 1:length(X);
            r                    = Xloading(sub,gene,comp); % observed 
            z                    = null(:,gene)';           % null
            pcorr(sub,gene,comp) = 1-(sum(z<r)/nperm);      % pcorr
        end
    end
    t = toc;
    disp(sprintf('Subject %g of 270 complete (%g seconds)',sub,t));
end
%% 3. Permuted PLS: parameterised matching
% initialise observed variables
Xloading      = zeros(nsub,length(X),ncomp);
Xscore        = zeros(nsub,nnode,ncomp);
var           = zeros(nsub,2,ncomp);
% initialise permutated variables
Xscore_perm   = zeros(nsub,nnode,ncomp,nperm);
var_perm      = zeros(nsub,2,ncomp,nperm);
% initialise corrected p values for gene loadings
pcorr  = zeros(nsub,length(X),ncomp);
% display estimated time for completion
disp('Estimated time for completion of approximately 30 minutes (~5s per subject)');
% loop through subjects
for sub = 1:270
    % clock the time it takes per subject
    tic;
    % response variable for the subject is its average parameterised matching
    Y = matching_Fksum_average(sub,:)';
    % take only the left hemisphere only
    Y = Y(35:end);
    % starting for this subject
    disp(sprintf('Subject %g of 270: Running PLSR...',sub));
    % run the observed PLS for a specific subject
    [Xloading(sub,:,:),~,Xscore(sub,:,:),~,~,var(sub,:,:)] = plsregress(X,Y,ncomp);
     disp(sprintf('Subject %g of 270: Permuting %g PLSRs...',sub,nperm));
    % run a permutation test for each subject
    for permutation = 1:nperm
        % randomise the response variable
        Y = Y(randperm(length(Y)));
        % run the permuted PLS
        [Xl,~,Xs,~,~,V] = plsregress(X,Y,ncomp);
        % keep permutations
        Xloading_perm(permutation,:,:) = Xl;
        Xscore_perm(sub,:,:,permutation) = Xs;
        var_perm(sub,:,:,permutation) = V;
    end
    % calculate the p vector for this subject for all components
    for comp = 1:ncomp
        null = squeeze(squeeze(Xloading_perm(:,:,comp,:))); % select the null matrix for this component  
        for gene = 1:length(X);
            r                    = Xloading(sub,gene,comp); % observed 
            z                    = null(:,gene)';           % null
            pcorr(sub,gene,comp) = 1-(sum(z<r)/nperm);      % pcorr
        end
    end
    t = toc;
    disp(sprintf('Subject %g of 270 complete (%g seconds)',sub,t));
end