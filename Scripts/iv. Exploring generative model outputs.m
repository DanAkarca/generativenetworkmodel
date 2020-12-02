%% Exploring the generative model outputs
%{
Akarca D, et al. A generative network model of neurodevelopment.
Code written by Dr Danyal Akarca, MRC Cognition and Brain Sciences Unit.
University of Cambridge: https://www.neuroscience.cam.ac.uk/directory/profile.php?da434
Email: danyal.akarca@mrc-cbu.cam.ac.uk
2nd December 2020.
%}
%% About this script
%{
This script analyses the outputs of generative models across each rule. 
It replicates Figure 2, parts of Supplementary Figure 2 and Supplementary Tables 1, 2.

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
% load example generative network model output data, ordered
load('generative10000_energy');      % subject by models by nparams matrix
load('generative10000_parameters');  % subject by models by nparams by parameter matrix
load('generative50000_energy');      % subject by nparams matrix
load('generative50000_parameters');  % subject by nparams by parameter matrix
load('generative50000_ks');          % subject by nparams by ks statistic matrix
%% 1. tabulate descriptives of best networks across rules for 10000 simulations [Supplementary Table 1]
% set model types
modeltype = string({...
    'Spatial',...     % spatial model
    'Neighbors',...   % number of common neighbors
    'Matching',...    % matching index
    'C-Avg',...       % average clustering coeff
    'C-Min',...       % minimum clustering coeff
    'C-Max',...       % maximum clustering coeff
    'C-Diff',...      % difference in clustering coeff
    'C-Prod',...      % product of clustering coeff
    'D-Avg',...       % average degree
    'D-Min',...       % minimum degree
    'D-Max',...       % maximum degree
    'D-Diff',...      % difference in degree
    'D-Prod'})';      % product of degree
% initialise
e  = [];
p1 = [];
p2 = [];
% take top performing network for each subject
for sub = 1:270
    % and for each rule
        for rule = 1:13;
            % top energy
            e(sub,:,rule)= generative10000_energy(sub,rule,1);
            % top eta
            p1(sub,:,rule)= generative10000_parameters(sub,rule,1,1);
            % top gamma
            p2(sub,:,rule)= generative10000_parameters(sub,rule,1,2);
        end   
end
% squeeze matrices
e  = squeeze(e);
p1 = squeeze(p1);
p2 = squeeze(p2);
% tablulate and display
networktable = table(modeltype,mean(e)',std(e)',mean(p1)',std(p1)',mean(p2)',std(p2)','VariableNames',...
    {'Rule','Mean E','Std E','Mean P1','Std P1','Mean P2','Std P2'});
disp(networktable);
%% 2. visualise energy landscapes across rules [Figure 2a-d]
% set generative rule to visualise
model = 3; 
% initialise
landscape = [];
% loop for each and keep
for sub = 1:270;
    se                 = squeeze(generative10000_energy(sub,model,:));
    sp                 = squeeze(generative10000_parameters(sub,model,:,:));
    [u x y]            = unique(sp,'rows');
    se                 = se(x);
    se                 = reshape(se,[100 100]);
    landscape(sub,:,:) = se;
end
% mean across participants
landscape = squeeze(mean(landscape,1));
% 3d landscape
figure; 
surf(landscape);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$'); zlabel('Energy');
title(sprintf('%s',modeltype(model)));
xticks(linspace(0,100,5)); xticklabels(-7:3.5:7);
yticks(linspace(0,100,5)); yticklabels(-7:3.5:7);
zticks([]);
grid off; shading flat;
set(gca,'TickLength',[0 0]);
set(0,'DefaultAxesColor','none');
% 2d landscape
figure;
imagesc(landscape);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
title(sprintf('%s',modeltype(model)));
xticks(linspace(0,100,5)); xticklabels(-7:3.5:7);
yticks(linspace(0,100,5)); yticklabels(-7:3.5:7);
grid off; 
c = colorbar;
caxis([0 1]);
c.Label.String = 'Energy'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';
c.Ticks = [];
set(gca,'TickLength',[0 0]);
%% 3. visualise best network distributions acorss rules [Figure 2e]
% rearrange the order of the model types for visualisation
rearrange = [3 2 4 6 7 8 5 9 11 12 13 10 1];
% Figure 2e
figure;  
set(gcf,'Position',[10 10 1800 650])
h = boxplot(e(:,rearrange),'plotstyle','compact','color',[0.3 0.3 0.3]); 
title('Energy distribution of the best network for each generative rule'); 
ylabel('Energy'); 
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]);
xticklabels(modeltype(rearrange));
a = get(get(gca,'children'),'children'); 
bi = [53:65];
li = [66:79];
% set colours
c = ['r' 'r' 'g' 'g' 'g' 'g' 'g' 'b' 'b' 'b' 'b' 'b' 'k']; 
c = flip(c);
for i = 1:length(c)
    set(a(bi(i)),'Color',c(i));
    set(a(li(i)),'Color',c(i));
end
set(gca,'TickLength',[0 0])
%% 4. tabulate descriptives in the narrow matching window [Supplementary Table 2]
% set how many networks to average across
n  = 1;
% initialise
e  = [];
p1 = [];
p2 = [];
% loop over subjects
for sub = 1:270
    % top n energy
    k = generative50000_energy(sub,1:n)';
    e(sub,:)= k;
    % top n eta
    k = generative50000_parameters(sub,1:n,1)';
    p1(sub,:)= k;
    % top n gamma
    k = generative50000_parameters(sub,1:n,2)';
    p2(sub,:)= k;       
end
% tabulate and display
networktable = table(mean(e(:)),std(e(:)),mean(p1(:)),std(p1(:)),mean(p2(:)),std(p2(:)),'VariableNames',...
    {'Mean E','Std E','Mean P1','Std P1','Mean P2','Std P2'});
disp(networktable);
% show parameter distribution in the space
figure;
set(gcf,'Position',[10 10 1800 650])
h = scatterhist(mean(p1,2),mean(p2,2),'kernel','on','direction','out','Marker','x','MarkerSize',30,'Color','k');
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
set(gca,'TickLength',[0 0]);
%% 5. visualise the narrow matching landscape [Figure 2f, Supplementary Figure 1a-d]
% initialis
elandscape = [];
klandscape = [];
ks         = string({'Degree','Clustering coefficient','Betweeness centrality','Edge length'});
% loop for each and keep
for sub = 1:270;
    se             = generative50000_energy(sub,:);
    sp             = squeeze(generative50000_parameters(sub,:,:));
    sk             = squeeze(generative50000_ks(sub,:,:));
    [u x y]        = unique(sp,'rows');
    se             = se(x);
    se             = reshape(se,[120 417]);
    sk             = sk(x,:);
    for i = 1:4
        ska(:,:,i) = reshape(sk(:,i),[120 417]); 
    end
    elandscape(sub,:,:)   = se;
    klandscape(sub,:,:,:) = ska;
end
% mean across participants
elandscape = squeeze(mean(elandscape,1));
klandscape = squeeze(mean(klandscape,1));
% 3d energy landscape
figure; 
surf(elandscape);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$'); zlabel('E');
title('Narrow matching window');
xticks([1 416]); xticklabels([-3.61 0.35]); 
yticks([1 119]); yticklabels([0.21 0.49]); 
grid off; shading flat;
set(gca,'TickLength',[0 0]);
set(0,'DefaultAxesColor','none')
% 2d energy landscape
figure;
imagesc(elandscape);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
title('Narrow matching window');
xticks([1 416]); xticklabels([-3.61 0.35]); 
yticks([1 119]); yticklabels([0.21 0.49]); 
grid off; 
c = colorbar;
c.Label.String = 'Energy'; 
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
set(gca,'TickLength',[0 0]);
% 2d ks landscape
figure;
for i = 1:4
    subplot(2,2,i);
    imagesc(squeeze(klandscape(:,:,i)));
    xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
    title(ks(i));
    xticks([1 416]); xticklabels([-3.61 0.35]); 
    yticks([1 119]); yticklabels([0.21 0.49]); 
    grid off; 
    c = colorbar;
    caxis([0 0.4]);
    c.Label.String = 'KS'; 
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    set(gca,'TickLength',[0 0]);
end
%% 6. Visualise the parameter trade off [Supplementary Figure 1e]
% set a vector of number of parameters to average over
n_set  = [1 10 100 500];
% initialise
e  = {};
p1 = {};
p2 = {};
% loop across sets of parameters
for n = 1:length(n_set)
    % take number of parameters to average over
    number = n_set(n);
    % loop over subjects
    for sub = 1:270
        % energy
        k = generative50000_energy(sub,1:number)';
        e{n}(sub,:) = k;
        % eta
        k = squeeze(generative50000_parameters(sub,1:number,1));
        p1{n}(sub,:) = k;
        % gamma
        k = squeeze(generative50000_parameters(sub,1:number,2));
        p2{n}(sub,:)  = k;  
    end
end
% initialise
mean_e  = [];
mean_p1 = [];
mean_p2 = [];
% average over parameters
for n = 1:length(n_set)
    % energy
    j = e{n};
    mean_e(:,n) = mean(j,2)';
    % eta
    j = p1{n};
    mean_p1(:,n) = mean(j,2)';
    % gamma
    j = p2{n};
    mean_p2(:,n) = mean(j,2');
end
% calculate correlations
r = []; p = [];
for n = 1:length(n_set)
[r(n) p(n)] = corr(mean_p1(:,n),mean_p2(:,n));
end
% form strings for legend
output = {};
for k = 1:length(n_set)
    output{k} = sprintf('%g: R = %g, P = %g',n_set(k),r(k),p(k));
end
% set colours
cols = lines(length(n_set));
% visualise
figure;
for n = 1:length(n_set)
    scatter(mean_p1(:,n),mean_p2(:,n),40,cols(n,:),'filled'); 
    hold on;
    set(gca,'TickLength',[0 0]);
end
legend(string(output));
