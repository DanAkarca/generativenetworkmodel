%% Running the homophily generative model
%{
Akarca D, et al. A generative network model of neurodevelopment.
Code written by Dr Danyal Akarca, MRC Cognition and Brain Sciences Unit.
University of Cambridge: https://www.neuroscience.cam.ac.uk/directory/profile.php?da434
Email: danyal.akarca@mrc-cbu.cam.ac.uk
2nd December 2020.
%}
%% About this script
%{
This script runs the homophily "matching" generative network model over a
user-defined parameter combinatoin. This script, beyond a simple grid search 
and parameter optimisation, keeps edge-wise and node-wise data 
for each developmental time step, including the network (A), values (K), 
parameterised matching (Fk), and wiring probability (P) over the simulation.

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
% load example euclidean distances
load('example_euclidean.mat');
%% 1. set the single parameter combination to run
% note, in our work this was previously tailored to each subject's optimally performing
% generative model. Here, we can select an arbitrary combination across the
% whole sample.
eta    = -1.5;
gam    = 0.25;
params = [eta gam];
%% 2. set hyperparameters and run generative simulations with matching trajectories
% nsub
nsub = 270;
% initialise 
matching_K     = {};   % edge-wise K
matching_Ksum  = {};   % node-wise K
matching_Fk    = {};   % edge-wise parameterised K
matching_Fksum = {};   % node-wise parameterised K
matching_P     = {};   % edge-wise wiring probability
matching_Psum  = {};   % node-wise wiring probability
matching_A = {};       % networks
% run generative models for each subject
for sub = 1:nsub;
    tic  
    % set target connectomes 
    Atgt_set = example_binarised_connectomes;
    % form seed network 
    proportion  = 0.2; % set proportion of edges in common across subjects [0 1] 
    connections = squeeze(mean(Atgt_set,1));
    index       = find(connections==proportion);
    A           = zeros(size(connections));
    A(index)    = 1; A = upper(A);
    % set costs
    D           = euclidean;
    % set target
    Atgt        = squeeze(Atgt_set(sub,:,:));
    % params were defined earlier
    params      = params;
    % number of bi-directional connections
    m           = nnz(Atgt)/2;
    % model var
    modelvar    = [{'powerlaw'},{'powerlaw'}];
    % minimum edge
    epsilon     = 1e-5;
    % run initial matching
    n           = length(D);
    nparams     = size(params,1);
    b           = zeros(m,nparams);
    K           = matching_ind(A);
    K           = K + K';
    % keep the K,Fk,P,A at each iteration, under the current parameter combination
    Kall        = [];
    Fkall       = [];
    Pall        = [];
    Aall        = [];
    % save the first K
    Kall(1,:,:)  = K;
    % save the first A
    Aall(1,:,:)  = A;
    for iparam = 1:nparams
        eta = params(iparam,1);
        gam = params(iparam,2);
        K = K + epsilon;
        n = length(D);
        mseed = nnz(A)/2;
        mv1 = modelvar{1};
        mv2 = modelvar{2};
        switch mv1
            case 'powerlaw'
                Fd = D.^eta;
            case 'exponential'
                Fd = exp(eta*D);
        end
        switch mv2
            case 'powerlaw'
                Fk = K.^gam;
            case 'exponential'
                Fk = exp(gam*K);
        end
        Ff = Fd.*Fk.*~A;
        [u,v] = find(triu(ones(n),1));
        indx = (v - 1)*n + u;
        P = Ff(indx);
        % save the first parameterised K in Fk
        FKall(1,:,:)  = Fk;
        % save the first probabilities
        Ff(isinf(Ff)) = 0;
        Pall(1,:,:)   = Ff;
        step = 2; 
        for ii = (mseed + 1):m
            C = [0; cumsum(P)];
            r = sum(rand*C(end) >= C);
            uu = u(r);
            vv = v(r);
            A(uu,vv) = 1;
            A(vv,uu) = 1;
            updateuu = find(A*A(:,uu));
            updateuu(updateuu == uu) = [];
            updateuu(updateuu == vv) = [];
            updatevv = find(A*A(:,vv));
            updatevv(updatevv == uu) = [];
            updatevv(updatevv == vv) = [];
            c1 = [A(:,uu)', A(uu,:)];
            for i = 1:length(updateuu)
                j = updateuu(i);
                c2 = [A(:,j)' A(j,:)];
                use = ~(~c1&~c2);
                use(uu) = 0;  use(uu+n) = 0;
                use(j) = 0;  use(j+n) = 0;
                ncon = sum(c1(use))+sum(c2(use));
                if (ncon==0)
                    K(uu,j) = epsilon;
                    K(j,uu) = epsilon;
                else
                    K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                    K(j,uu) = K(uu,j);
                end
            end
            c1 = [A(:,vv)', A(vv,:)];
            for i = 1:length(updatevv)
                j = updatevv(i);
                c2 = [A(:,j)' A(j,:)];
                use = ~(~c1&~c2);
                use(vv) = 0;  use(vv+n) = 0;
                use(j) = 0;  use(j+n) = 0;
                ncon = sum(c1(use))+sum(c2(use));
                if (ncon==0)
                    K(vv,j) = epsilon;
                    K(j,vv) = epsilon;
                else
                    K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                    K(j,vv) = K(vv,j);
                end
            end
            switch mv2
                case 'powerlaw'
                    Fk = K.^gam;
                case 'exponential'
                    Fk = exp(gam*K);
            end
            Kall(step,:,:)  = K;                   
            Fkall(step,:,:) = Fk;                  
            Aall(step,:,:)  = A;                     
            Ff = Fd.*Fk.*~A;
            P = Ff(indx);
            Ff(isinf(Ff))   = 0; 
            Pall(step,:,:)  = Ff;                
            % change the step
            step = step+1;
        end
        b(:,iparam)         = find(triu(A,1));
        % keep data
        matching_K{sub}     = Kall;                  
        matching_Fk{sub}    = Fkall;                 
        matching_Ksum{sub}  = squeeze(sum(Kall,2)); 
        matching_Fksum{sub} = squeeze(sum(Fkall,2)); 
        matching_P{sub}     = Pall;                 
        matching_Psum{sub}  = squeeze(sum(Pall,2)); 
        matching_A{sub}     = Aall;                 
    end
    time = toc;
    disp(sprintf('Generative network model computed: Subject %g took %g seconds',sub,time));
end
%% 3. explore trajectory outputs
% below are a range of ways to explore the derivatives of the generative
% model outputs. Some are taken forward to the next step of analyses.
% 1. network over time
matching_A;
% 2. edge-wise matching over time
matching_K;
% 3. node-wise matching over time
matching_Ksum;
% 4. calculate grouped node-wise matching over time
matching_Ksum_group   = [];
    % find average developmental time
    h = [];
    for i = 1:length(matching_Ksum)
        h(i) = size(matching_Ksum{i},1);
    end
    av = floor(mean(h));
    % scale all the data to this into one matrix
    for sub = 1:nsub;
        matching_Ksum_group(sub,:,:) = imresize(matching_Ksum{sub},[av 68]);
    end
% 5. calculate individual node-wise matching averaged over time
matching_Ksum_average = [];
for sub = 1:nsub;
    matching_Ksum_average(sub,:) = mean(matching_Ksum{sub});
end
% 6. calculate grouped node-wise matching averaged over time
matching_Ksum_group_average = squeeze(mean(matching_Ksum_group));
% 7. edge-wise parameterised matching over time
matching_Fk;
% 8. node-wise parameterised matching over time
matching_Fksum;
% 9. calculate grouped node-wise parameterised matching over time
matching_Fksum_group   = [];
    % find average developmental time
    h = [];
    for i = 1:length(matching_Fksum)
        h(i) = size(matching_Fksum{i},1);
    end
    av = floor(mean(h));
    % scale all the data to this into one matrix
    for sub = 1:nsub;
        matching_Fksum_group(sub,:,:) = imresize(matching_Fksum{sub},[av 68]);
    end
% 10. calculate individual node-wise parameterised matching averaged over time
matching_Fksum_average = [];
for sub = 1:nsub;
    matching_Fksum_average(sub,:) = mean(matching_Fksum{sub});
end
% 11. calculate grouped node-wise parameterised matching averaged over time
matching_Fksum_group_average = squeeze(mean(matching_Fksum_group));
% 12. edge-wise matching wiring probabilities over time
matching_P;
% 13. node-wise matching wiring probabilities over time
matching_Psum;
% 14. calculate grouped node-wise matching wiring probabilities over time
matching_Psum_group   = [];
    % find average developmental time
    h = [];
    for i = 1:length(matching_Psum)
        h(i) = size(matching_Psum{i},1);
    end
    av = floor(mean(h));
    % scale all the data to this into one matrix
    for sub = 1:nsub;
        matching_Psum_group(sub,:,:) = imresize(matching_Psum{sub},[av 68]);
    end
% 15. calculate individual node-wise matching wiring probabilities averaged over time
matching_Psum_average = [];
for sub = 1:nsub;
    matching_Psum_average(sub,:) = mean(matching_Psum{sub});
end
% 16. calculate grouped node-wise matching wiring probabilities averaged over time
matching_Psum_group_average = squeeze(mean(matching_Psum_group));