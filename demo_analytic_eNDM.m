%% Demo eNDM using analytic (closed form) solutions
%
% The goal of this demo is to apply eNDM to mouse data. 
% Uses analytic (closed form) solutions for models in objective functions. 
% Several consistency checks are provided. 
%
% Written by: Pedro D. Maia
% Created: Feb/6/2020
% Last Modified: 

%%  0.  Setup   

% Initialization
    clear all; close all; clc; 

% Add folder with raw data to path  
    addpath([pwd '/raw_data_mouse'])

% Add library with eNDM analytic functions to path
    addpath([pwd '/lib_eNDM_analytic'])

% Load dataset of interest
   load celltypedata_mrx3_540.mat
   load MouseDataForPedro.mat

% Select connectome C (426x426)
    C = Networks.nd;        % symmmetric  
    %C = Networks.ret;      % non-symmetric  
    %C = Networks.ant;     % non-symmetric      

% Normalize C
    cmax = max(max(C));
    cmin = min(min(C)); 
    C = (C - cmin)./(cmax-cmin);
        
% Define number of regions of interest (nroi) 
    nroi = size(C,1);

% Load time stamps, pathology measurements, and the seed_location 
   
        % Hippocampus Injection   
           time_stamps = [1,3,6];
           pathology = data426.IbaHippInj;
           seed_location = seed426.IbaHippInj;  

        % Striatal Injection 
%            time_stamps = [1,3,6];
%            pathology = data426.IbaStrInj;
%            seed_location = seed426.IbaStrInj;
         


%% 4. eNDM 

% Select cell type from list 
for j = 1:size(celltypedata,2)

disp([' Results for Cell Type #' num2str(j)])
u = celltypedata(:,j); 
u = (u - min(u))./(max(u) - min(u));

% Parameters to fit
    % param(1) = x0_value
    % param(2) = beta1
    % param(3) = alpha1  = linear growth/clearance rate term in x
    % param(4) = alpha_2 = microglia reweighting for alpha_1  
    % param(5) = beta2 = microglia reweighting for beta1

% Extra input to required by objective function
    % seed_location 
    % pathology
    % time_stamps
    % C
    
% Set initial guess for [x0_value, beta1, alpha1, alpha2, beta2]
    init_guess_params(1) = nansum(pathology(:,1));
    init_guess_params(2) = 1.8;
    init_guess_params(3) = .5;
    init_guess_params(4) = 0;
    init_guess_params(5)= 0;    
        
% Set lower bounds (lb) for [x0_value, beta1, alpha1, alpha2, beta2]
    lb = [0,0,0,-5,-5];

% Set upper bounds (ub) for [x0_value, beta1, alpha1, alpha2, beta2]
    ub = [.4*nansum(pathology(:,1)),10,5,5,5];

%% 4.1 Optimize parameters minimizing square errror on all time stamps 

% Apply fmincon
    [param, fval] = fmincon(@(param)objfun_eNDM_analytic(param,seed_location,pathology,time_stamps,C,u),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDMwS with the optimal parameters
    x0 = seed_location*param(1);
    beta1 = param(2);
    alpha1 = param(3);
    alpha2 = param(4);
    beta2= param(5);

    y = eNDM_analytic(x0,time_stamps,C,u,beta1,alpha1,alpha2,beta2);
     
 % Evaluate corr between prediction and data at every time stamp
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(y(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('eNDM minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['x0 value = ' num2str(param(1)) ', beta1 = ' num2str(param(2))...
            ', alpha1 = ' num2str(param(3)),  ', alpha2 = ' num2str(param(4)) ', beta2 = ' num2str(param(5))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(y,pathology,time_stamps)
clearvars y

                   
end