%% Demo eNDM 
%
% The goal of this demo is to apply eNDM to mouse data. 
% Uses analytic (closed form) solutions for models in objective functions. 
% Several consistency checks are provided. 
%
% Written by: Pedro D. Maia
% Created: Feb/6/2020
% Last Modified: Feb/19/2020

%%  0.  Setup   

% Initialization
    clear all; close all; clc; 

% Add folder with raw data to path  
    addpath([pwd '/raw_data_mouse'])

% Add library with eNDM functions to path
    addpath([pwd '/lib_eNDM_analytic'])

% Load dataset of interest
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

disp('Results using analytical (closed form) solutions')

%% 1.  NDM with exhaustive search maximizing correlation 
% 
% * Must specify range for parameter search 
% * Objective function will maximize correlation 

% No need to specify value at IC. 
    x0_value = 1; 
    x0 = x0_value*seed_location;

% Define range of diffusivity parameter (beta)
    beta_range = 1:.1:3;

%----------------------------------------------------------    
% Solve NDM for every value of beta in specified range    
% Choose best value from exhaustive search

% Start Loop 
    for j = 1:length(beta_range)
        
        % Set value of beta
        beta = beta_range(j);
        
        % Solve NDM; 
        [y] = NDM_analytic(x0,time_stamps,C,beta);
        
        % Evaluate corr between prediction and data at every time stamp
        % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,j) = corr(y(:,jj),pathology(:,jj), 'rows','complete');
            end
    end
%----------------------------------------------------------    

      %% 1.1 Optimize beta using exhaustive search on sum of Rvalues 
    
    sum_Rvalues = sum(Rvalues);
    [sum_Rmax,opt_index] = find(sum(Rvalues) == max(sum_Rvalues));
    best_Rvalues = Rvalues(:,opt_index);
    sum_Rmax = sum_Rvalues(opt_index);            
    opt_beta = beta_range(opt_index);

% Calculate Square error
    [y] = NDM_analytic(x0,time_stamps,C,opt_beta);
    SE = nansum(nansum((y - pathology).^2));

% Display results 
    disp('--------------------------------------------------')
    disp(['NDM maximizing sum of R values over all ' num2str(length(time_stamps)) ' time stamps']);   
    disp(' ')
    disp('Optimal beta parameter value')
    disp(opt_beta)
    disp(['R values at each time stamp'])
    disp(best_Rvalues);
    disp('Sum of R values')
    disp(sum_Rmax)
    disp('Square error') 
    disp(SE)
    
      %% 1.2 Optimize beta using exhaustive search on last-time Rvalue alone 

    last_Rvalues = Rvalues(end,:);
    [lastRmax,opt_index] = find(last_Rvalues == max(last_Rvalues));
    best_Rvalues = Rvalues(1:3,opt_index);
    lastRmax = Rvalues(3,opt_index);
    opt_beta = beta_range(opt_index);
    
% Calculate Square error
    [y] = NDM_analytic(x0,time_stamps,C,opt_beta);
    SE = nansum(nansum((y - pathology).^2));    

    % Display results 
    disp('--------------------------------------------------')
    disp(['NDM maximizing R value at time stamp ' num2str(length(time_stamps))]);           
    disp(' ')
    disp('Optimal beta parameter value')
    disp(opt_beta)
    disp(['R values at each time stamp'])
    disp(best_Rvalues);
    disp(['R value at time stamp ' num2str(length(time_stamps))])
    disp(lastRmax)
    disp('Square error') 
    disp(SE)

    
%% 2.  NDM with fmincon minimizing MSE 

% Parameters to fit
    % param(1) = beta
    % param(2) = x0_value

% Extra input to required by objective function
    % seed_location 
    % pathology
    % time_stamps
    % C
    
% Set initial guess for beta
    init_guess_params(1) = 1;

% Set initial guess for x0_value
    init_guess_params(2) = nansum(pathology(:,1));

% Set lower bounds (lb) for parameters
    lb = [0,0];

% Set upper bounds (ub) for parameters
%    ub = [10,nansum(pathology(:,end))];
ub = [3,3]
clearvars  x0_value beta Rvalues

     %% 2.1 Optimize parameters minimizing square errror on all time stamps 

% Apply fmincon
    [param, fval] = fmincon(@(param)objfun_NDM_analytic(param,seed_location,pathology,time_stamps,C),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDM with the optimal parameters
    beta = param(1);
    x0 = seed_location*param(2);
    y = NDM_analytic(x0,time_stamps,C,beta);

 % Evaluate corr between prediction and data at every time stamp
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(y(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDM minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['Beta = ' num2str(param(1)) ', x0 value = ' num2str(param(2))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(y,pathology,time_stamps)
clearvars y 

return



















     %% 2.2 Optimize parameters minimizing square errror on last time stamps alone 

% Apply fmincon
    [param, fval] = fmincon(@(param)objfun_NDM_analytic(param,seed_location,pathology(:,end),time_stamps(end),C),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDM with the optimal parameters
    beta = param(1);
    x0 = seed_location*param(2);
    y = NDM_analytic(x0,time_stamps,C,beta);

 % Evaluate corr between prediction and data at every time stamp
 % Save Rvalues in a matrix 
                 Rvalues(end,:) = corr(y(:,end),pathology(:,end), 'rows','complete');
        
 % Display results
    disp('--------------------------------------------------')
    disp('NDM minimizing quadratic error at last time stamp only with fmincon');   
    disp(' ')
    disp(['Beta = ' num2str(param(1)) ', x0 value = ' num2str(param(2))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval)
 
 % Plot prediction vs data using optimal parameters
 plot_pred_vs_data(y(:,end),pathology(:,end),time_stamps(end))
 
 
 %% 3. NDM with source (NDMwS) minimizing MSE
 % Parameters to fit
    % param(1) = beta
    % param(2) = x0_value
    % param(3) = alpha1 = linear growth/clearance term in x

% Extra input to required by objective function
    % seed_location 
    % pathology
    % time_stamps
    % C
    
% Set initial guess for beta
    init_guess_params(1) = 1;

% Set initial guess for x0_value
    init_guess_params(2) = nansum(pathology(:,1));

% Set initial guess for alpha1 and alpha2
    init_guess_params(3) = .5;
    
    
% Set lower bounds (lb) for parameters
    lb = [0,0,-5];

% Set upper bounds (ub) for parameters
    ub = [10,nansum(pathology(:,end)),5];

%% 3.1 Optimize parameters minimizing square errror on all time stamps 

% Apply fmincon
    [param, fval] = fmincon(@(param)objfun_NDMwS_analytic(param,seed_location,pathology,time_stamps,C),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDMwS with the optimal parameters
    beta = param(1);
    x0 = seed_location*param(2);
    alpha1 = param(3);

    y = NDMwS_analytic(x0,time_stamps,C,beta,alpha1);
     
 % Evaluate corr between prediction and data at every time stamp
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(y(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDMwS minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['Beta = ' num2str(param(1)) ', x0 value = ' num2str(param(2)) ...
            ' ,alpha1 = ' num2str(param(3))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(y,pathology,time_stamps)
clearvars y


%% 4. NDM with source and with microglia (NDMwSwM) minimizing MSE
 
% Load Microglia and do min-max normalization 
load microglia_density.mat
u = microdata; 
u = (u - min(u))./(max(u) - min(u));
%u = (u - mean(u))/std(u);
%u = 1+u;

% Parameters to fit
    % param(1) = beta
    % param(2) = x0_value
    % param(3) = alpha1  = linear growth/clearance rate term in x
    % param(4) = alpha_2 = microglia reweighting for alpha_1  

% Extra input to required by objective function
    % seed_location 
    % pathology
    % time_stamps
    % C
    
% Set initial guess for beta
    init_guess_params(1) = 1.8;

% Set initial guess for x0_value
    %init_guess_params(2) = nansum(pathology(:,1));
    init_guess_params(2) = 0.35;

    % Set initial guess for alpha1 and alpha2
    init_guess_params(3) = .5;
    init_guess_params(4) = 0;
    
% Set lower bounds (lb) for parameters
    lb = [0,0,-5,-5];

% Set upper bounds (ub) for parameters
    ub = [10,nansum(pathology(:,end)),5,5];

%% 4.1 Optimize parameters minimizing square errror on all time stamps 

% Apply fmincon
    [param, fval] = fmincon(@(param)objfun_NDMwSwM_analytic(param,seed_location,pathology,time_stamps,C,u),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDMwS with the optimal parameters
    beta = param(1);
    x0 = seed_location*param(2);
    alpha1 = param(3);
    alpha2 = param(4);

    y = NDMwSwM_analytic(x0,time_stamps,C,u,beta,alpha1,alpha2);
     
 % Evaluate corr between prediction and data at every time stamp
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(y(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDMwSwM minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['Beta = ' num2str(param(1)) ', x0 value = ' num2str(param(2)) ...
            ' ,alpha1 = ' num2str(param(3)),  ' ,alpha2 = ' num2str(param(4))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(y,pathology,time_stamps)
clearvars y

                   
