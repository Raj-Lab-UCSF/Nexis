%% Demo [Analytics vs Numeric]
%
% *****NB: Only the numeric function called is modified to include multiple costfun
% options while we investigate the proper naming and reorganization of the
% analytic functions*****

% The goal of this demo is to apply eNDM to mouse data. 
% Uses numerical solutions for models using ode solvers in objective functions. 
% Several consistency checks are provided. 
%
% Written by: Pedro D. Maia
%
% Created: Feb/6/2020
%
% Last Modified: Feb/19/2020

%%  0.  Setup   

% Initialization
    clear all; close all; clc; 

% Add folder with raw data to path  
    addpath([pwd '/raw_data_mouse'])

% Add library with eNDM functions to path
    addpath([pwd '/lib_eNDM_numeric'])
    addpath([pwd '/lib_eNDM_analytic'])

% Load dataset of interest
   load Tasic_CTMaps.mat %switch out for other cell type data
   load eNDM_mousedata.mat
   
% Specify cost function, options: 'sse_sum', 'sse_end', 'rval_sum',
% 'rval_end', 'LinR'
    costfun = 'sse_sum';
    
% Defining LinR and type of correlation to display, R_c or R
    LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);
    corrtype = 'R_c';

% Select connectome C (426x426), use only ND (symmetric) for now
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
   
        % Mouse pathology data inputs based on datsetname     
           datsetname = 'IbaHippInj';
   
        % Mouse pathology data inputs based on datsetname     
           time_stamps = tpts.(datsetname);
           pathology = data426.(datsetname);
           
           pathology = pathology./norm(pathology,2);
           seed_location = seed426.(datsetname);  
           
           base_location = base426.(datsetname);
           
%% 1.  [Analytics vs Numeric] NDM, fmincon minimizing MSE  
    
% Set initial guess for beta
    init_guess_params(1) = 1;

% Set initial guess for x0_value
    init_guess_params(2) = nansum(pathology(:,1));

% Set lower bounds (lb) for parameters
    lb = [0,0];

% Set upper bounds (ub) for parameters
%    ub = [10,nansum(pathology(:,1))];
ub = [3,3];
% Apply fmincon w/ numeric 
%     [param_num, fval_num] = fmincon(@(param)objfun_NDM_numeric(param,seed_location,pathology,time_stamps,C),...
%                            init_guess_params,[],[],[],[],lb,ub,[]);
[param_num, fval_num] = fmincon(@(param)objfun_NDM_numeric_costopts(param,seed_location,pathology,time_stamps,C,costfun),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDM with the optimal parameters
    beta_num = param_num(1);
    x0_num = seed_location*param_num(2);
    ynum = NDM_numeric(x0_num,time_stamps,C,beta_num);

 %% Display Numeric Results
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 if strcmp(corrtype,'R')
                    Rvalues(jj,:) = corr(ynum(:,jj),pathology(:,jj), 'rows','complete');
                elseif strcmp(corrtype,'R_c')
                    naninds = isnan(pathology(:,1));
                    newxt = ynum; newxt(naninds,:) = [];
                    newpath = pathology; newpath(naninds,:) = [];
                    Rvalues(jj,:) = LinRcalc(newxt(:,jj),newpath(:,jj));
                end
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDM minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['NUMERIC Beta = ' num2str(param_num(1)) ', x0 value = ' num2str(param_num(2))])
    disp(' ')
    disp([corrtype ' values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp(costfun) 
    disp(fval_num)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(ynum,pathology,time_stamps)

clearvars y 


