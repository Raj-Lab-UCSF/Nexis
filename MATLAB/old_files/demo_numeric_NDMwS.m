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

    datsetname = 'IbaHippInj';
   
        % Mouse pathology data inputs based on datsetname  
           time_stamps = tpts.(datsetname);
           pathology = data426.(datsetname);
           
           pathology = pathology./norm(pathology,2);
           seed_location = seed426.(datsetname);  
           
           base_location = base426.(datsetname);

%% 2.  [Analytics vs Numeric] NDMwS, fmincon minimizing MSE 
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
%    ub = [10,nansum(pathology(:,end)),5];
    ub = [3,3,3];

% Apply fmincon w/ numeric 
    [param_num, fval_num] = fmincon(@(param)objfun_NDMwS_numeric_costopts(param,seed_location,pathology,time_stamps,C,costfun),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDMwS with the optimal parameters
    beta_num = param_num(1);
    x0_num = seed_location*param_num(2);
    alpha1_num = param_num(3);
    ynum = NDMwS_numeric(x0_num,time_stamps,C,beta_num,alpha1_num);
    
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
    disp('NDMwS minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['NUMERIC Beta = ' num2str(param_num(1)) ', x0 value = ' num2str(param_num(2)) ...
        ', alpha1 = ' num2str(param_num(3))])
    disp(' ')
    disp([corrtype ' values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp(costfun) 
    disp(fval_num)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(ynum,pathology,time_stamps)
clearvars y 
