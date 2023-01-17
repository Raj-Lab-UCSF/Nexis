%%  0.  Setup   

% Initialization
    clear; close all; clc; 

% Add folder with raw data to path  
    addpath([pwd '/raw_data_mouse'])

% Add library with eNDM functions to path
    addpath([pwd '/lib_eNDM_numeric'])
    addpath([pwd '/lib_eNDM_analytic'])

% Load dataset of interest
    study = 'Tasic';
    
    if strcmp(study,'Tasic')
        load('Tasic_CTMaps.mat');
        load('classkey_tasic.mat','classkey_tasic');
        classkey = classkey_tasic;
    else
        load('Zeisel_CTMaps.mat');
        load('classkey_zeisel.mat','classkey_zeisel');
        classkey = classkey_zeisel;
    end
    load eNDM_mousedata.mat

% Specify cost function, options: 'sse_sum', 'sse_end', 'rval_sum',
% 'rval_end', 'LinR'
    costfun = 'LinR';
    
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

% Define cell type indices (from classkey)
    type_inds = [9,10,24];

% Define cell type matrix, U
    U = celltypedata.densities(:,type_inds);
    U = U./nansum(U);

% Analytic vs. numeric ODE implementation
    solvetype = 'analytic';

% Load time stamps, pathology measurements, and the seed_location 

    datsetname = 'IbaHippInj';

    % Mouse pathology data inputs based on datsetname  
       time_stamps = tpts.(datsetname);
       pathology = data426.(datsetname);
%        pathology = pathology./norm(pathology,2);
       pathology = pathology./norm(pathology(~isnan(pathology)),2);
       seed_location = seed426.(datsetname);  

       base_location = base426.(datsetname);

%% B.  Generalized eNDM

% Parameters to fit
    % param(1) = seed rescale factor
    % param(2) = alpha
    % param(3) = beta
    % param(4:(n_types+3)) = a
    % param((n_types+4):(2*n_types+3)) = b
    % param((2*n_types+4):(3*n_types+3)) = p

% Extra input to required by objective function
    % seed_location 
    % pathology
    % time_stamps
    % C
    % solvetype

% Set initial guesses and bounds for parameters
    n_types = length(type_inds);
    init_guess_params = zeros((3*n_types+3),1);
    ub = init_guess_params; lb = ub;

    % seed rescale value
        init_guess_params(1) = nansum(pathology(:,1));
        lb(1) = 0;
        ub(1) = 3;
        
    % alpha
        init_guess_params(2) = 0.5;
        lb(2) = -3;
        ub(2) = 3;
        
    % beta
        init_guess_params(3) = 1;
        lb(3) = 0; % beta must be strictly non-negative
        ub(3) = 3;
    
    % a
        init_guess_params(4:(n_types+3)) = 0; % must be set to 0 for solvetype = 'analytic'
        lb(4:(n_types+3)) = 0; % must be set to 0 for solvetype = 'analytic'
        ub(4:(n_types+3)) = 0; % must be set to 0 for solvetype = 'analytic'

    % b
        init_guess_params((n_types+4):(2*n_types+3)) = 0.5;
        lb((n_types+4):(2*n_types+3)) = -1000;
        ub((n_types+4):(2*n_types+3)) = 1000;

    % p
        init_guess_params((2*n_types+4):(3*n_types+3)) = 0.5;
        lb((2*n_types+4):(3*n_types+3)) = -1000;
        ub((2*n_types+4):(3*n_types+3)) = 1000;
        
% Apply fmincon
    objfun_handle = @(param) objfun_eNDM_general_costopts(param,...
                    seed_location,pathology,time_stamps,C,U,solvetype,costfun);
    options = optimoptions('fmincon','MaxFunctionEvaluations',10000);
    [param_num, fval_num] = fmincon(objfun_handle,init_guess_params,[],[],[],[],lb,ub,[],options);
 
% Solve eNDM with the optimal parameters
    x0_num = seed_location*param_num(1);
    alpha_num = param_num(2);
    beta_num = param_num(3);
    a_num = param_num(4:(n_types+3));
    b_num = param_num((n_types+4):(2*n_types+3));
    p_num = param_num((2*n_types+4):(3*n_types+3));
    ynum = eNDM_general(x0_num,time_stamps,C,U,alpha_num,beta_num,a_num,b_num,p_num,solvetype);
      
 %% Display Quantitative Results  
 % Save Rvalues in a matrix
    Rvalues = zeros(1,length(time_stamps));
    for jj = 1:length(time_stamps)
        if strcmp(corrtype,'R')
            Rvalues(jj) = corr(ynum(:,jj),pathology(:,jj), 'rows','complete');
        elseif strcmp(corrtype,'R_c')
            naninds = isnan(pathology(:,1));
            newxt = ynum; newxt(naninds,:) = [];
            newpath = pathology; newpath(naninds,:) = [];
            Rvalues(jj) = LinRcalc(newxt(:,jj),newpath(:,jj));
        end
    end

 % Display results
    disp('--------------------------------------------------')
    disp('General eNDM minimizing quadratic error at all time stamps with fmincon')   
    disp(' ')
    disp(['Optimal seed rescale value = ' num2str(param_num(1))])
    disp(['Optimal alpha = ' num2str(param_num(2))])
    disp(['Optimal beta = ' num2str(param_num(3))])
    disp(' ')
    for i = 1:length(type_inds)      
        disp(['Cell Type - ' classkey{type_inds(i)}])
        disp(['Optimal a = ' num2str(param_num(3 + i))])
        disp(['Optimal b = ' num2str(param_num(3 + n_types + i))])
        disp(['Optimal p = ' num2str(param_num(3 + 2*n_types + i))])
        disp(' ')
    end
    disp([corrtype ' values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    if strcmp(costfun,'LinR')
        disp(['Cost Function = ' num2str(length(time_stamps)) ' - sum(LinR)']) 
    else
        disp(costfun)
    end
    disp(fval_num)
    disp(' ')
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(ynum,pathology,time_stamps)
clearvars y 

