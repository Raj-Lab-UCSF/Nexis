%% 1. eNDM_mousedata.mat datasets

% Add library with eNDM functions to path
addpath('/Users/justintorok/Documents/MATLAB/Nexis/MATLAB/lib_eNDM_general')
datapath = '/Users/justintorok/Documents/MATLAB/Nexis/raw_data_mouse';

% Load dataset of interest
% load([datapath filesep 'Tasic_606.mat'],'outstruct_606','classkey');
load([datapath filesep 'eNDM_mousedata.mat']);

% Specify cost function and others, options: 'sse_sum', 'sse_end', 'rval_sum',
% 'rval_end', 'LinR'
costfun = 'LinR';
volcorrect = 1;
exclseed_costfun = 0;
excltpts_costfun = [];

% Defining LinR and type of correlation to display, R_c or R
LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);
corrtype = 'R_c';

% Use standard directional connectome but will symmetrize
C = Networks.ret;      % non-symmetric

% Normalize C
cmax = max(max(C));
cmin = min(min(C));
C = (C - cmin)./(cmax-cmin);

% Define number of regions of interest (nroi)
nroi = size(C,1);
U = zeros(nroi,1); % no cell type information

% Analytic vs. numeric ODE implementation
solvetype = 'analytic';
nexis_global_tau_struct = struct;
datsetnames = {'BolundaDSAD','BolundaCBD','IbaHippInj','IbaStrInj','Clavaguera','Hurtado'};
for i = 1:length(datsetnames)
    datsetname = datsetnames{i};
    fprintf('eNDM Mouse Dataset %d of %d\n',i,length(datsetnames))
    % Mouse pathology data inputs based on datsetname
    time_stamps = tpts.(datsetname);
    pathology = data426.(datsetname);

    pathology = pathology/sum(pathology(~isnan(pathology(:,1))));
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
    n_types = size(U,2);
    init_guess_params = zeros((3*n_types+3),1);
    ub = init_guess_params; lb = ub;

    % seed rescale value
    if strcmp(datsetname,'Hurtado')
        init_guess_params(1) = (nansum(pathology(:,1)))/nnz(seed_location);
        lb(1) = 0;
        ub(1) = Inf;    
    else
        init_guess_params(1) = (nansum(pathology(:,1)))/nnz(seed_location);
        lb(1) = 0;
        ub(1) = Inf;    
    end


    % alpha
    alpha_0 = log(nansum(pathology(:,end))/nansum(pathology(:,1)))/(time_stamps(end) - time_stamps(1));
    init_guess_params(2) = alpha_0;
    lb(2) = 0;
    ub(2) = 2*alpha_0;

    % beta
    init_guess_params(3) = 1;
    lb(3) = 0; % beta must be strictly non-negative
    ub(3) = Inf;
    
    % s - fixed at 0.5
    init_guess_params(4) = 0.5;
    lb(4) = 0.5;
    ub(4) = 0.5;

    % a
    init_guess_params(5:(n_types+4)) = 0; % must be set to 0 for solvetype = 'analytic'
    lb(5:(n_types+4)) = 0; % must be set to 0 for solvetype = 'analytic'
    ub(5:(n_types+4)) = 0; % must be set to 0 for solvetype = 'analytic'

    % b
    init_guess_params((n_types+5):(2*n_types+4)) = 0;
    lb((n_types+5):(2*n_types+4)) = 0;
    ub((n_types+5):(2*n_types+4)) = 0;

    % p
    init_guess_params((2*n_types+5):(3*n_types+4)) = 0;
    lb((2*n_types+5):(3*n_types+4)) = 0;
    ub((2*n_types+5):(3*n_types+4)) = 0;

    % Apply fmincon
    objfun_handle = @(param) objfun_eNDM_general_dir_costopts(param,...
        seed_location,pathology,time_stamps,C,U,solvetype,volcorrect,costfun,excltpts_costfun,exclseed_costfun);
    %   options = optimoptions(@fmincon, 'Display', 'iter','Algorithm', 'sqp','MaxFunctionEvaluations',10000,'OptimalityTolerance',1e-10);
    options = optimoptions(@fmincon,'Display', 'iter', 'MaxFunctionEvaluations',10000,'OptimalityTolerance',1e-8);
    [param_num, fval_num] = fmincon(objfun_handle,init_guess_params,[],[],[],[],lb,ub,[],options);

    % Solve eNDM with the optimal parameters
    x0_num = seed_location*param_num(1);
    alpha_num = param_num(2);
    beta_num = param_num(3);
    s_num = param_num(4);
    a_num = param_num(5:(n_types+4));
    b_num = param_num((n_types+5):(2*n_types+4));
    p_num = param_num((2*n_types+5):(3*n_types+4));
    ynum = eNDM_general_dir(x0_num,time_stamps,C,U,alpha_num,beta_num,s_num,a_num,b_num,p_num,solvetype);

    %convert the matrix of data and timepoints to a single column vector
    P = reshape(pathology, [], 1);
    Y = reshape(ynum, [], 1);

    lm_endm = fitlm(Y, P);
    endm_p = lm_endm.Coefficients.pValue(1);
    endm_Rsqr_adj = lm_endm.Rsquared.Adjusted; 
    
    % Store relevant output
    nexis_global_tau_struct.(datsetname).t = time_stamps;
    nexis_global_tau_struct.(datsetname).seed = seed_location;
    nexis_global_tau_struct.(datsetname).data = pathology;
    nexis_global_tau_struct.(datsetname).predicted = ynum;
    nexis_global_tau_struct.(datsetname).seed_rescale = param_num(1);
    nexis_global_tau_struct.(datsetname).alpha = param_num(2);
    nexis_global_tau_struct.(datsetname).beta = param_num(3);
    nexis_global_tau_struct.(datsetname).s = param_num(4);
    nexis_global_tau_struct.(datsetname).seed_rescale_0 = init_guess_params(1);
    nexis_global_tau_struct.(datsetname).alpha_0 = init_guess_params(2);
    nexis_global_tau_struct.(datsetname).beta_0 = init_guess_params(3);
    nexis_global_tau_struct.(datsetname).s_0 = param_num(4);
    nexis_global_tau_struct.(datsetname).lm_R2 = endm_Rsqr_adj;
    nexis_global_tau_struct.(datsetname).lm_pval = endm_p;
end
clearvars -except nexis_global_tau_struct

%% 2. Diamond DS datasets

% Add library with eNDM functions to path
addpath('/Users/justintorok/Documents/MATLAB/Nexis/MATLAB/lib_eNDM_general')
datapath = '/Users/justintorok/Documents/MATLAB/Nexis/raw_data_mouse';

% Load dataset of interest
% load([datapath filesep 'Tasic_606.mat'],'outstruct_606','classkey');
load([datapath filesep 'eNDM_mousedata.mat'],'Networks');
load([datapath filesep 'KaufmanDiamond_datasets_dat&seed.mat']);

% Specify cost function and others, options: 'sse_sum', 'sse_end', 'rval_sum',
% 'rval_end', 'LinR'
costfun = 'LinR';
volcorrect = 1;
exclseed_costfun = 0;
excltpts_costfun = [];

% Defining LinR and type of correlation to display, R_c or R
LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);
corrtype = 'R_c';

% Use standard directional connectome but will symmetrize
C = Networks.ret;      % non-symmetric

% Normalize C
cmax = max(max(C));
cmin = min(min(C));
C = (C - cmin)./(cmax-cmin);

% Define number of regions of interest (nroi)
nroi = size(C,1);
U = zeros(nroi,1); % no cell type information

% Analytic vs. numeric ODE implementation
solvetype = 'analytic';
datsetnames = fieldnames(data426);
for i = 1:length(datsetnames)
    datsetname = datsetnames{i};
    fprintf('DS Dataset %d of %d\n',i,length(datsetnames))
    % Mouse pathology data inputs based on datsetname
    time_stamps = tpts.(datsetname);
    pathology = data426.(datsetname);

    pathology = pathology/sum(pathology(~isnan(pathology(:,1))));
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
    n_types = size(U,2);
    init_guess_params = zeros((3*n_types+3),1);
    ub = init_guess_params; lb = ub;

    % seed rescale value
    if strcmp(datsetname,'Hurtado')
        init_guess_params(1) = (nansum(pathology(:,1)))/nnz(seed_location);
        lb(1) = 0;
        ub(1) = Inf;    
    else
        init_guess_params(1) = (nansum(pathology(:,1)))/nnz(seed_location);
        lb(1) = 0;
        ub(1) = Inf;    
    end


    % alpha
    alpha_0 = log(nansum(pathology(:,end))/nansum(pathology(:,1)))/(time_stamps(end) - time_stamps(1));
    init_guess_params(2) = alpha_0;
    lb(2) = 0;
    ub(2) = 2*alpha_0;

    % beta
    init_guess_params(3) = 1;
    lb(3) = 0; % beta must be strictly non-negative
    ub(3) = Inf;
    
    % s - fixed at 0.5
    init_guess_params(4) = 0.5;
    lb(4) = 0.5;
    ub(4) = 0.5;

    % a
    init_guess_params(5:(n_types+4)) = 0; % must be set to 0 for solvetype = 'analytic'
    lb(5:(n_types+4)) = 0; % must be set to 0 for solvetype = 'analytic'
    ub(5:(n_types+4)) = 0; % must be set to 0 for solvetype = 'analytic'

    % b
    init_guess_params((n_types+5):(2*n_types+4)) = 0;
    lb((n_types+5):(2*n_types+4)) = 0;
    ub((n_types+5):(2*n_types+4)) = 0;

    % p
    init_guess_params((2*n_types+5):(3*n_types+4)) = 0;
    lb((2*n_types+5):(3*n_types+4)) = 0;
    ub((2*n_types+5):(3*n_types+4)) = 0;

    % Apply fmincon
    objfun_handle = @(param) objfun_eNDM_general_dir_costopts(param,...
        seed_location,pathology,time_stamps,C,U,solvetype,volcorrect,costfun,excltpts_costfun,exclseed_costfun);
    %   options = optimoptions(@fmincon, 'Display', 'iter','Algorithm', 'sqp','MaxFunctionEvaluations',10000,'OptimalityTolerance',1e-10);
    options = optimoptions(@fmincon,'Display', 'iter', 'MaxFunctionEvaluations',10000,'OptimalityTolerance',1e-8);
    [param_num, fval_num] = fmincon(objfun_handle,init_guess_params,[],[],[],[],lb,ub,[],options);

    % Solve eNDM with the optimal parameters
    x0_num = seed_location*param_num(1);
    alpha_num = param_num(2);
    beta_num = param_num(3);
    s_num = param_num(4);
    a_num = param_num(5:(n_types+4));
    b_num = param_num((n_types+5):(2*n_types+4));
    p_num = param_num((2*n_types+5):(3*n_types+4));
    ynum = eNDM_general_dir(x0_num,time_stamps,C,U,alpha_num,beta_num,s_num,a_num,b_num,p_num,solvetype);

    %convert the matrix of data and timepoints to a single column vector
    P = reshape(pathology, [], 1);
    Y = reshape(ynum, [], 1);

    lm_endm = fitlm(Y, P);
    endm_p = lm_endm.Coefficients.pValue(1);
    endm_Rsqr_adj = lm_endm.Rsquared.Adjusted; 
    
    % Store relevant output
    nexis_global_tau_struct.(datsetname).t = time_stamps;
    nexis_global_tau_struct.(datsetname).seed = seed_location;
    nexis_global_tau_struct.(datsetname).data = pathology;
    nexis_global_tau_struct.(datsetname).predicted = ynum;
    nexis_global_tau_struct.(datsetname).seed_rescale = param_num(1);
    nexis_global_tau_struct.(datsetname).alpha = param_num(2);
    nexis_global_tau_struct.(datsetname).beta = param_num(3);
    nexis_global_tau_struct.(datsetname).s = param_num(4);
    nexis_global_tau_struct.(datsetname).seed_rescale_0 = init_guess_params(1);
    nexis_global_tau_struct.(datsetname).alpha_0 = init_guess_params(2);
    nexis_global_tau_struct.(datsetname).beta_0 = init_guess_params(3);
    nexis_global_tau_struct.(datsetname).s_0 = param_num(4);
    nexis_global_tau_struct.(datsetname).lm_R2 = endm_Rsqr_adj;
    nexis_global_tau_struct.(datsetname).lm_pval = endm_p;
end
save('nexis_global_tau_struct.mat','nexis_global_tau_struct');