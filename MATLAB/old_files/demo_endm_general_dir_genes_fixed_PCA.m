%%  0.  Setup

% Initialization
clear; close all; clc;

% Add folder with raw data to path
addpath('C:/Users/chait/Documents/NOTES/CHAI/UCSF/Raj_Lab/NDM/eNDM-master_genes_new/raw_data_mouse')

% Add library with eNDM functions to path
addpath('C:/Users/chait/Documents/NOTES/CHAI/UCSF/Raj_Lab/NDM/eNDM-master_genes_new/lib_eNDM_analytic')
addpath('C:/Users/chait/Documents/NOTES/CHAI/UCSF/Raj_Lab/NDM/eNDM-master_genes_new/lib_eNDM_numeric')

% Load dataset of interest
study = 'Gene';
load('Regional_Gene_Data.mat');
load('gene_names_trans.mat','gene_names_trans');
classkey = gene_names_trans;
load eNDM_mousedata.mat

% Specify cost function, options: 'sse_sum', 'sse_end', 'rval_sum',
% 'rval_end', 'LinR'
costfun = 'LinR';

% Defining LinR and type of correlation to display, R_c or R
LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);
corrtype = 'R_c';

% Select connectome C (426x426)
% C = Networks.nd;        % symmmetric
C = Networks.ret;      % non-symmetric, for use with directionality
%C = Networks.ant;     % non-symmetric

% Normalize C
cmax = max(max(C));
cmin = min(min(C));
C = (C - cmin)./(cmax-cmin);

% Define number of regions of interest (nroi)
nroi = size(C,1);

% Define cell/gene type indices (from classkey)
%  type_inds = [1463];

% Define multiple cell/gene type indices (from classkey)
type_inds = [2344, 801, 410, 1162, 2305];
% Trem2 = 3578 (activated microglia (+ve regulator of Aβ phagocytosis))
% CD33 = 519 (activated microglia (-ve regulator of Aβ phagocytosis))
% P2ry12 = 2344 (microglial marker
% Cx3cr1 = 801 (microglial marker)
% C1qa = 410 (microglial marker)
% Fcrl5 = 1162 (microglial marker)
% Olfml3 = 2305 (microglial marker)
% Tspo = 3629 (myeloid cell marker, PET ligand)
% Spp1 = 3260 (upregulated in microglia that surround plaques in AD models)
% Apoe = 213 (upregulated in microglia that surround plaques in AD models)
% Lpl = 1896 (upregulated in microglia that surround plaques in AD models)

% Define cell type matrix, U_gene, and normalize it
U_gene = regvgene_mean(:,type_inds);
U = U_gene./nanmean(U_gene);

%U = U - repmat(nanmean(U),size(U,1),1);

% mean of each row of U matrix thus giving one vector
U_mean = mean(U,2);

% Z score (of what?)
% U_Z = (U - U_mean)/std(U);

%%PCA%%%
[coeff, score, ~, ~, explained] = pca(U);
U = score(:,1);
if corr(score(:,1),U_mean) < 0
    U = -U;
end
U = (U - min(U))/(max(U) - min(U));
%U = U - min(U);
%U = U / norm(U);
%U = U - U_mean;
type_inds = 1;

% Analytic vs. numeric ODE implementation
solvetype = 'analytic';

%% Load datasets, time stamps, pathology measurements, and the seed_location

% the following arrays are initial guesses per global parameter for each
% dataset as achieved by running the stdNDM

seed_rescale_array = [0.6269, 0.0572, 0.0535, 0.0405];
alpha_array = [0.4238, 0.3956, 0.1487, 0.5254];
beta_array = [2.7803, 2.6757, 2.9797, 0.6614];
s_array = [0.5, 0.5, 0.5, 0.5]; %CHANGE THESE

dataset = {'IbaHippInj','IbaStrInj','Clavaguera','Hurtado'};
out_final = cell(length(dataset)+1, 19);

for dataset_itr = 1:length(dataset)
    datsetname = dataset{dataset_itr};
    % Mouse pathology data inputs based on datsetname
    time_stamps = tpts.(datsetname);
    pathology = data426.(datsetname);
    %      pathology = pathology./norm(pathology,2);
    
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
    n_types = length(type_inds);
    init_guess_params = zeros((3*n_types+4),1);
    ub = init_guess_params; lb = ub;
    
    % seed rescale value
    init_guess_params(1) = seed_rescale_array(dataset_itr); %(nansum(pathology(:,1)))/nnz(seed_location);
    lb(1) = seed_rescale_array(dataset_itr);
    ub(1) = seed_rescale_array(dataset_itr);
    
    % alpha (ub and lb 10% of initial guess)
    init_guess_params(2) = alpha_array(dataset_itr); %0.5
    lb(2) = 0.9*alpha_array(dataset_itr); %0
    ub(2) = 1.1*alpha_array(dataset_itr); %Inf
    
    % beta (beta must be strictly non-negative)(ub and lb 10% of initial guess)
    init_guess_params(3) = beta_array(dataset_itr); %1
    lb(3) = 0.9*beta_array(dataset_itr);  %0
    ub(3) = 1.1*beta_array(dataset_itr); %Inf
    
    % s (s must be bounded between 0 and 1)
    init_guess_params(4) = s_array(dataset_itr); %1
    lb(4) = s_array(dataset_itr);  %0
    ub(4) = s_array(dataset_itr); %Inf
    
    % a
    init_guess_params(5:(n_types+4)) = 0; % must be set to 0 for solvetype = 'analytic'
    lb(5:(n_types+4)) = 0; % must be set to 0 for solvetype = 'analytic'
    ub(5:(n_types+4)) = 0; % must be set to 0 for solvetype = 'analytic'
    
    % b
    init_guess_params((n_types+5):(2*n_types+4)) = 0;
    lb((n_types+5):(2*n_types+4)) = -Inf;
    ub((n_types+5):(2*n_types+4)) = Inf;
    
    % p
    init_guess_params((2*n_types+5):(3*n_types+4)) = 0;
    lb((2*n_types+5):(3*n_types+4)) = -Inf;
    ub((2*n_types+5):(3*n_types+4)) = Inf;
    
    
    % Apply fmincon
    objfun_handle = @(param) objfun_eNDM_general_dir_costopts(param,...
        seed_location,pathology,time_stamps,C,U,solvetype,costfun);
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
    
    %convert the matrix of data and timepoints to a single column vector
    P = reshape(pathology, [], 1);
    Y = reshape(ynum, [], 1);
    numObs1 = length(P(~isnan(P)))
    
    % morder has to be 1 + the number of estimated parameters in the model
    morder = 1 + 4 + 2*length(type_inds);
    
    %fit linear model and extract model criteria to calculate AIC and BIC
    lm_endm = fitlm(Y, P);
    lm_endm.ModelCriterion;
    logL = lm_endm.LogLikelihood;
    endm_aic = -2*logL + 2*morder;
    endm_bic = -2*logL + log(numObs1)*morder;
    endm_intercept = lm_endm.Coefficients.Estimate(1);
    endm_p = lm_endm.Coefficients.pValue(1);
    endm_Rsqr_ord = lm_endm.Rsquared.Ordinary;
    endm_Rsqr_adj = lm_endm.Rsquared.Adjusted;
    
    % Display results
    output = cell(length(type_inds)+1, 19);
    output(1,:) = {'Gene', 'Seed rescale', 'Alpha', 'Beta', 's', 'a', 'b', 'p', 'Cost Func', 'R_time_pt_1', 'R_time_pt_2', 'R_time_pt_3',...
        'logL', 'endm_aic', 'endm_bic', 'endm_intercept', 'endm_p', 'endm_Rsqr_ord', 'endm_Rsqr_adj'};
    disp('--------------------------------------------------')
    disp('General eNDM minimizing quadratic error at all time stamps with fmincon')
    disp(' ')
    disp(['Optimal seed rescale value = ' num2str(param_num(1))])
    disp(['Optimal alpha = ' num2str(param_num(2))])
    disp(['Optimal beta = ' num2str(param_num(3))])
    disp(['Optimal s = ' num2str(param_num(4))])
    disp(' ')

    for i = 1:length(type_inds)
        disp(['Cell Type - ' classkey{type_inds(i)}])
        disp(['Optimal a = ' num2str(param_num(4 + i))])
        disp(['Optimal b = ' num2str(param_num(4 + n_types + i))])
        disp(['Optimal p = ' num2str(param_num(4 + 2*n_types + i))])
        disp(' ')

        output(i+1,:)={classkey{type_inds(i)}, param_num(1), param_num(2), param_num(3), param_num(4), param_num(4 + i),...
            param_num(4 + n_types + i),param_num(4 + 2*n_types + i), fval_num, Rvalues(1), Rvalues(2), Rvalues(3), ...
            logL, endm_aic, endm_bic, endm_intercept, endm_p, endm_Rsqr_ord, endm_Rsqr_adj};
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
    disp(['AIC = ' num2str(endm_aic)])
    disp(['BIC = ' num2str(endm_bic)])
    disp(['Intercept = ' num2str(endm_intercept)])
    disp(['pValue = ' num2str(endm_p)])
    disp(['Rsqr_ord = ' num2str(endm_Rsqr_ord)])
    disp(['Rsqr_adj = ' num2str(endm_Rsqr_adj)])
    disp(' ')
    
    % Plot prediction vs data using optimal parameters
    plot_pred_vs_data_corr(ynum,pathology,time_stamps, datsetname)
    clearvars y
    lm_endm;
    numObs1 = length(P(~isnan(P)))
    Correlation = corrcoef(Y, P,'Rows','complete')
    logL;
    
    lm_endm
    out_final{dataset_itr+1,1} = datsetname;
    out_final(dataset_itr+1,2:end) = output(2,:);
    % xlswrite('output.xlsx', output);
%     [~,indx]=ismember('Tspo',gene_names_trans)
end
out_final(1,2:end) = output(1,:);
xlswrite('output.xlsx', out_final);