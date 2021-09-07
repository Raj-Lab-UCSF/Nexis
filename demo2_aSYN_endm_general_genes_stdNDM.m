%%  0.  Setup

% Initialization
clear; close all; clc;

% Add folder with raw data to path
% addpath('/Users/samuelteshome/eNDM/raw_data_mouse')
% % Add library with eNDM functions to path
% addpath('/Users/samuelteshome/eNDM/lib_eNDM_analytic')
% addpath('/Users/samuelteshome/eNDM/lib_eNDM_numeric')
% addpath('/Users/samuelteshome/eNDM/lib_eNDM_general')
% addpath('/Users/samuelteshome/Desktop/eNDM_aSynData')
% Load dataset of interest
%study = 'Gene';
load('mouse_aSynData_426.mat');
%load('eNDM_mousedata');
% load('CellNames.mat');
% load('CellTypeData.mat');
load('Regional_Gene_Data.mat');
% classkey = classkey_Tasic;
%load mouse_aSynData_426.mat
count = 0;
% Specify cost function, options: 'sse_sum', 'sse_end', 'rval_sum',
% 'rval_end', 'LinR', 'LinR_end
array = ["LinR", "LinR_end"]; 
filename = 'init parameter changes.xlsx'; %didn't use this

%T = table('Size',[10 8],'VariableTypes',"single");
newtable = array2table(zeros(150,9), 'VariableNames',{'costfunctionname', 'datasetname', 'initsr', 'initA', 'initB', 'SR', 'alpha', 'beta', 'rsq'});
%instantiated an array with defined rows and columns 
%had a lot of trouble doing this without making it an array first
for kj = 1 %nested for loop for costfunctions
disp(convertStringsToChars(array(kj)));
costfun = convertStringsToChars(array(kj));

% Defining LinR and type of correlation to display, R_c orR
LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);
corrtype = 'R_c';

% Select connectome C (426x426), use only ND (symmetric) for now
%C = Networks.nd;        % symmmetric
C = Networks.ret;      % non-symmetric
%C = Networks.ant;     % non-symmetric

% Normalize C
cmax = max(max(C));
cmin = min(min(C));
C = (C - cmin)./(cmax-cmin);

% Define number of regions of interest (nroi)
nroi = size(C,1);

% Define cell/gene type indices (from classkey)
type_inds = 1

% Define multiple cell/gene type indices (from classkey)
%     type_inds = [2344, 801, 410, 1162, 2305];
% Trem2 = 3578 (activated microglia (+ve regulator of Aβ phagocytosis))
% CD33 = 519 (activated microglia (-ve regulator of Aβ phagocytosis))
% P2ry12 = 2344 (microglial marker
% Cx3cr1 = 801 (microglial marker)
% C1qa = 410 (microglial marker)
% Fcrls = 1162 (microglial marker)
% Olfml3 = 2305 (microglial marker)
% Spp1 = 3260 (upregulated in microglia that surround plaques in AD models)
% Apoe = 213 (upregulated in microglia that surround plaques in AD models)
% Lpl = 1896 (upregulated in microglia that surround plaques in AD models)

% Define cell type matrix, U
U_gene = regvgene_mean(:,type_inds);
U = U_gene./nanmean(U_gene);
U = zeros(426,1);
%U = U - repmat(nanmean(U),size(U,1),1);

% mean of U
U_mean = mean(U,2);

% Z score
U_Z = (U - mean(U, 2))/std(U);

%
%PCA
%     [coeff, score, ~, ~, explained] = pca(U);
%     U = score(:,1);
%     %U = U - min(U);
%     %U = U / norm(U);
%     %U = U - U_mean;
%     type_inds = 1;

% Analytic vs. numeric ODE implementation
solvetype = 'analytic';

% Load time stamps, pathology measurements, and the seed_location

%Human or Mouse
array2 = ["mouse", "human"];

for kj2 = 1 %second nested for loop for dataset type
disp(convertStringsToChars(array2(kj2)));
datsetname = convertStringsToChars(array2(kj2));
%datsetname = 'mouse';
% Mouse pathology data inputs based on datsetname
time_stamps = tpts.(datsetname);
pathology = data426.(datsetname);
%      pathology = pathology./norm(pathology,2);

%pathology = pathology/sum(pathology(~isnan(pathology(:,1))));
%pathology = pathology+1;
%pathology = log(pathology);
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

%Set initial guesses and bounds for parameters
n_types = length(type_inds);
init_guess_params = zeros((3*n_types+3),1);
ub = init_guess_params; lb = ub;
seedrescalevalues = [(nansum(pathology(:,1)))/nnz(seed_location), .25, .5 , 2 ]; %(nansum(pathology(:,1)))/nnz(seed_location)
Alphavalues = [0.5, .1, .25 .75];
betavalues = [1, .25, .5, 2];
%cycled through every combination of 4 different init values for SR,A, & B
%I used 3 nested for loops to do this, might've been an easier way but
%wasn't sure how to do it another way. 5 nested for loops added overall

for srlength = 1
% seed rescale value
init_guess_params(1) = seedrescalevalues(srlength); 
lb(1) = 0;
ub(1) = inf;
sprintf('init seedrescale');
disp(seedrescalevalues(srlength));
sprintf("init seedrescale:")
initsr = seedrescalevalues(srlength); %value added to table
%nansum(pathology(:,1)); % orginally inf
%(.25, .5 , 2)


for alength = 1
% alpha
%init_guess_params(2) = log(nansum(pathology(:,end))/nansum(pathology(:,1)))/(tpts.(datsetname)(end)- tpts.(datsetname)(1));%oringally 0.5
init_guess_params(2) = Alphavalues(alength);
lb(2) = 0;
%ub(2) = 2 * log(nansum(pathology(:,end))/nansum(pathology(:,1)))/(tpts.(datsetname)(end)- tpts.(datsetname)(1));  % originally inf
ub(2) = inf;
sprintf("init alpha:")
disp(Alphavalues(alength));
initA = Alphavalues(alength); %value added to table
%(.1, .25 .75)


for blength = 1
% beta
%init_guess_params(3) = rand * randi(10);
init_guess_params(3) = betavalues(blength);
lb(3) = 0; % beta must be strictly non-negative
ub(3) = Inf;
%(.25, .5, 2)
sprintf("init beta:")
disp(betavalues(blength));
initB = betavalues(blength); %value added to table


% a
init_guess_params(4:(n_types+3)) = 0; % must be set to 0 for solvetype = 'analytic'
lb(4:(n_types+3)) = 0; % must be set to 0 for solvetype = 'analytic'
ub(4:(n_types+3)) = 0; % must be set to 0 for solvetype = 'analytic'

% b
init_guess_params((n_types+4):(2*n_types+3)) = 0;
lb((n_types+4):(2*n_types+3)) = 0;
ub((n_types+4):(2*n_types+3)) = 0;

% p
init_guess_params((2*n_types+4):(3*n_types+3)) = 0;
lb((2*n_types+4):(3*n_types+3)) = 0;
ub((2*n_types+4):(3*n_types+3)) = 0;

% Apply fmincon
objfun_handle = @(param) objfun_eNDM_general_costopts(param,...
    seed_location,pathology,time_stamps,C,U,solvetype,costfun);
%   options = optimoptions(@fmincon, 'Display', 'iter','Algorithm', 'sqp','MaxFunctionEvaluations',10000,'OptimalityTolerance',1e-10);
options = optimoptions(@fmincon,'Display', 'iter', 'MaxFunctionEvaluations',10000,'OptimalityTolerance',1e-8);
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

%convert the matrix of data and timepoints to a single column vector
P = reshape(pathology, [], 1);
Y = reshape(ynum, [], 1);
numObs1 = length(P(~isnan(P)))

% morder has to be 1 + the number of estimated parameters in the model
morder = 1 + 3;

%fit linear model and extract model criteria to calculate AIC and BIC
%lm_endm = fitlm(Y, P, 'Intercept', false);
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
output = cell(length(type_inds)+1, 18);
output(1,:) = {'Gene', 'Seed rescale', 'Alpha', 'Beta', 'a', 'b', 'p', 'Cost Func', 'R_time_pt_1', 'R_time_pt_2', 'R_time_pt_3',...
    'logL', 'endm_aic', 'endm_bic', 'endm_intercept', 'endm_p', 'endm_Rsqr_ord', 'endm_Rsqr_adj'};
disp('--------------------------------------------------')
disp('General eNDM minimizing quadratic error at all time stamps with fmincon')
disp(' ')
disp(['Optimal seed rescale value = ' num2str(param_num(1))])
disp(['Optimal alpha = ' num2str(param_num(2))])
disp(['Optimal beta = ' num2str(param_num(3))])
disp(' ')

for i = 1:length(type_inds)
%     disp(['Cell Type - ' classkey{type_inds(i)}])
    disp(['Optimal a = ' num2str(param_num(3 + i))])
    disp(['Optimal b = ' num2str(param_num(3 + n_types + i))])
    disp(['Optimal p = ' num2str(param_num(3 + 2*n_types + i))])
    disp(' ')
    
    output(i+1,:)={'test', param_num(1), param_num(2), param_num(3), param_num(3 + i),...
        param_num(3 + n_types + i),param_num(3 + 2*n_types + i), fval_num, Rvalues(1), Rvalues(2), Rvalues(3), ...
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
count = count + 1;
SR = param_num(1);
alpha = param_num(2);
beta = param_num(3);
rsq = endm_Rsqr_ord;
%organized variables that would go into table above

if array(kj) == "LinR"
costfunctionname = 1000;
else 
    costfunctionname = 1001;
end

if datsetname == "mouse"
    datasetname = 1002;
    
else
    datasetname = 1003;
end

%T = table(costfunctionname, initsr, initA, initB, SR, alpha, beta, rsq);
newtable(count,:) = table(costfunctionname, datasetname, initsr, initA, initB, SR, alpha, beta, rsq);
%this is where the table is updated
% Plot prediction vs data using optimal parameters
initSeed = num2str(initsr);
initAlpha = num2str(initA);
initBeta = num2str(initB);
%converted to string so it would display on graph ^

% plot_pred_vs_data1(ynum, pathology, time_stamps, costfun, datsetname, initSeed, initAlpha, initBeta);
% plot_pred_vs_data(predicted,data,time_stamps,costfun, initsr, initA, initB)
clearvars y
numObs1 = length(P(~isnan(P)))
Correlation = corrcoef(Y, P,'Rows','complete')
logL;
lm_endm
%xlswrite('output.xlsx', output);
% [~,indx]=ismember('Trem2',gene_names_trans)

flow = FlowCalculator(ynum,C,beta_num,1);
%flow(flow < prctile(nonzeros(flow),90)) = 0;


for i = 1:size(flow,3)
flow_ = flow(:,:,i);
flow_(flow_ < prctile(nonzeros(flow),99.9)) = 0;
flow(:,:,i) = flow_;
end
end
end
end
end
end