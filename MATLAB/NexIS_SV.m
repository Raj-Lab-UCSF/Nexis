function outputs = NexIS_SV(varargin)

% Parameters to fit
% param(1) = seed rescale factor
% param(2) = alpha
% param(3) = beta
% param(4) = s
% param(5) = b 
% param(6) = p

% Define defaults and set inputs
study_ = 'User-specified'; % double-check this
C_ = [];
data_ = [];
tpts_ = [];
seed_ = [];
use_dataspace_ = 0;
matdir_ = [cd filesep 'raw_data_mouse'];

% Define defaults and set inputs
costfun_ = 'linr'; % see /lib_NexIS/CostFunction_NexIS.m for options
lambda_ = 0; % see /lib_NexIS/CostFunction_NexIS.m for options
solvetype_ = 'analytic';
volcorrect_ = 1;
exclseed_costfun_ = 0;
excltpts_costfun_ = [];
logtrans_ = 'none';
Cnormtype_ = 'minmax';
normtype_ = 'sum';
w_dir_ = 0; 
param_init_ = [NaN,0.5,1,0.5];
ub_ = [Inf,Inf,Inf,1];
lb_ = zeros(1,4);
algo_ = 'sqp';
opttol_ = 1e-8;
fxntol_ = 1e-8;
steptol_ = 1e-12;
maxeval_ = 10000;
bootstrapping_ = 0;
resample_rate_ = 0.8;
niters_ = 100;
verbose_ = 0;
fmindisplay_ = 0;
outputs_nexisglobal_ = [];
bounds_type_nexis_sv_ = 'old'; % 'old', 'none', 'CI_X' where X is the percent
bootstrapping_nexis_sv_ = 0;
resample_rate_nexis_sv_ = 0.8;
niters_nexis_sv_ = 100;
verbose_nexis_sv_ = 0;
fmindisplay_nexis_sv_ = 0;
datatype_nexis_sv_ = 'gene'; % 'gene', 'Yao', 'Tasic', 'Zeisel', 
% 'Zhuang_Class', 'Zhuang_Subclass'
datalist_nexis_sv_ = {'Trem2'}; % if names are used, requires a cell array
% even for one element; can also supply numeric indices
datapca_nexis_sv_ = 0;
flowthresh_ = 99.93;

ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
validBoolean = @(x) isscalar(x) && (x==0 || x==1);
validChar = @(x) ischar(x);
validST = @(x) ismember(x,{'analytic','numeric'});
validBoundsType = @(x) strcmp(x,'old') || strcmp(x(1:2),'CI');
validParam = @(x) (length(x) == 4);
validDataTypenexis_sv = @(x) ismember(x,{'gene','Yao','Tasic','Zeisel',...
    'Zhuang_Class','Zhuang_Subclass'});

addParameter(ip, 'study', study_, validChar);
addParameter(ip, 'C', C_);
addParameter(ip, 'data',data_);
addParameter(ip, 'tpts',tpts_);
addParameter(ip, 'seed', seed_);
addParameter(ip, 'matdir', matdir_);
addParameter(ip, 'use_dataspace', use_dataspace_, validBoolean);
addParameter(ip, 'costfun', costfun_, validChar);
addParameter(ip, 'lambda', lambda_, validScalar);
addParameter(ip, 'solvetype', solvetype_, validST);
addParameter(ip, 'volcorrect', volcorrect_, validBoolean);
addParameter(ip, 'exclseed_costfun', exclseed_costfun_, validBoolean);
addParameter(ip, 'excltpts_costfun', excltpts_costfun_);
addParameter(ip, 'logtrans', logtrans_);
addParameter(ip, 'normtype', normtype_, validChar);
addParameter(ip, 'Cnormtype', Cnormtype_, validChar);
addParameter(ip, 'w_dir', w_dir_, validBoolean);
addParameter(ip, 'param_init', param_init_, validParam);
addParameter(ip, 'ub', ub_, validParam);
addParameter(ip, 'lb', lb_, validParam);
addParameter(ip, 'opttol', opttol_, validScalar);
addParameter(ip, 'steptol', steptol_, validScalar);
addParameter(ip, 'fxntol', fxntol_, validScalar);
addParameter(ip, 'algo', algo_, validChar);
addParameter(ip, 'maxeval', maxeval_, validScalar);
addParameter(ip, 'bootstrapping', bootstrapping_, validBoolean);
addParameter(ip, 'resample_rate', resample_rate_, validScalar);
addParameter(ip, 'niters', niters_, validScalar);
addParameter(ip, 'verbose', verbose_, validBoolean);
addParameter(ip, 'fmindisplay', fmindisplay_, validBoolean);
addParameter(ip, 'outputs_nexisglobal', outputs_nexisglobal_);
addParameter(ip, 'bounds_type_nexis_sv', bounds_type_nexis_sv_, validBoundsType);
addParameter(ip, 'bootstrapping_nexis_sv', bootstrapping_nexis_sv_, validBoolean);
addParameter(ip, 'resample_rate_nexis_sv', resample_rate_nexis_sv_, validScalar);
addParameter(ip, 'niters_nexis_sv', niters_nexis_sv_, validScalar);
addParameter(ip, 'verbose_nexis_sv', verbose_nexis_sv_, validBoolean);
addParameter(ip, 'fmindisplay_nexis_sv', fmindisplay_nexis_sv_, validBoolean);
addParameter(ip, 'datatype_nexis_sv', datatype_nexis_sv_, validDataTypenexis_sv);
addParameter(ip, 'datalist_nexis_sv', datalist_nexis_sv_);
addParameter(ip, 'datapca_nexis_sv', datapca_nexis_sv_, validBoolean);
addParameter(ip, 'flowthresh', flowthresh_, validScalar);

parse(ip, varargin{:});
ipR = ip.Results;

outputs = ipR.outputs_nexisglobal;
if isempty(outputs)
   outputs = NexIS_global('study', ipR.study,...
                          'C', ipR.C,...
                          'data', ipR.data,...
                          'tpts', ipR.tpts,...
                          'seed', ipR.seed,...
                          'use_dataspace', ipR.use_dataspace,...
                          'matdir', ipR.matdir,...
                          'lambda', ipR.lambda,...
                          'logtrans', ipR.logtrans,...
                          'costfun', ipR.costfun,...
                          'solvetype', ipR.solvetype,...
                          'volcorrect', ipR.volcorrect,...
                          'exclseed_costfun',ipR.exclseed_costfun,...
                          'excltpts_costfun',ipR.excltpts_costfun,...
                          'normtype', ipR.normtype,...
                          'w_dir', ipR.w_dir,...
                          'param_init', ipR.param_init,...
                          'ub', ipR.ub,...
                          'lb', ipR.lb,...
                          'bootstrapping', ipR.bootstrapping,...
                          'resample_rate', ipR.resample_rate,...
                          'niters', ipR.niters,...
                          'verbose', ipR.verbose,...
                          'fmindisplay', ipR.fmindisplay);
end

% Load in data from NexIS/raw_data_mouse directory if needed. Copied from
% NexIS_global.m so should be self-consistent, though could load these from
% 'outputs' above. Would require a bit of careful coding because 'data' in
% 'outputs' are normalized, though.
if isempty(ipR.C) && ~strcmp(ipR.study,'User_specified')
    if (length(ipR.study) > 3) && strcmp(ipR.study(1:4),'asyn')
        load([ipR.matdir filesep 'Mouse_Synuclein_Data.mat'],...
            'mousedata_struct');
        load([ipR.matdir filesep 'Connectomes.mat'],'Connectomes');
        C = Connectomes.default;
        data426 = struct; data426.(ipR.study) = mousedata_struct.(ipR.study).data;
        seed426 = struct; seed426.(ipR.study) = mousedata_struct.(ipR.study).seed;
        tpts = struct; tpts.(ipR.study) = mousedata_struct.(ipR.study).time_stamps;
    elseif strcmp(ipR.study,'Henderson')
        load([ipR.matdir filesep 'Mouse_Synuclein_Data.mat'],...
            'mousedata_struct');
        load([ipR.matdir filesep 'Connectomes.mat'],'Connectomes');
        C = Connectomes.Henderson;
        data426 = struct; data426.(ipR.study) = mousedata_struct.(ipR.study).data;
        seed426 = struct; seed426.(ipR.study) = mousedata_struct.(ipR.study).seed;
        tpts = struct; tpts.(ipR.study) = mousedata_struct.(ipR.study).time_stamps;
    elseif strcmp(ipR.study(1:3),'GCI') || strcmp(ipR.study(1:3),'PFF')
        load([ipR.matdir filesep 'Mouse_Synuclein_Data.mat'],...
            'mousedata_struct');
        load([ipR.matdir filesep 'Connectomes.mat'],'Connectomes');
        C = Connectomes.Peng;
        data426 = struct; data426.(ipR.study) = mousedata_struct.(ipR.study).data;
        seed426 = struct; seed426.(ipR.study) = mousedata_struct.(ipR.study).seed;
        tpts = struct; tpts.(ipR.study) = mousedata_struct.(ipR.study).time_stamps;
    else
        load([ipR.matdir filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],...
            'mousedata_struct');
        load([ipR.matdir filesep 'Connectomes.mat'],'Connectomes');
        data426 = struct; data426.(ipR.study) = mousedata_struct.(ipR.study).data;
        seed426 = struct; seed426.(ipR.study) = mousedata_struct.(ipR.study).seed;
        tpts = struct; tpts.(ipR.study) = mousedata_struct.(ipR.study).time_stamps;
        C = Connectomes.default;
        if strcmp(ipR.study(1:2),'DS')
            tpts.(ipR.study) = tpts.(ipR.study)/4; % weeks to months for Kaufman studies
        end
    end
else
    data426 = struct; tpts = struct; seed426 = struct;
    data426.(ipR.study) = ipR.data; tpts.(ipR.study) = ipR.tpts; seed426.(ipR.study) = ipR.seed;
    C = ipR.C;
end

% Normalize C
if strcmp(ipR.Cnormtype,'eig')
    C = C/max(eig(C));
else
    cmax = max(max(C));
    cmin = min(min(C));
    C = (C - cmin)./(cmax-cmin);
end

% Define cell type matrix, U
if ~isequal(ipR.datalist_nexis_sv,{'random'})
    if ~isnumeric(ipR.datalist_nexis_sv)
        ind_nexis_sv = NameIndex(ipR.datalist_nexis_sv,ipR.datatype_nexis_sv);
    else
        ind_nexis_sv = ipR.datalist_nexis_sv;
    end

    if strcmp(ipR.datatype_nexis_sv,'gene')
        datstruct_gene = load([cd filesep 'raw_data_mouse' filesep ...
            'GeneExpressionMaps.mat'],'GeneExpressionMaps');
        U = datstruct_gene.GeneExpressionMaps.All.expression_426(:,ind_nexis_sv);        
        % load([cd filesep 'raw_data_mouse' filesep 'Regional_Gene_Data.mat'],'regvgene_mean');
        % U = regvgene_mean(:,ind_nexis_sv);
    else
        datstruct_ct = load([cd filesep 'raw_data_mouse' filesep 'CellTypeMaps.mat'],'CellTypeMaps');
        U = datstruct_ct.GeneExpressionMaps.(ipR.datatype_nexis_sv).maps(:,ind_nexis_sv);  
    end
else
    U = rand(size(C,1),1);
end

% Reorder U data if needed (IGNORING FOR NOW; should not use this code for
% Yuanxi's data at this point - 12/30/24)
if strcmp(ipR.study,'Henderson')
    error('No gene data on Henderson regional atlas at this time!')
elseif strcmp(ipR.study(1:3),'GCI') || strcmp(ipR.study(1:3),'PFF')
    load([cd filesep 'raw_data_mouse' filesep 'GCI_PFF_Data.mat'],...
        'BrainRegionReorderMat');
    load([cd filesep 'raw_data_mouse' filesep 'regionvoxels.mat'],...
        'voxels');
    BrainRegionReorderMat_2h = [BrainRegionReorderMat; ...
        (BrainRegionReorderMat+213)]; % add second hemisphere
    voxels_2h = [voxels; voxels];
    U_ = nan(size(BrainRegionReorderMat_2h,1),1);
    voxels_ = U_;
    for i = 1:length(U_)
        CCF_inds = BrainRegionReorderMat_2h(i,:);
        CCF_inds(isnan(CCF_inds)) = [];
        if ~isempty(CCF_inds)
            voxels_i = voxels_2h(CCF_inds);
            U_(i) = (U(CCF_inds).' * voxels_i) / sum(voxels_i);
            voxels_(i) = sum(voxels_i);
        end
    end
    U_(isnan(U_)) = (U_(~isnan(U_)).' * voxels_(~isnan(U_))) / ...
        sum(voxels_(~isnan(U_)));
    U = U_;
end

% Mean-normalize gene/cell-type data, works best empirically
U = U ./ mean(U,'omitmissing');
if logical(ipR.datapca_nexis_sv) && (length(ipR.datalist_nexis_sv) > 1)
    U_mean = mean(U,2,'omitmissing');
    [~, score, ~, ~, ~] = pca(U);
    U = score(:,1);
    if corr(U,U_mean) < 0
        U = -U;
    end
end

% Solve and store results
outputs.nexis_sv = struct;
%
%
% No bootstrapping
%
%
if ~logical(ipR.bootstrapping_nexis_sv)
    fprintf('Creating Optimal NexIS:SV Model\n');
    time_stamps = tpts.(ipR.study);
    if size(C,1) ~= size(data426.(ipR.study),1) % Convert all to CCF space for simulation/comparison
        pathology_raw = DataToCCF(data426.(ipR.study),ipR.study,ipR.matdir);
        pathology = normalizer(pathology_raw,ipR.normtype);
        pathology_orig = pathology;
        time_stamps_orig = time_stamps;
        if ~isnan(seed426.(ipR.study)) % Convert not NaN seed to CCF space for model init.
            seed_location = DataToCCF(seed426.(ipR.study),ipR.study,ipR.matdir);
        else % Use timepoint 1 pathology as init. for NaN seed (e.g., Hurtado et al.)
            seed_location = pathology(:,1);
            pathology = pathology(:,2:end);
            time_stamps = time_stamps(2:end) - time_stamps(1); % offset model times from t1
            ipR.param_init(1) = 1; ipR.ub(1) = 1; ipR.lb(1) = 1; % don't fit gamma
        end
    else % Use pathology and seed as is
        pathology_raw = data426.(ipR.study);
        pathology = normalizer(pathology_raw,ipR.normtype);
        if isnan(seed426.(ipR.study)) % Use timepoint 1 pathology as init. for NaN seed (i.e., run from baseline)
            seed_location = pathology(:,1);
            pathology = pathology(:,2:end);
            time_stamps = time_stamps(2:end) - time_stamps(1); % offset model times from t1
            ipR.param_init(1) = 1; ipR.ub(1) = 1; ipR.lb(1) = 1; % don't fit gamma
        else
            seed_location = seed426.(ipR.study);
        end
        pathology_orig = pathology;
        time_stamps_orig = time_stamps;
    end

    n_types = size(U,2);
    ndmflds = fieldnames(outputs.nexis_global);
    if length(ndmflds) == 1
        param_inits = outputs.nexis_global.Full.param_fit; 
    else
        param_inits = zeros((length(ndmflds)-1),length(outputs.nexis_global.Full.param_fit));
        for k = 1:(length(ndmflds)-1)
            fld = ndmflds{k};
            param_inits(k,:) = outputs.nexis_global.(fld).param_fit;
        end    
    end

    if strcmp(ipR.bounds_type_nexis_sv,'none') && (size(param_inits,1) > 1)
        param_init = median(param_inits); % Was originally mean, median probably better estimator
        ub = param_init; % fix all global parameters
        lb = param_init; % fix all global parameters
    elseif strcmp(ipR.bounds_type_nexis_sv,'unconstrained') && (size(param_inits,1) > 1)
        param_init = median(param_inits); % Was originally mean, median probably better estimator
        ub = [param_init(1),Inf,Inf,1,0,0]; % fix gamma, unconstrain others
        lb = [param_init(1),0,0,0,0,0]; % fix gamma, unconstrain others;
    elseif strcmp(ipR.bounds_type_nexis_sv,'old') && (size(param_inits,1) > 1)
        param_init = mean(param_inits); % Was originally mean, median probably better estimator
        ub = [1.3*param_init(1:4),0,0]; % Vary within +/- 30%
        lb = [0.7*param_init(1:4),0,0]; % Vary within +/- 30%
    elseif (size(param_inits,1) == 1)
        param_init = param_inits;
        ub = [param_init(1),10*param_init(2),10*param_init(3),1,0,0]; % fix gamma, very loosely constrain others;
        lb = [param_init(1),0.1*param_init(2),0.1*param_init(3),0,0,0]; % fix gamma, very loosely constrain others;
    else
        prct = str2double(ipR.bounds_type_nexis_sv(4:end));
        param_init = median(param_inits);
        ub = prctile(param_inits,((100-prct)/2)+prct,1); ub(1) = param_init(1); % use CI to bound all but gamma
        lb = prctile(param_inits,((100-prct)/2),1); lb(1) = param_init(1); % use CI to bound all but gamma
    end

    if ~logical(ipR.w_dir)
        param_init(4) = 0.5;
        ub(4) = 0.5;
        lb(4) = 0.5;
    end
    
    mordervec = zeros(1,length(ipR.param_init));
    for m = 1:length(param_init)
        inclparam = (lb(m) ~= ub(m));
        mordervec(m) = inclparam;
    end
    morder = 1 + sum(mordervec) + 2*n_types;
    
    param_init = [param_init(1:4),zeros(1,n_types),zeros(1,n_types)];
    lb = [lb(1:4),-Inf(1,n_types),-Inf(1,n_types)];
    ub = [ub(1:4),Inf(1,n_types),Inf(1,n_types)];
    objfun_handle = @(param) CostFunction_NexIS(param,C,U,time_stamps,...
        seed_location,pathology,ipR.solvetype,ipR.volcorrect,ipR.costfun,...
        ipR.excltpts_costfun,ipR.exclseed_costfun,ipR.use_dataspace,ipR.study,...
        ipR.logtrans,ipR.lambda,ipR.matdir);
    if logical(ipR.fmindisplay_nexis_sv)
        options = optimoptions(@fmincon,'Display','final-detailed','Algorithm',ipR.algo,...
            'MaxFunctionEvaluations',ipR.maxeval,'OptimalityTolerance',ipR.opttol,...
            'FunctionTolerance',ipR.fxntol,'StepTolerance',ipR.steptol);
    else
        options = optimoptions(@fmincon,'Algorithm',ipR.algo,...
            'MaxFunctionEvaluations',ipR.maxeval,'OptimalityTolerance',ipR.opttol,...
            'FunctionTolerance',ipR.fxntol,'StepTolerance',ipR.steptol);
    end
    [param_num, fval_num] = fmincon(objfun_handle,param_init,[],[],[],[],lb,ub,[],options);
    % tinds = setdiff(1:length(time_stamps), ipR.excltpts_costfun);
    ynum = NexIS_fun(C,U,time_stamps,seed_location,param_num,ipR.solvetype,ipR.volcorrect,ipR.matdir);
    
    % Store all outputs
    if ipR.use_dataspace && (size(C,1) ~= size(data426.(ipR.study),1)) % Convert back to data space
        pathology = CCFToData(pathology,ipR.study,ipR.matdir);
        % pathology = CCFToData(pathology(:,tinds),ipR.study,ipR.matdir);
        pathology_orig = CCFToData(pathology_orig,ipR.study,ipR.matdir); 
        ynum = CCFToData(ynum,ipR.study,ipR.matdir);
        % ynum = CCFToData(ynum(:,tinds),ipR.study,ipR.matdir);
        baseline = pathology_orig(:,1);
        if isnan(seed426.(ipR.study))
            seed_save = NaN;
        else
            seed_save = CCFToData(seed_location,ipR.study,ipR.matdir);
        end
    else
        if isnan(seed426.(ipR.study))
            seed_save = NaN;
        else
            seed_save = seed_location;
        end
        baseline = pathology_orig(:,1);
        % pathology = pathology(:,tinds);
        % ynum = ynum(:,tinds);
    end
    outputs.nexis_sv.Full.data = pathology; % this has been normalized
    outputs.nexis_sv.Full.baseline = baseline; % this has been normalized
    % if isnan(seed426.(ipR.study))
    %     outputs.nexis_sv.Full.time_stamps = time_stamps_orig(tinds+1);
    % else
    %     outputs.nexis_sv.Full.time_stamps = time_stamps_orig(tinds);
    % end
    outputs.nexis_global.Full.time_stamps = time_stamps_orig;
    outputs.nexis_sv.Full.predicted = ynum;
    outputs.nexis_sv.Full.param_fit = param_num;
    outputs.nexis_sv.Full.fval = fval_num;
    outputs.nexis_sv.Full.init.seed = seed_save;
    outputs.nexis_sv.Full.init.C = C;
    outputs.nexis_sv.Full.init.U_norm = U;
    outputs.nexis_sv.Full.init.study = ipR.study;
    outputs.nexis_sv.Full.init.solvetype = ipR.solvetype;
    outputs.nexis_sv.Full.init.volcorrect = ipR.volcorrect;
    outputs.nexis_sv.Full.init.normtype = ipR.normtype;
    outputs.nexis_sv.Full.init.costfun = ipR.costfun;
    outputs.nexis_sv.Full.init.lambda = ipR.lambda;
    outputs.nexis_sv.Full.init.exclseed_costfun = ipR.exclseed_costfun;
    outputs.nexis_sv.Full.init.excltpts_costfun = ipR.excltpts_costfun;
    outputs.nexis_sv.Full.init.logtrans = ipR.logtrans;
    outputs.nexis_sv.Full.init.w_dir = ipR.w_dir;
    outputs.nexis_sv.Full.init.param_init = param_init;
    outputs.nexis_sv.Full.init.ub = ub;
    outputs.nexis_sv.Full.init.lb = lb;
    outputs.nexis_sv.Full.init.bootstrapping_nexis_sv = ipR.bootstrapping_nexis_sv;
    outputs.nexis_sv.Full.init.resample_rate_nexis_sv = ipR.resample_rate_nexis_sv;
    outputs.nexis_sv.Full.init.bounds_type_nexis_sv = ipR.bounds_type_nexis_sv;
    outputs.nexis_sv.Full.init.datatype_nexis_sv = ipR.datatype_nexis_sv;
    outputs.nexis_sv.Full.init.datalist_nexis_sv = ipR.datalist_nexis_sv;
    outputs.nexis_sv.Full.init.datapca_nexis_sv = ipR.datapca_nexis_sv;
    outputs.nexis_sv.Full.init.niters_nexis_sv = ipR.niters_nexis_sv;
    outputs.nexis_sv.Full.fmincon.optimality_tolerance = ipR.opttol;
    outputs.nexis_sv.Full.fmincon.function_tolerance = ipR.fxntol;
    outputs.nexis_sv.Full.fmincon.step_tolerance = ipR.steptol;
    outputs.nexis_sv.Full.fmincon.algorithm = ipR.algo;
    outputs.nexis_sv.Full.fmincon.max_evaluations = ipR.maxeval;
    
    % Calculate per-timepoint R values
    ynum_fitassess = ynum;
    pathology_fitassess = pathology;
    if ipR.exclseed_outputs && all(~isnan(seed426.(ipR.study)))
        seedinds = find(seed_save);
        ynum_fitassess(seedinds,:) = [];
        pathology_fitassess(seedinds,:) = [];
    end

    Rvalues = zeros(1,length(time_stamps));
    for jj = 1:length(time_stamps)
        Rvalues(jj) = corr(ynum_fitassess(:,jj),pathology_fitassess(:,jj),'rows','complete');
    end
    outputs.nexis_sv.Full.results.Corrs = Rvalues; 

    % Rvalues = zeros(1,length(tinds));
    % for jj = 1:length(tinds)
    %     Rvalues(jj) = corr(ynum(:,jj),pathology(:,jj),'rows','complete');
    % end
    % outputs.nexis_sv.Full.results.Corrs = Rvalues; % NOTE: not corrected for seed

    P = reshape(pathology, [], 1);
    Y = reshape(ynum, [], 1);
    numObs1 = length(P(~isnan(P)));
    lm_nexis = fitlm(Y, P);
    logL = lm_nexis.LogLikelihood;
    outputs.nexis_sv.Full.results.lm_LogL = logL;
    outputs.nexis_sv.Full.results.lm_AIC = -2*logL + 2*morder;
    outputs.nexis_sv.Full.results.lm_BIC = -2*logL + log(numObs1)*morder;
    outputs.nexis_sv.Full.results.lm_intercept = lm_nexis.Coefficients.Estimate(1);
    outputs.nexis_sv.Full.results.lm_pval = lm_nexis.Coefficients.pValue(2);
    outputs.nexis_sv.Full.results.lm_Rsquared_ord = lm_nexis.Rsquared.Ordinary;
    outputs.nexis_sv.Full.results.lm_Rsquared_adj = lm_nexis.Rsquared.Adjusted;

    % flow = FlowCalculator(ynum,C,beta_num,0,U,b_num);
    % for i = 1:size(flow,3)
    %     flow_ = flow(:,:,i);
    %     flow_(flow_ < prctile(nonzeros(flow),ipR.flowthresh)) = 0;
    %     flow(:,:,i) = flow_;
    % end
    % outputs.nexis_sv.Full.flow = flow;

    if logical(ipR.verbose_nexis_sv)
        disp('--------------------------------------------------')
        disp('General nexis_sv minimizing quadratic error at all time stamps with fmincon')
        disp(' ')
        disp(['Optimal seed rescale value = ' num2str(param_num(1))])
        disp(['Optimal alpha = ' num2str(param_num(2))])
        disp(['Optimal beta = ' num2str(param_num(3))])
        if logical(ipR.w_dir)
            disp(['Optimal s = ' num2str(param_num(4))])
        end

        for i = 1:length(type_inds)
            disp(['Gene/Cell Type  ' num2str(i)])
            disp(['Optimal b = ' num2str(param_num(4 + n_types + i))])
            disp(['Optimal p = ' num2str(param_num(4 + 2*n_types + i))])
            disp(' ')
        end

        disp('R values at each time stamp')
        disp(Rvalues)
        disp(' ')
        if strcmp(ipR.costfun,'LinR')
            disp(['Cost Function = ' num2str(length(time_stamps)) ' - sum(LinR)'])
        else
            disp(ipR.costfun)
        end
        disp(fval_num)
        disp(['AIC = ' num2str(outputs.nexis_sv.Full.results.lm_AIC)])
        disp(['BIC = ' num2str(outputs.nexis_sv.Full.results.lm_BIC)])
        disp(['Intercept = ' num2str(outputs.nexis_sv.Full.results.lm_intercept)])
        disp(['pValue = ' num2str(outputs.nexis_sv.Full.results.lm_pval)])
        disp(['Rsqr_ord = ' num2str(outputs.nexis_sv.Full.results.lm_Rsquared_ord)])
        disp(['Rsqr_adj = ' num2str(outputs.nexis_sv.Full.results.lm_Rsquared_adj)])
        disp(' ')
    end
%
%
% With bootstrapping
%
%
else
    rng(1); % Different from NexIS_global.m seed
    for i = 1:ipR.niters_nexis_sv
        fldname = sprintf('Iter_%d',i);
        time_stamps = tpts.(ipR.study);
        if size(C,1) ~= size(data426.(ipR.study),1) % Convert all to CCF space for simulation/comparison
            pathology_raw = data426.(ipR.study);
            notnaninds = find(~isnan(pathology_raw(:,1)));
            settonansize = round((1-ipR.resample_rate)*length(notnaninds));
            settonaninds = randperm(length(notnaninds));
            settonaninds = notnaninds(settonaninds(1:settonansize));
            pathology = pathology_raw;
            pathology_raw(settonaninds,:) = NaN;
            pathology_raw = DataToCCF(pathology_raw,ipR.study,ipR.matdir);
            pathology = DataToCCF(pathology,ipR.study,ipR.matdir);
            pathology = normalizer(pathology,ipR.normtype);
            pathology_orig = pathology;
            pathology(isnan(pathology_raw(:,1)),:) = NaN;
            time_stamps_orig = time_stamps;
            if all(~isnan(seed426.(ipR.study))) % Convert not NaN seed to CCF space for model init.
                seed_location = DataToCCF(seed426.(ipR.study),ipR.study,ipR.matdir);
            else % Use timepoint 1 pathology as init. for NaN seed (e.g., Hurtado et al.)
                seed_location = pathology(:,1);
                pathology = pathology(:,2:end);
                % pathology_orig = pathology_orig(:,2:end);
                time_stamps = time_stamps(2:end) - time_stamps(1); % offset model times from t1
                ipR.param_init(1) = 1; ipR.ub(1) = 1; ipR.lb(1) = 1; % don't fit gamma
            end
        else % Use pathology and seed as is
            pathology_raw = data426.(ipR.study);
            notnaninds = find(~isnan(pathology_raw(:,1)));
            settonansize = round((1-ipR.resample_rate)*length(notnaninds));
            settonaninds = randperm(length(notnaninds));
            settonaninds = notnaninds(settonaninds(1:settonansize));
            pathology = pathology_raw;
            pathology_raw(settonaninds,:) = NaN;
            pathology = normalizer(pathology,ipR.normtype);
            pathology_orig = pathology;
            pathology(isnan(pathology_raw(:,1)),:) = NaN;
            time_stamps_orig = time_stamps;
            if isnan(seed426.(ipR.study)) % Use timepoint 1 pathology as init. for NaN seed (i.e., run from baseline)
                seed_location = pathology(:,1);
                pathology = pathology(:,2:end);
                time_stamps = time_stamps(2:end) - time_stamps(1); % offset model times from t1
                ipR.param_init(1) = 1; ipR.ub(1) = 1; ipR.lb(1) = 1; % don't fit gamma
            else
                seed_location = seed426.(ipR.study);
            end
        end

        fprintf('NexIS:SV Bootstrapping Iteration %d/%d\n',i,ipR.niters_nexis_sv);        
        n_types = size(U,2);
        ndmflds = fieldnames(outputs.nexis_sv);
        if length(ndmflds) == 1
            param_inits = outputs.nexis_sv.Full.param_fit; 
        else
            param_inits = zeros((length(ndmflds)-1),length(outputs.nexis_sv.Full.param_fit));
            for k = 1:(length(ndmflds)-1)
                fld = ndmflds{k};
                param_inits(k,:) = outputs.nexis_sv.(fld).param_fit;
            end    
        end

        if strcmp(ipR.bounds_type_nexis_sv,'none') && (size(param_inits,1) > 1)
            param_init = median(param_inits); % Was originally mean, median probably better estimator
            ub = param_init; % fix all global parameters
            lb = param_init; % fix all global parameters
        elseif strcmp(ipR.bounds_type_nexis_sv,'unconstrained') && (size(param_inits,1) > 1)
            param_init = median(param_inits); % Was originally mean, median probably better estimator
            ub = [param_init(1),Inf,Inf,1,0,0]; % fix gamma, unconstrain others
            lb = [param_init(1),0,0,0,0,0]; % fix gamma, unconstrain others;
        elseif strcmp(ipR.bounds_type_nexis_sv,'old') && (size(param_inits,1) > 1)
            param_init = mean(param_inits); % Was originally mean, median probably better estimator
            ub = [1.3*param_init(1:4),0,0]; % Vary within +/- 30%
            lb = [0.7*param_init(1:4),0,0]; % Vary within +/- 30%
        elseif (size(param_inits,1) == 1)
            param_init = param_inits;
            ub = [param_init(1),10*param_init(2),10*param_init(3),1,0,0]; % fix gamma, very loosely constrain others;
            lb = [param_init(1),0.1*param_init(2),0.1*param_init(3),0,0,0]; % fix gamma, very loosely constrain others;
        else
            prct = str2double(ipR.bounds_type_nexis_sv(4:end));
            param_init = median(param_inits);
            ub = prctile(param_inits,((100-prct)/2)+prct,1); ub(1) = param_init(1); % use CI to bound all but gamma
            lb = prctile(param_inits,((100-prct)/2),1); lb(1) = param_init(1); % use CI to bound all but gamma
        end

        if ~logical(ipR.w_dir)
            param_init(4) = 0.5;
            ub(4) = 0.5;
            lb(4) = 0.5;
        end
        
        mordervec = zeros(1,length(ipR.param_init));
        for m = 1:length(param_init)
            inclparam = (lb(m) ~= ub(m));
            mordervec(m) = inclparam;
        end
        morder = 1 + sum(mordervec) + 2*n_types;
        
        param_init = [param_init(1:4),zeros(1,n_types),zeros(1,n_types)];
        lb = [lb(1:4),-Inf(1,n_types),-Inf(1,n_types)];
        ub = [ub(1:4),Inf(1,n_types),Inf(1,n_types)];

        objfun_handle = @(param) CostFunction_NexIS(param,C,U,time_stamps,...
            seed_location,pathology,ipR.solvetype,ipR.volcorrect,ipR.costfun,...
            ipR.excltpts_costfun,ipR.exclseed_costfun,ipR.use_dataspace,ipR.study,...
            ipR.logtrans,ipR.lambda,ipR.matdir);
        if logical(ipR.fmindisplay)
            options = optimoptions(@fmincon,'Display','final-detailed','Algorithm',ipR.algo,...
                'MaxFunctionEvaluations',ipR.maxeval,'OptimalityTolerance',ipR.opttol,...
                'FunctionTolerance',ipR.fxntol,'StepTolerance',ipR.steptol);
        else
            options = optimoptions(@fmincon,'Algorithm',ipR.algo,...
                'MaxFunctionEvaluations',ipR.maxeval,'OptimalityTolerance',ipR.opttol,...
                'FunctionTolerance',ipR.fxntol,'StepTolerance',ipR.steptol);
        end
        try 
            [param_num, fval_num] = fmincon(objfun_handle,param_init,[],[],[],[],lb,ub,[],options);
        catch ME
            fprintf('Error: %s\n',ME.message);
            fprintf('Trying initial parameter tweak\n')
            try 
                param_init = lb + rand(1,length(param_init)).*(ub - lb);
                param_init(isnan(param_init)) = 0;
                [param_num, fval_num] = fmincon(objfun_handle,param_init,[],[],[],[],lb,ub,[],options);
            catch ME
                fprintf('Error: %s\n',ME.message);
                fprintf('Skipping iteration\n');
                continue;
            end
        end

        % Solve NexIS SV with the optimized parameters
        ynum = NexIS_fun(C,U,time_stamps,seed_location,param_num,ipR.solvetype,ipR.volcorrect,ipR.matdir);

        % Store all outputs
        if ipR.use_dataspace && (size(C,1) ~= size(data426.(ipR.study),1)) % Convert back to data space
            pathology = CCFToData(pathology,ipR.study,ipR.matdir);
            pathology_orig = CCFToData(pathology_orig,ipR.study,ipR.matdir); 
            ynum = CCFToData(ynum,ipR.study,ipR.matdir);
            baseline = pathology_orig(:,1);
            ynum_save = ynum;
            if isnan(seed426.(ipR.study))
                seed_save = NaN;
                % ynum_save = [pathology_orig(:,1), ynum]; % save ynum with baseline to keep consistent with pathology_orig
            else
                seed_save = CCFToData(seed_location,ipR.study,ipR.matdir);
                % ynum_save = ynum;
            end
        else
            baseline = pathology_orig(:,1);
            ynum_save = ynum;
            if isnan(seed426.(ipR.study))
                seed_save = NaN;
                % ynum_save = [pathology_orig(:,1), ynum]; % save ynum with baseline to keep consistent with pathology_orig
            else
                seed_save = seed_location;
                % ynum_save = ynum;
            end
        end

        outputs.nexis_sv.(fldname).data = pathology;
        outputs.nexis_sv.(fldname).baseline = baseline;
        outputs.nexis_sv.(fldname).time_stamps = time_stamps_orig;
        outputs.nexis_sv.(fldname).predicted = ynum_save;
        outputs.nexis_sv.(fldname).param_fit = param_num;
        outputs.nexis_sv.(fldname).fval = fval_num;
        outputs.nexis_sv.(fldname).init.seed = seed_save;
        outputs.nexis_sv.(fldname).init.C = C;
        outputs.nexis_sv.(fldname).init.study = ipR.study;
        outputs.nexis_sv.(fldname).init.solvetype = ipR.solvetype;
        outputs.nexis_sv.(fldname).init.volcorrect = ipR.volcorrect;
        outputs.nexis_sv.(fldname).init.normtype = ipR.normtype;
        outputs.nexis_sv.(fldname).init.costfun = ipR.costfun;
        outputs.nexis_sv.(fldname).init.exclseed_costfun = ipR.exclseed_costfun;
        outputs.nexis_sv.(fldname).init.excltpts_costfun = ipR.excltpts_costfun;
        outputs.nexis_sv.(fldname).init.w_dir = ipR.w_dir;
        outputs.nexis_sv.(fldname).init.param_init = param_init;
        outputs.nexis_sv.(fldname).init.ub = ub;
        outputs.nexis_sv.(fldname).init.lb = lb;
        outputs.nexis_sv.(fldname).init.bootstrapping_nexis_sv = ipR.bootstrapping_nexis_sv;
        outputs.nexis_sv.(fldname).init.resample_rate_nexis_sv = ipR.resample_rate_nexis_sv;
        outputs.nexis_sv.(fldname).init.bounds_type_nexis_sv = ipR.bounds_type_nexis_sv;
        outputs.nexis_sv.(fldname).init.datatype_nexis_sv = ipR.datatype_nexis_sv;
        outputs.nexis_sv.(fldname).init.datalist_nexis_sv = ipR.datalist_nexis_sv;
        outputs.nexis_sv.(fldname).init.datapca_nexis_sv = ipR.datapca_nexis_sv;
        outputs.nexis_sv.(fldname).init.niters_nexis_sv = ipR.niters_nexis_sv;
        outputs.nexis_sv.(fldname).fmincon.optimality_tolerance = ipR.opttol;
        outputs.nexis_sv.(fldname).fmincon.function_tolerance = ipR.fxntol;
        outputs.nexis_sv.(fldname).fmincon.step_tolerance = ipR.steptol;
        outputs.nexis_sv.(fldname).fmincon.algorithm = ipR.algo;
        outputs.nexis_sv.(fldname).fmincon.max_evaluations = ipR.maxeval;

        Rvalues = zeros(1,length(time_stamps));
        for jj = 1:length(time_stamps)
            % if strcmp(ipR.corrtype,'R')
            Rvalues(jj) = corr(ynum(:,jj),pathology(:,jj), 'rows','complete');
            % elseif strcmp(ipR.corrtype,'R_c')
            %    naninds = isnan(pathology(:,1));
            %    newxt = ynum; newxt(naninds,:) = [];
            %    newpath = pathology; newpath(naninds,:) = [];
            %    Rvalues(jj) = LinRcalc(newxt(:,jj),newpath(:,jj));
            % end
        end
        outputs.nexis_sv.(fldname).results.Corrs = Rvalues;

        P = reshape(pathology, [], 1);
        Y = reshape(ynum, [], 1);
        numObs1 = length(P(~isnan(P)));
        lm_nexis_sv = fitlm(Y, P);
        logL = lm_nexis_sv.LogLikelihood;
        outputs.nexis_sv.(fldname).results.lm_LogL = logL;
        outputs.nexis_sv.(fldname).results.lm_AIC = -2*logL + 2*morder;
        outputs.nexis_sv.(fldname).results.lm_BIC = -2*logL + log(numObs1)*morder;
        outputs.nexis_sv.(fldname).results.lm_intercept = lm_nexis_sv.Coefficients.Estimate(1);
        outputs.nexis_sv.(fldname).results.lm_pval = lm_nexis_sv.Coefficients.pValue(1);
        outputs.nexis_sv.(fldname).results.lm_Rsquared_ord = lm_nexis_sv.Rsquared.Ordinary;
        outputs.nexis_sv.(fldname).results.lm_Rsquared_adj = lm_nexis_sv.Rsquared.Adjusted;
        if logical(ipR.verbose_nexis_sv)
            disp('--------------------------------------------------')
            disp('General nexis_sv minimizing quadratic error at all time stamps with fmincon')
            disp(' ')
            disp(['Fit seed rescale value = ' num2str(param_num(1))])
            disp(['Fit alpha = ' num2str(param_num(2))])
            disp(['Fit beta = ' num2str(param_num(3))])
            if logical(ipR.w_dir)
                disp(['Fit s = ' num2str(param_num(4))])
            end

            for k = 1:length(type_inds)
                disp(['Gene/Cell Type  ' num2str(k)])
                disp(['Fit b = ' num2str(param_num(4 + n_types + k))])
                disp(['Fit p = ' num2str(param_num(4 + 2*n_types + k))])
                disp(' ')
            end

            disp('R values at each time stamp')
            disp(Rvalues)
            disp(' ')
            if strcmp(ipR.costfun,'LinR')
                disp(['Cost Function = ' num2str(length(time_stamps)) ' - sum(LinR)'])
            else
                disp(ipR.costfun)
            end
            disp(fval_num)
            disp(['AIC = ' num2str(outputs.nexis_sv.(fldname).results.lm_AIC)])
            disp(['BIC = ' num2str(outputs.nexis_sv.(fldname).results.lm_BIC)])
            disp(['Intercept = ' num2str(outputs.nexis_sv.(fldname).results.lm_intercept)])
            disp(['pValue = ' num2str(outputs.nexis_sv.(fldname).results.lm_pval)])
            disp(['Rsqr_ord = ' num2str(outputs.nexis_sv.(fldname).results.lm_Rsquared_ord)])
            disp(['Rsqr_adj = ' num2str(outputs.nexis_sv.(fldname).results.lm_Rsquared_adj)])
            disp(' ')
        end
        clear ub lb param_init
    end
    %
    fprintf('Creating Optimal NexIS:SV Model\n');
    time_stamps = tpts.(ipR.study);
    if size(C,1) ~= size(data426.(ipR.study),1) % Convert all to CCF space for simulation/comparison
        pathology_raw = DataToCCF(data426.(ipR.study),ipR.study,ipR.matdir);
        pathology = normalizer(pathology_raw,ipR.normtype);
        pathology_orig = pathology;
        time_stamps_orig = time_stamps;
        if all(~isnan(seed426.(ipR.study))) % Convert not NaN seed to CCF space for model init.
            seed_location = DataToCCF(seed426.(ipR.study),ipR.study,ipR.matdir);
        else % Use timepoint 1 pathology as init. for NaN seed (e.g., Hurtado et al.)
            seed_location = pathology(:,1);
            pathology = pathology(:,2:end);
            time_stamps = time_stamps(2:end) - time_stamps(1); % offset model times from t1
            ipR.param_init(1) = 1; ipR.ub(1) = 1; ipR.lb(1) = 1; % don't fit gamma
        end
    else % Use pathology and seed as is
        pathology_raw = data426.(ipR.study);
        pathology = normalizer(pathology_raw,ipR.normtype);
        time_stamps_orig = time_stamps;
        if isnan(seed426.(ipR.study)) % Use timepoint 1 pathology as init. for NaN seed (i.e., running from baseline)
            seed_location = pathology(:,1);
            pathology = pathology(:,2:end);
            time_stamps = time_stamps(2:end) - time_stamps(1); % offset model times from t1
            ipR.param_init(1) = 1; ipR.ub(1) = 1; ipR.lb(1) = 1; % don't fit gamma
        else
            seed_location = seed426.(ipR.study);
        end
        pathology_orig = pathology;
    end

    if any(isnan(seed_location))
        seed_location(isnan(seed_location)) = 0;
    end

    fldnames = fieldnames(outputs.nexis_global);
    param_fits = zeros(length(fldnames),length(outputs.nexis_global.(fldnames{1}).param_fit));
    for i = 1:length(fldnames)
        fldname = fldnames{i};
        param_fits(i,:) = outputs.nexis_global.(fldname).param_fit;
    end
    
    % Evaluate NexIS:global with best estimate of parameters
    param_opt = median(param_fits); % Used mean before, median should be better estimator
    yopt = NexIS_fun(C,U,time_stamps,seed_location,param_opt,ipR.solvetype,ipR.volcorrect,ipR.matdir);

    % Store all outputs
    if ipR.use_dataspace && (size(C,1) ~= size(data426.(ipR.study),1)) 
        pathology_fvalcalc = pathology; % need pathology in CCF space for cost function input
        pathology = CCFToData(pathology,ipR.study,ipR.matdir);% Convert back to data space
        pathology_orig = CCFToData(pathology_orig,ipR.study,ipR.matdir); 
        yopt = CCFToData(yopt,ipR.study,ipR.matdir);
        yopt_save = yopt;
        if isnan(seed426.(ipR.study))
            seed_save = NaN;
            % yopt_save = [pathology_orig(:,1), yopt]; % save yopt with baseline to keep consistent with pathology_orig
        else
            seed_save = CCFToData(seed_location,ipR.study,ipR.matdir);
            % yopt_save = yopt;
        end
    else
        pathology_fvalcalc = pathology; % need pathology in CCF space for cost function input
        seed_save = seed_location;
        yopt_save = yopt;
    end

    outputs.nexis_sv.Full.data = pathology;
    outputs.nexis_sv.Full.baseline = pathology_orig(:,1);
    outputs.nexis_sv.Full.time_stamps = time_stamps_orig;
    outputs.nexis_sv.Full.predicted = yopt_save;
    outputs.nexis_sv.Full.param_fit = param_opt;
    outputs.nexis_sv.Full.fval = CostFunction_NexIS(param_opt,...
                                                        C,...
                                                        U,...
                                                        time_stamps,...
                                                        seed_location,...
                                                        pathology_fvalcalc,...
                                                        ipR.solvetype,...
                                                        ipR.volcorrect,...
                                                        ipR.costfun,...
                                                        ipR.excltpts_costfun,...
                                                        ipR.exclseed_costfun,...
                                                        ipR.use_dataspace,...
                                                        ipR.study,...
                                                        ipR.logtrans,...
                                                        ipR.lambda,...
                                                        ipR.matdir);
    outputs.nexis_sv.Full.init = outputs.nexis_sv.(fldnames{1}).init;
    outputs.nexis_sv.Full.init.seed = seed_save;
    outputs.nexis_sv.Full.fmincon = [];
    Rvalues = zeros(1,length(time_stamps));
    for jj = 1:length(time_stamps)
        % if strcmp(ipR.corrtype,'R')
        Rvalues(jj) = corr(yopt(:,jj),pathology(:,jj), 'rows','complete');
        % elseif strcmp(ipR.corrtype,'R_c')
        %    naninds = isnan(pathology(:,1));
        %    newxt = ynum; newxt(naninds,:) = [];
        %    newpath = pathology; newpath(naninds,:) = [];
        %    Rvalues(jj) = LinRcalc(newxt(:,jj),newpath(:,jj));
        % end
    end
    outputs.nexis_sv.Full.results.Corrs = Rvalues;
    P = reshape(pathology, [], 1);
    Y = reshape(yopt, [], 1);
    numObs1 = length(P(~isnan(P)));
    lm_nexis_sv = fitlm(Y, P);
    logL = lm_nexis_sv.LogLikelihood;
    outputs.nexis_sv.Full.results.lm_LogL = logL;
    outputs.nexis_sv.Full.results.lm_AIC = -2*logL + 2*morder;
    outputs.nexis_sv.Full.results.lm_BIC = -2*logL + log(numObs1)*morder;
    outputs.nexis_sv.Full.results.lm_intercept = lm_nexis_sv.Coefficients.Estimate(1);
    outputs.nexis_sv.Full.results.lm_pval = lm_nexis_sv.Coefficients.pValue(1);
    outputs.nexis_sv.Full.results.lm_Rsquared_ord = lm_nexis_sv.Rsquared.Ordinary;
    outputs.nexis_sv.Full.results.lm_Rsquared_adj = lm_nexis_sv.Rsquared.Adjusted;
    
    % flow = FlowCalculator(yopt,C,beta_opt,0,U,b_opt);
    % for i = 1:size(flow,3)
    %     flow_ = flow(:,:,i);
    %     flow_(flow_ < prctile(nonzeros(flow),ipR.flowthresh)) = 0;
    %     flow(:,:,i) = flow_;
    % end
    % outputs.nexis_sv.Full.flow = flow;
    
    if logical(ipR.verbose_nexis_sv)
        disp('--------------------------------------------------')
        disp('Optimal model from bootstrapping')
        disp(' ')
        disp(['Optimal seed rescale value = ' num2str(param_opt(1))])
        disp(['Optimal alpha = ' num2str(param_opt(2))])
        disp(['Optimal beta = ' num2str(param_opt(3))])
        if logical(ipR.w_dir)
            disp(['Optimal s = ' num2str(param_opt(4))])
        end
        disp(' ')

        for k = 1:length(type_inds)
            disp(['Gene/Cell Type  ' num2str(k)])
            disp(['Optimal b = ' num2str(param_opt(4 + n_types + k))])
            disp(['Optimal p = ' num2str(param_opt(4 + 2*n_types + k))])
            disp(' ')
        end

        disp('R values at each time stamp')
        disp(Rvalues)
        disp(' ')

        disp(['AIC = ' num2str(outputs.nexis_sv.Full.results.lm_AIC)])
        disp(['BIC = ' num2str(outputs.nexis_sv.Full.results.lm_BIC)])
        disp(['Intercept = ' num2str(outputs.nexis_sv.Full.results.lm_intercept)])
        disp(['pValue = ' num2str(outputs.nexis_sv.Full.results.lm_pval)])
        disp(['Rsqr_ord = ' num2str(outputs.nexis_sv.Full.results.lm_Rsquared_ord)])
        disp(['Rsqr_adj = ' num2str(outputs.nexis_sv.Full.results.lm_Rsquared_adj)])
        disp(' ')
    end
end

    function indices = NameIndex(names,dattypenexis_sv)
        if strcmp(dattypenexis_sv,'gene')
            datstruct = load([cd filesep 'raw_data_mouse' filesep 'GeneExpressionMaps.mat'],'GeneExpressionMaps');
            namescell = datstruct.GeneExpressionMaps.All.gene_names;
        else
            datstruct = load([cd filesep 'raw_data_mouse' filesep 'CellTypeMaps.mat'],'CellTypeMaps');
            namescell = datstruct.CellTypeMaps.(dattypenexis_sv).classkey;
        % elseif strcmp(dattypenexis_sv,'ct_tasic')
        %     load([cd filesep 'raw_data_mouse' filesep 'classkey_tasic.mat'],'classkey_tasic');
        %     namescell = classkey_tasic;
        % elseif strcmp(dattypenexis_sv,'ct_zeisel')
        %     load([cd filesep 'raw_data_mouse' filesep 'classkey_zeisel.mat'],'classkey_zeisel');
        %     namescell = classkey_zeisel;
        % elseif strcmp(dattypenexis_sv,'ct_yao')
        %     load([ipR.yaoctdir filesep 'Yao_Dependencies.mat'],'classkey');
        %     namescell = classkey;
        end
        indices = zeros(1,length(names));
        for n_i = 1:length(names)
            indices(n_i) = find(ismember(namescell,names(n_i)));
        end
    end

    function normdata = normalizer(data,ntype)
        if strcmp(ntype,'sum')
            normdata = data/sum(data(:,1),'omitnan');
        elseif strcmp(ntype,'masssum')
            load([ipR.matdir filesep 'DefaultAtlas.mat'], 'DefaultAtlas');
            voxels_2hem = DefaultAtlas.volumes;
            data_1 = data(:,1);
            nonnans = isnan(data_1);
            data_1(nonnans) = [];
            voxels_2hem(nonnans) = [];
            masssum = data_1.' * voxels_2hem;
            normdata = data / masssum;
        elseif strcmp(ntype,'mean')
            normdata = data/mean(data(:,1),'omitnan');
        elseif strcmp(ntype,'norm2')
            normdata = data/norm(data(:,1),2);
        elseif strcmp(ntype,'log')
            normdata = log(data+1);
        else
            normdata = data; 
        end
    end

end