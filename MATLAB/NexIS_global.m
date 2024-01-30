function outputs = NexIS_global(varargin)
% Parameters to fit
% param(1) = gamma (seed rescale factor)
% param(2) = alpha
% param(3) = beta
% param(4) = s

% Define defaults and set inputs
study_ = 'User-specified'; 
C_ = [];
data_ = [];
tpts_ = [];
seed_ = [];
use_dataspace_ = 0;
matdir_ = [cd filesep 'raw_data_mouse'];

costfun_ = 'linr'; % see /lib_NexIS/CostFunction_NexIS.m for options
solvetype_ = 'analytic';
volcorrect_ = 1;
exclseed_costfun_ = 0;
excltpts_costfun_ = [];
logtrans_ = 'none';
normtype_ = 'sum';
Cnormtype_ = 'minmax';
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
flowthresh_ = 99.93;

ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
validBoolean = @(x) isscalar(x) && (x==0 || x==1);
validChar = @(x) ischar(x);
validST = @(x) ismember(x,{'analytic','numeric'});
validParam = @(x) (length(x) == 4);

addParameter(ip, 'study', study_, validChar);
addParameter(ip, 'C', C_);
addParameter(ip, 'data',data_);
addParameter(ip, 'tpts',tpts_);
addParameter(ip, 'seed', seed_);
addParameter(ip, 'matdir', matdir_);
addParameter(ip, 'use_dataspace', use_dataspace_, validBoolean);
addParameter(ip, 'costfun', costfun_, validChar);
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
addParameter(ip, 'flowthresh', flowthresh_, validScalar);
parse(ip, varargin{:});
ipR = ip.Results;

% Load in data from NexIS/raw_data_mouse directory if needed
if isempty(ipR.C) && ~strcmp(ipR.study,'User-specified')
    if (length(ipR.study) > 3) && strcmp(ipR.study(1:4),'asyn')
        load([ipR.matdir filesep 'mouse_aSynData_426.mat'],...
            'data426','seed426','tpts');
        ipR.study = ipR.study(6:end);
        load([ipR.matdir filesep 'Connectomes.mat'],'Connectomes');
        C = Connectomes.default;
    elseif strcmp(ipR.study,'Henderson')
        load([ipR.matdir filesep 'Henderson_Asyn_Data.mat'],...
            'tpts','Henderson_Asyn_Seed_Data','Henderson_Asyn_Pathology_Data');
        load([ipR.matdir filesep 'Henderson_Asyn_Data.mat'],...
            'Connection');    
        C = Connection;
        tpts_ = struct; tpts_.(ipR.study) = tpts.NTG; tpts = tpts_;
        seed426 = struct; seed426.(ipR.study) = Henderson_Asyn_Seed_Data;
        data426 = struct; data426.(ipR.study) = Henderson_Asyn_Pathology_Data.NTG;
    elseif strcmp(ipR.study(1:3),'GCI') || strcmp(ipR.study(1:3),'PFF')
        load([ipR.matdir filesep 'GCI_PFF_Data.mat'],...
            'tpts','GCI_PFF_Pathology_Data','GCI_PFF_Seed_Data');
        if strcmp(ipR.study,'GCI_Average_New')
            tpts.(ipR.study) = tpts.GCI_Average;
        end
        seed426 = struct; seed426.(ipR.study) = GCI_PFF_Seed_Data;
        data426 = struct; data426.(ipR.study) = GCI_PFF_Pathology_Data.(ipR.study);
    else
        load([ipR.matdir filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],...
            'mousedata_struct');
        load([ipR.matdir filesep 'Connectomes.mat'],'Connectomes');
        data426 = struct; data426.(ipR.study) = mousedata_struct.(ipR.study).data;
        seed426 = struct; seed426.(ipR.study) = mousedata_struct.(ipR.study).seed;
        tpts = struct; tpts.(ipR.study) = mousedata_struct.(ipR.study).time_stamps;
        C = Connectomes.default;
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

% Solve and store results
outputs.nexis_global = struct;
if ~logical(ipR.bootstrapping) % No bootstrapping of parameters
    fprintf('Creating Optimal NexIS:global Model\n');
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
        seed_location = seed426.(ipR.study);
        pathology = normalizer(pathology_raw,ipR.normtype);
        pathology_orig = pathology;
        time_stamps_orig = time_stamps;
    end
    U = zeros(size(C,1),1);
    if any(isnan(seed_location))
        seed_location(isnan(seed_location)) = 0;
    end
    if isnan(ipR.param_init(1))
        ipR.param_init(1) = sum(pathology(:,1),'omitnan')/nnz(seed_location); % heuristic default, study-dependent
    end

    if ~logical(ipR.w_dir)
        ipR.param_init(4) = 0.5;
        ipR.ub(4) = 0.5;
        ipR.lb(4) = 0.5;
    end
    mordervec = zeros(1,length(ipR.param_init));
    for m = 1:length(ipR.param_init)
        inclparam = (ipR.lb(m) ~= ipR.ub(m));
        mordervec(m) = inclparam;
    end
    morder = 1 + sum(mordervec);
    
    param_init = [ipR.param_init,zeros(1,2)]; % dummy values for b and p
    lb = [ipR.lb,zeros(1,2)]; % dummy values for b and p
    ub = [ipR.ub,zeros(1,2)]; % dummy values for b and p
    
    objfun_handle = @(param) CostFunction_NexIS(param,C,U,time_stamps,...
        seed_location,pathology,ipR.solvetype,ipR.volcorrect,ipR.costfun,...
        ipR.excltpts_costfun,ipR.exclseed_costfun,ipR.use_dataspace,ipR.study,...
        ipR.logtrans,ipR.matdir);
    if logical(ipR.fmindisplay)
        options = optimoptions(@fmincon,'Display','final-detailed','Algorithm',ipR.algo,...
            'MaxFunctionEvaluations',ipR.maxeval,'OptimalityTolerance',ipR.opttol,...
            'FunctionTolerance',ipR.fxntol,'StepTolerance',ipR.steptol);
    else
        options = optimoptions(@fmincon,'Algorithm',ipR.algo,...
            'MaxFunctionEvaluations',ipR.maxeval,'OptimalityTolerance',ipR.opttol,...
            'FunctionTolerance',ipR.fxntol,'StepTolerance',ipR.steptol);
    end
    [param_num, fval_num] = fmincon(objfun_handle,param_init,[],[],[],[],lb,ub,[],options);

    % Solve NexIS global with the optimal parameters
    ynum = NexIS_fun(C,U,time_stamps,seed_location,param_num,ipR.solvetype,ipR.volcorrect,ipR.matdir);
    
    % Store all outputs
    if ipR.use_dataspace && (size(C,1) ~= size(data426.(ipR.study),1)) % Convert back to data space
        pathology = CCFToData(pathology,ipR.study,ipR.matdir);
        pathology_orig = CCFToData(pathology_orig,ipR.study,ipR.matdir); 
        ynum = CCFToData(ynum,ipR.study,ipR.matdir);
        if isnan(seed426.(ipR.study))
            seed_save = NaN;
            ynum_save = [pathology_orig(:,1), ynum]; % save ynum with baseline to keep consistent with pathology_orig
        else
            seed_save = CCFToData(seed_location,ipR.study,ipR.matdir);
            ynum_save = ynum;
        end
    else
        seed_save = seed_location;
        ynum_save = ynum;
    end
    outputs.nexis_global.Full.data = pathology_orig; % this has been normalized
    outputs.nexis_global.Full.time_stamps = time_stamps_orig;
    outputs.nexis_global.Full.predicted = ynum_save;
    outputs.nexis_global.Full.param_fit = param_num;
    outputs.nexis_global.Full.fval = fval_num;
    outputs.nexis_global.Full.init.seed = seed_save;
    outputs.nexis_global.Full.init.C = C;
    if ismember(ipR.study,{'human','mouse'})
        outputs.nexis_global.Full.init.study = ['asyn ' ipR.study];
    else
        outputs.nexis_global.Full.init.study = ipR.study;
    end
    outputs.nexis_global.Full.init.solvetype = ipR.solvetype;
    outputs.nexis_global.Full.init.volcorrect = ipR.volcorrect;
    outputs.nexis_global.Full.init.normtype = ipR.normtype;
    outputs.nexis_global.Full.init.costfun = ipR.costfun;
    outputs.nexis_global.Full.init.exclseed_costfun = ipR.exclseed_costfun;
    outputs.nexis_global.Full.init.excltpts_costfun = ipR.excltpts_costfun;
    outputs.nexis_global.Full.init.w_dir = ipR.w_dir;
    outputs.nexis_global.Full.init.param_init = ipR.param_init;
    outputs.nexis_global.Full.init.ub = ipR.ub;
    outputs.nexis_global.Full.init.lb = ipR.lb;
    outputs.nexis_global.Full.init.bootstrapping = ipR.bootstrapping;
    outputs.nexis_global.Full.init.resample_rate = ipR.resample_rate;
    outputs.nexis_global.Full.init.niters = ipR.niters;
    outputs.nexis_global.Full.fmincon.optimality_tolerance = ipR.opttol;
    outputs.nexis_global.Full.fmincon.function_tolerance = ipR.fxntol;
    outputs.nexis_global.Full.fmincon.step_tolerance = ipR.steptol;
    outputs.nexis_global.Full.fmincon.algorithm = ipR.algo;
    outputs.nexis_global.Full.fmincon.max_evaluations = ipR.maxeval;
    
    % Calculate per-timepoint R values
    Rvalues = zeros(1,length(time_stamps));
    for jj = 1:length(time_stamps)
        Rvalues(jj) = corr(ynum(:,jj),pathology(:,jj),'rows','complete');
    end
    outputs.nexis_global.Full.results.Corrs = Rvalues;

    P = reshape(pathology, [], 1);
    Y = reshape(ynum, [], 1);
    numObs1 = length(P(~isnan(P)));
    lm_nexis = fitlm(Y, P);
    logL = lm_nexis.LogLikelihood;
    outputs.nexis_global.Full.results.lm_LogL = logL;
    outputs.nexis_global.Full.results.lm_AIC = -2*logL + 2*morder;
    outputs.nexis_global.Full.results.lm_BIC = -2*logL + log(numObs1)*morder;
    outputs.nexis_global.Full.results.lm_intercept = lm_nexis.Coefficients.Estimate(1);
    outputs.nexis_global.Full.results.lm_pval = lm_nexis.Coefficients.pValue(2);
    outputs.nexis_global.Full.results.lm_Rsquared_ord = lm_nexis.Rsquared.Ordinary;
    outputs.nexis_global.Full.results.lm_Rsquared_adj = lm_nexis.Rsquared.Adjusted;
    
    % flow = FlowCalculator(ynum,C,param_num(2),1,U,param_num(5));
    % for i = 1:size(flow,3)
    %     flow_ = flow(:,:,i);
    %     flow_(flow_ < prctile(nonzeros(flow),ipR.flowthresh)) = 0;
    %     flow(:,:,i) = flow_;
    % end
    % outputs.nexis_global.Full.flow = flow;
    
    if logical(ipR.verbose)
        disp('--------------------------------------------------')
        disp('NexIS:global minimizing quadratic error at all time stamps with fmincon')
        disp(' ')
        disp(['Optimal seed rescale value = ' num2str(param_num(1))])
        disp(['Optimal alpha = ' num2str(param_num(2))])
        disp(['Optimal beta = ' num2str(param_num(3))])
        if logical(ipR.w_dir)
            disp(['Optimal s = ' num2str(param_num(4))])
        end
        disp(' ')

        disp('R values at each time stamp')
        disp(Rvalues)
        disp(' ')
        if strcmp(ipR.costfun,'LinR')
            disp(['Cost Function = ' num2str(length(time_stamps)) ' - sum(LinR)'])
        else
            disp(ipR.costfun)
        end
        disp(fval_num)
        disp(['AIC = ' num2str(outputs.nexis_global.Full.results.lm_AIC)])
        disp(['BIC = ' num2str(outputs.nexis_global.Full.results.lm_BIC)])
        disp(['Intercept = ' num2str(outputs.nexis_global.Full.results.lm_intercept)])
        disp(['pValue = ' num2str(outputs.nexis_global.Full.results.lm_pval)])
        disp(['Rsqr_ord = ' num2str(outputs.nexis_global.Full.results.lm_Rsquared_ord)])
        disp(['Rsqr_adj = ' num2str(outputs.nexis_global.Full.results.lm_Rsquared_adj)])
        disp(' ')
    end
else % With bootstrapping of parameters
    rng(0);
    for i = 1:ipR.niters
        fldname = sprintf('Iter_%d',i);
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
                pathology = pathology(2:end);
                time_stamps = time_stamps(2:end) - time_stamps(1); % offset model times from t1
                ipR.param_init(1) = 1; ipR.ub(1) = 1; ipR.lb(1) = 1; % don't fit gamma
            end
        else % Use pathology and seed as is
            pathology_raw = data426.(ipR.study);
            seed_location = seed426.(ipR.study);
            pathology = normalizer(pathology_raw,ipR.normtype);
            pathology_orig = pathology;
            time_stamps_orig = time_stamps;
        end
        U = zeros(size(C,1),1);
        fprintf('NexIS:global Bootstrapping Iteration %d/%d\n',i,ipR.niters);
        notnaninds = find(~isnan(pathology(:,1)));
        settonansize = round((1-ipR.resample_rate)*length(notnaninds));
        settonaninds = randperm(length(notnaninds));
        settonaninds = notnaninds(settonaninds(1:settonansize));
        pathology(settonaninds,:) = NaN;
        pathology = normalizer(pathology,ipR.normtype);

        if isnan(ipR.param_init(1))
            ipR.param_init(1) = sum(pathology(:,1),'omitnan')/nnz(seed_location); % heuristic default, study-dependent
        end

        if ~logical(ipR.w_dir)
            ipR.param_init(4) = 0.5;
            ipR.ub(4) = 0.5;
            ipR.lb(4) = 0.5;
        end
        mordervec = zeros(1,length(ipR.param_init));
        for m = 1:length(ipR.param_init)
            inclparam = (ipR.lb(m) ~= ipR.ub(m));
            mordervec(m) = inclparam;
        end
        morder = 1 + sum(mordervec);
        param_init = [ipR.param_init,zeros(1,3)]; % dummy values for a, b, and p
        lb = [ipR.lb,zeros(1,3)]; % dummy values for a, b, and p
        ub = [ipR.ub,zeros(1,3)]; % dummy values for a, b, and p
        
        objfun_handle = @(param) objfun_eNDM_general_dir_costopts(param,...
            seed_location,pathology,time_stamps,C,U,ipR.solvetype,ipR.volcorrect,...
            ipR.costfun,ipR.excltpts_costfun,ipR.exclseed_costfun);
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
        
        % Solve eNDM with the optimal parameters
        x0_num = seed_location*param_num(1);
        alpha_num = param_num(2);
        beta_num = param_num(3);
        s_num = param_num(4); % not fit if directionality is turned off
        a_num = param_num(5); % not fit
        b_num = param_num(6); % not fit
        p_num = param_num(7); % not fit
        ynum = eNDM_general_dir(x0_num,time_stamps,C,U,alpha_num,beta_num,s_num,a_num,b_num,p_num,ipR.solvetype,ipR.volcorrect);
        
        outputs.nexis_global.(fldname).data = pathology;
        outputs.nexis_global.(fldname).time_stamps = time_stamps;
        outputs.nexis_global.(fldname).predicted = ynum;
        outputs.nexis_global.(fldname).param_fit = param_num;
        outputs.nexis_global.(fldname).fval = fval_num;
        outputs.nexis_global.(fldname).init.C = C;
        if ismember(ipR.study,{'human','mouse'})
            outputs.nexis_global.(fldname).init.study = ['asyn ' ipR.study];
        else
            outputs.nexis_global.(fldname).init.study = ipR.study;
        end
        outputs.nexis_global.(fldname).init.solvetype = ipR.solvetype;
        outputs.nexis_global.(fldname).init.volcorrect = ipR.volcorrect;
        outputs.nexis_global.(fldname).init.normtype = ipR.normtype;
        outputs.nexis_global.(fldname).init.costfun = ipR.costfun;
        outputs.nexis_global.(fldname).init.exclseed_costfun = ipR.exclseed_costfun;
        outputs.nexis_global.(fldname).init.excltpts_costfun = ipR.excltpts_costfun;
        outputs.nexis_global.(fldname).init.w_dir = ipR.w_dir;
        outputs.nexis_global.(fldname).init.param_init = ipR.param_init;
        outputs.nexis_global.(fldname).init.ub = ipR.ub;
        outputs.nexis_global.(fldname).init.lb = ipR.lb;
        outputs.nexis_global.(fldname).init.bootstrapping = ipR.bootstrapping;
        outputs.nexis_global.(fldname).init.resample_rate = ipR.resample_rate;
        outputs.nexis_global.(fldname).init.niters = ipR.niters;
        outputs.nexis_global.(fldname).fmincon.optimality_tolerance = ipR.opttol;
        outputs.nexis_global.(fldname).fmincon.function_tolerance = ipR.fxntol;
        outputs.nexis_global.(fldname).fmincon.step_tolerance = ipR.steptol;
        outputs.nexis_global.(fldname).fmincon.algorithm = ipR.algo;
        outputs.nexis_global.(fldname).fmincon.max_evaluations = ipR.maxeval;
        
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
        outputs.nexis_global.(fldname).results.Corrs = Rvalues;
        P = reshape(pathology, [], 1);
        Y = reshape(ynum, [], 1);
        numObs1 = length(P(~isnan(P)));
        lm_nexis = fitlm(Y, P);
        logL = lm_nexis.LogLikelihood;
        outputs.nexis_global.(fldname).results.lm_LogL = logL;
        outputs.nexis_global.(fldname).results.lm_AIC = -2*logL + 2*morder;
        outputs.nexis_global.(fldname).results.lm_BIC = -2*logL + log(numObs1)*morder;
        outputs.nexis_global.(fldname).results.lm_intercept = lm_nexis.Coefficients.Estimate(1);
        outputs.nexis_global.(fldname).results.lm_pval = lm_nexis.Coefficients.pValue(2);
        outputs.nexis_global.(fldname).results.lm_Rsquared_ord = lm_nexis.Rsquared.Ordinary;
        outputs.nexis_global.(fldname).results.lm_Rsquared_adj = lm_nexis.Rsquared.Adjusted;
        if logical(ipR.verbose)
            disp('--------------------------------------------------')
            disp('General eNDM minimizing quadratic error at all time stamps with fmincon')
            disp(' ')
            disp(['Fit seed rescale value = ' num2str(param_num(1))])
            disp(['Fit alpha = ' num2str(param_num(2))])
            disp(['Fit beta = ' num2str(param_num(3))])
            if logical(ipR.w_dir)
                disp(['Fit s = ' num2str(param_num(4))])
            end
            disp(' ')

            disp('R values at each time stamp')
            disp(Rvalues)
            disp(' ')
            if strcmp(ipR.costfun,'LinR')
                disp(['Cost Function = ' num2str(length(time_stamps)) ' - sum(LinR)'])
            else
                disp(ipR.costfun)
            end
            disp(fval_num)
            disp(['AIC = ' num2str(outputs.nexis_global.(fldname).results.lm_AIC)])
            disp(['BIC = ' num2str(outputs.nexis_global.(fldname).results.lm_BIC)])
            disp(['Intercept = ' num2str(outputs.nexis_global.(fldname).results.lm_intercept)])
            disp(['pValue = ' num2str(outputs.nexis_global.(fldname).results.lm_pval)])
            disp(['Rsqr_ord = ' num2str(outputs.nexis_global.(fldname).results.lm_Rsquared_ord)])
            disp(['Rsqr_adj = ' num2str(outputs.nexis_global.(fldname).results.lm_Rsquared_adj)])
            disp(' ')
        end
    end
    
    fprintf('Creating Optimal NDM Model\n');
    time_stamps = tpts.(ipR.study);
    pathology = normalizer(data426.(ipR.study),ipR.normtype);   
    % Yuanxi's comment: for testing the program
%     pathology = pathology/nansum(pathology(:,1));


    seed_location = seed426.(ipR.study);
    fldnames = fieldnames(outputs.nexis_global);
    param_fits = zeros(length(fldnames),length(outputs.nexis_global.(fldnames{1}).param_fit));
    for i = 1:length(fldnames)
        fldname = fldnames{i};
        param_fits(i,:) = outputs.nexis_global.(fldname).param_fit;
    end
    param_opt = mean(param_fits);
    x0_opt = seed_location*param_opt(1);
    alpha_opt = param_opt(2);
    beta_opt = param_opt(3);
    s_opt = param_opt(4); % not fit if directionality is turned off
    a_opt = param_opt(5); % not fit
    b_opt = param_opt(6); % not fit
    p_opt = param_opt(7); % not fit
    yopt = eNDM_general_dir(x0_opt,time_stamps,C,U,alpha_opt,beta_opt,s_opt,a_opt,b_opt,p_opt,ipR.solvetype,ipR.volcorrect);
    outputs.nexis_global.Full.data = pathology;
    outputs.nexis_global.Full.time_stamps = time_stamps;
    outputs.nexis_global.Full.predicted = yopt;
    outputs.nexis_global.Full.param_fit = param_opt;
    outputs.nexis_global.Full.fval = objfun_eNDM_general_dir_costopts(param_opt,...
        seed_location,pathology,time_stamps,C,U,ipR.solvetype,ipR.volcorrect,...
        ipR.costfun,ipR.excltpts_costfun,ipR.exclseed_costfun);
    outputs.nexis_global.Full.init = outputs.nexis_global.(fldnames{1}).init;
    outputs.nexis_global.Full.fmincon = [];
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
    outputs.nexis_global.Full.results.Corrs = Rvalues;
    P = reshape(pathology, [], 1);
    Y = reshape(yopt, [], 1);
    numObs1 = length(P(~isnan(P)));
    lm_nexis = fitlm(Y, P);
    logL = lm_nexis.LogLikelihood;
    outputs.nexis_global.Full.results.lm_LogL = logL;
    outputs.nexis_global.Full.results.lm_AIC = -2*logL + 2*morder;
    outputs.nexis_global.Full.results.lm_BIC = -2*logL + log(numObs1)*morder;
    outputs.nexis_global.Full.results.lm_intercept = lm_nexis.Coefficients.Estimate(1);
    outputs.nexis_global.Full.results.lm_pval = lm_nexis.Coefficients.pValue(1);
    outputs.nexis_global.Full.results.lm_Rsquared_ord = lm_nexis.Rsquared.Ordinary;
    outputs.nexis_global.Full.results.lm_Rsquared_adj = lm_nexis.Rsquared.Adjusted;
    
    % flow = FlowCalculator(yopt,C,beta_opt,1,U,b_opt);
    % for i = 1:size(flow,3)
    %     flow_ = flow(:,:,i);
    %     flow_(flow_ < prctile(nonzeros(flow),ipR.flowthresh)) = 0;
    %     flow(:,:,i) = flow_;
    % end
    % outputs.nexis_global.Full.flow = flow;    
    
    if logical(ipR.verbose)
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

        disp('R values at each time stamp')
        disp(Rvalues)
        disp(' ')

        disp(['AIC = ' num2str(outputs.nexis_global.Full.results.lm_AIC)])
        disp(['BIC = ' num2str(outputs.nexis_global.Full.results.lm_BIC)])
        disp(['Intercept = ' num2str(outputs.nexis_global.Full.results.lm_intercept)])
        disp(['pValue = ' num2str(outputs.nexis_global.Full.results.lm_pval)])
        disp(['Rsqr_ord = ' num2str(outputs.nexis_global.Full.results.lm_Rsquared_ord)])
        disp(['Rsqr_adj = ' num2str(outputs.nexis_global.Full.results.lm_Rsquared_adj)])
        disp(' ')
    end

end

    function normdata = normalizer(data,ntype)
        if strcmp(ntype,'sum')
            normdata = data/nansum(data(:,1));
        elseif strcmp(ntype,'mean')
            normdata = data/nanmean(data(:,1));
        elseif strcmp(ntype,'norm2')
            normdata = data/norm(data(:,1),2);
        elseif strcmp(ntype,'log')
            normdata = log(data+1);
        else
            normdata = data; 
        end
    end
end