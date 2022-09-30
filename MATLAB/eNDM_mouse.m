function outputs = eNDM_mouse(varargin)

% Parameters to fit
% param(1) = seed rescale factor
% param(2) = alpha
% param(3) = beta
% param(4) = s
% param(5) = a (not explored or fit)
% param(6) = b 
% param(7) = p

% Define defaults and set inputs
study_ = 'IbaHippInj'; 
costfun_ = 'LinR';
solvetype_ = 'analytic';
volcorrect_ = 0;
exclseed_costfun_ = 0;
excltpts_costfun_ = [];
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
outputs_ndm_ = [];
bounds_type_endm_ = 'old'; % 'old', 'CI_X' where X is the percent
bootstrapping_endm_ = 0;
resample_rate_endm_ = 0.8;
niters_endm_ = 100;
verbose_endm_ = 0;
fmindisplay_endm_ = 0;
datatype_endm_ = 'gene'; % 'gene', 'ct_tasic', 'ct_zeisel'
datalist_endm_ = 3578; % index for Trem2
datapca_endm_ = 0;
flowthresh_ = 99.93;

ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
% validNonnegative = @(x) isnumeric(x) && all(x(:) >= 0);
validBoolean = @(x) isscalar(x) && (x==0 || x==1);
validChar = @(x) ischar(x);
% validStudy = @(x) ismember(x,{'IbaHippInj','IbaStrInj','Clavaguera','Hurtado',...
%                                 'BolundaDSAD','BolundaCBD','DS4','DS6','DS9',...
%                                 'DS9_110','DS6_110','asyn_human','asyn_mouse'});
validST = @(x) ismember(x,{'analytic','numeric'});
validBoundsType = @(x) strcmp(x,'old') || strcmp(x(1:2),'CI');
validParam = @(x) (length(x) == 4);
validDataTypeENDM = @(x) ismember(x,{'gene','ct_tasic','ct_zeisel'});

addParameter(ip, 'study', study_, validChar);
addParameter(ip, 'costfun', costfun_, validChar);
addParameter(ip, 'solvetype', solvetype_, validST);
addParameter(ip, 'volcorrect', volcorrect_, validBoolean);
addParameter(ip, 'exclseed_costfun', exclseed_costfun_, validBoolean);
addParameter(ip, 'excltpts_costfun', excltpts_costfun_);
addParameter(ip, 'normtype', normtype_, validChar);
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
addParameter(ip, 'outputs_ndm', outputs_ndm_);
addParameter(ip, 'bounds_type_endm', bounds_type_endm_, validBoundsType);
addParameter(ip, 'bootstrapping_endm', bootstrapping_endm_, validBoolean);
addParameter(ip, 'resample_rate_endm', resample_rate_endm_, validScalar);
addParameter(ip, 'niters_endm', niters_endm_, validScalar);
addParameter(ip, 'verbose_endm', verbose_endm_, validBoolean);
addParameter(ip, 'fmindisplay_endm', fmindisplay_endm_, validBoolean);
addParameter(ip, 'datatype_endm', datatype_endm_, validDataTypeENDM);
addParameter(ip, 'datalist_endm', datalist_endm_);
addParameter(ip, 'datapca_endm', datapca_endm_, validBoolean);
addParameter(ip, 'flowthresh', flowthresh_, validScalar);

parse(ip, varargin{:});
ipR = ip.Results;

outputs = ipR.outputs_ndm;
if isempty(outputs)
   outputs = stdNDM_mouse('study', ipR.study,...
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

% Load in data from eNDM/raw_data_mouse directory
if strcmp(ipR.study(1:2),'DS')
    load([cd filesep 'raw_data_mouse' filesep 'KaufmanDiamond_datasets_dat&seed.mat'],...
        'data426','seed426','tpts');
elseif strcmp(ipR.study(1:4),'asyn')
    load([cd filesep 'raw_data_mouse' filesep 'mouse_aSynData_426.mat'],...
        'data426','seed426','tpts');
    ipR.study = ipR.study(6:end);
else
    load([cd filesep 'raw_data_mouse' filesep 'eNDM_mousedata.mat'],...
        'data426','seed426','tpts');
end
load([cd filesep 'raw_data_mouse' filesep 'eNDM_mousedata.mat'],'Networks');

% Define connectome
C = Networks.ret;

% Normalize C (minmax)
cmax = max(max(C));
cmin = min(min(C));
C = (C - cmin)./(cmax-cmin);

% Define cell type matrix, U
if ~isequal(ipR.datalist_endm,{'random'})
    if ~isnumeric(ipR.datalist_endm)
        ind_endm = NameIndex(ipR.datalist_endm,ipR.datatype_endm);
    else
        ind_endm = ipR.datalist_endm;
    end

    if strcmp(ipR.datatype_endm,'gene')
        load([cd filesep 'raw_data_mouse' filesep 'Regional_Gene_Data.mat'],'regvgene_mean');
        U = regvgene_mean(:,ind_endm);
    elseif strcmp(ipR.datatype_endm,'ct_tasic')
        load([cd filesep 'raw_data_mouse' filesep 'Tasic_CTMaps.mat'],'Tasic_ng606');
        U = Tasic_ng606(:,ind_endm);
    elseif strcmp(ipR.datatype_endm,'ct_zeisel')
        load([cd filesep 'raw_data_mouse' filesep 'Zeisel_CTMaps.mat'],'Zeisel_ng1360');
        U = Zeisel_ng1360(:,ind_endm);
    end
else
    U = rand(size(C,1),1);
end

U = U ./ nanmean(U);
if logical(ipR.datapca_endm) && (length(ipR.datalist_endm) > 1)
    U_mean = nanmean(U,2);
    [~, score, ~, ~, ~] = pca(U);
    U = score(:,1);
    if corr(U,U_mean) < 0
        U = -U;
    end
end
minU = repmat(min(U),size(U,1),1);
maxU = repmat(max(U),size(U,1),1);
U = (U - minU) ./ (maxU - minU);

% Solve and store results
outputs.endm = struct;
if ~logical(ipR.bootstrapping_endm)
    fprintf('Creating Optimal eNDM Model\n');
    time_stamps = tpts.(ipR.study);
    pathology = normalizer(data426.(ipR.study),ipR.normtype);
    seed_location = seed426.(ipR.study);
    n_types = size(U,2);
    ndmflds = fieldnames(outputs.ndm);
    if length(ndmflds) == 1
        param_inits = outputs.ndm.Full.param_fit; 
    else
        param_inits = zeros((length(ndmflds)-1),length(outputs.ndm.Full.param_fit));
        for k = 1:(length(ndmflds)-1)
            fld = ndmflds{k};
            param_inits(k,:) = outputs.ndm.(fld).param_fit;
        end    
    end

    if strcmp(ipR.bounds_type_endm,'old') && (size(param_inits,1) > 1)
        param_init = mean(param_inits);
        ub = 1.3*param_init;
        lb = 0.7*param_init;
    elseif (size(param_inits,1) == 1)
        param_init = param_inits;
        ub = 1.3*param_init;
        lb = 0.7*param_init;
    else
        prct = str2double(ipR.bounds_type_endm(4:end));
        param_init = median(param_inits);
        ub = prctile(param_inits,((100-prct)/2)+prct,1);
        lb = prctile(param_inits,((100-prct)/2),1);
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
    
    param_init = [param_init(1:4),zeros(1,n_types),zeros(1,n_types),zeros(1,n_types)];
    lb = [lb(1:4),zeros(1,n_types),-Inf(1,n_types),-Inf(1,n_types)];
    ub = [ub(1:4),zeros(1,n_types),Inf(1,n_types),Inf(1,n_types)];
    
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
    [param_num, fval_num] = fmincon(objfun_handle,param_init,[],[],[],[],lb,ub,[],options);

    % Solve eNDM with the optimal parameters
    x0_num = seed_location*param_num(1);
    alpha_num = param_num(2);
    beta_num = param_num(3);
    s_num = param_num(4); % not fit if directionality is turned off
    a_num = param_num(5:(n_types+4)); % not fit
    b_num = param_num((n_types+5):(2*n_types+4));
    p_num = param_num((2*n_types+5):(3*n_types+4));
    ynum = eNDM_general_dir(x0_num,time_stamps,C,U,alpha_num,beta_num,s_num,a_num,b_num,p_num,ipR.solvetype,ipR.volcorrect);
    outputs.endm.Full.data = pathology;
    outputs.endm.Full.time_stamps = time_stamps;
    outputs.endm.Full.predicted = ynum;
    outputs.endm.Full.param_fit = param_num;
    outputs.endm.Full.fval = fval_num;
    outputs.endm.Full.init.C = C;
    if ismember(ipR.study,{'human','mouse'})
        outputs.endm.Full.init.study = ['asyn ' ipR.study];
    else
        outputs.endm.Full.init.study = ipR.study;
    end
    outputs.endm.Full.init.solvetype = ipR.solvetype;
    outputs.endm.Full.init.volcorrect = ipR.volcorrect;
    outputs.endm.Full.init.normtype = ipR.normtype;
    outputs.endm.Full.init.costfun = ipR.costfun;
    outputs.endm.Full.init.exclseed_costfun = ipR.exclseed_costfun;
    outputs.endm.Full.init.excltpts_costfun = ipR.excltpts_costfun;
    outputs.endm.Full.init.w_dir = ipR.w_dir;
    outputs.endm.Full.init.param_init = param_init;
    outputs.endm.Full.init.ub = ub;
    outputs.endm.Full.init.lb = lb;
    outputs.endm.Full.init.bootstrapping_endm = ipR.bootstrapping_endm;
    outputs.endm.Full.init.resample_rate_endm = ipR.resample_rate_endm;
    outputs.endm.Full.init.bounds_type_endm = ipR.bounds_type_endm;
    outputs.endm.Full.init.datatype_endm = ipR.datatype_endm;
    outputs.endm.Full.init.datalist_endm = ipR.datalist_endm;
    outputs.endm.Full.init.datapca_endm = ipR.datapca_endm;
    outputs.endm.Full.init.niters_endm = ipR.niters_endm;
    outputs.endm.Full.fmincon.optimality_tolerance = ipR.opttol;
    outputs.endm.Full.fmincon.function_tolerance = ipR.fxntol;
    outputs.endm.Full.fmincon.step_tolerance = ipR.steptol;
    outputs.endm.Full.fmincon.algorithm = ipR.algo;
    outputs.endm.Full.fmincon.max_evaluations = ipR.maxeval;

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
    outputs.endm.Full.results.Corrs = Rvalues;

    P = reshape(pathology, [], 1);
    Y = reshape(ynum, [], 1);
    numObs1 = length(P(~isnan(P)));
    lm_endm = fitlm(Y, P);
    logL = lm_endm.LogLikelihood;
    outputs.endm.Full.results.lm_LogL = logL;
    outputs.endm.Full.results.lm_AIC = -2*logL + 2*morder;
    outputs.endm.Full.results.lm_BIC = -2*logL + log(numObs1)*morder;
    outputs.endm.Full.results.lm_intercept = lm_endm.Coefficients.Estimate(1);
    outputs.endm.Full.results.lm_pval = lm_endm.Coefficients.pValue(1);
    outputs.endm.Full.results.lm_Rsquared_ord = lm_endm.Rsquared.Ordinary;
    outputs.endm.Full.results.lm_Rsquared_adj = lm_endm.Rsquared.Adjusted;
    flow = FlowCalculator(ynum,C,beta_num,0,U,b_num);
    for i = 1:size(flow,3)
        flow_ = flow(:,:,i);
        flow_(flow_ < prctile(nonzeros(flow),ipR.flowthresh)) = 0;
        flow(:,:,i) = flow_;
    end
    outputs.endm.Full.flow = flow;
    if logical(ipR.verbose_endm)
        disp('--------------------------------------------------')
        disp('General eNDM minimizing quadratic error at all time stamps with fmincon')
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
        disp(['AIC = ' num2str(outputs.endm.Full.results.lm_AIC)])
        disp(['BIC = ' num2str(outputs.endm.Full.results.lm_BIC)])
        disp(['Intercept = ' num2str(outputs.endm.Full.results.lm_intercept)])
        disp(['pValue = ' num2str(outputs.endm.Full.results.lm_pval)])
        disp(['Rsqr_ord = ' num2str(outputs.endm.Full.results.lm_Rsquared_ord)])
        disp(['Rsqr_adj = ' num2str(outputs.endm.Full.results.lm_Rsquared_adj)])
        disp(' ')
    end
else
    rng(1);
    for i = 1:ipR.niters_endm
        fldname = sprintf('Iter_%d',i);
        fprintf('eNDM Bootstrapping Iteration %d/%d\n',i,ipR.niters_endm);        
        time_stamps = tpts.(ipR.study);
        pathology = data426.(ipR.study);
        seed_location = seed426.(ipR.study);
        n_types = size(U,2);

        notnaninds = find(~isnan(pathology(:,1)));
        settonansize = round((1-ipR.resample_rate)*length(notnaninds));
        settonaninds = randperm(length(notnaninds));
        settonaninds = notnaninds(settonaninds(1:settonansize));
        pathology(settonaninds,:) = NaN;
        pathology = normalizer(pathology,ipR.normtype);
        ndmflds = fieldnames(outputs.ndm);
        if length(ndmflds) == 1
            param_inits = outputs.ndm.Full.param_fit; 
        else
            param_inits = zeros((length(ndmflds)-1),length(outputs.ndm.Full.param_fit));
            for k = 1:(length(ndmflds)-1)
                fld = ndmflds{k};
                param_inits(k,:) = outputs.ndm.(fld).param_fit;
            end    
        end

        if strcmp(ipR.bounds_type_endm,'old') && (size(param_inits,1) > 1)
            param_init = mean(param_inits);
            ub = 1.3*param_init;
            lb = 0.7*param_init;
        elseif (size(param_inits,1) == 1)
            param_init = param_inits;
            ub = 1.3*param_init;
            lb = 0.7*param_init;
        else
            prct = str2double(ipR.bounds_type_endm(4:end)); 
            param_init = median(param_inits);
            ub = prctile(param_inits,((100-prct)/2)+prct,1);
            lb = prctile(param_inits,((100-prct)/2),1);
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
        
        param_init = [param_init(1:4),zeros(1,n_types),zeros(1,n_types),zeros(1,n_types)];
        lb = [lb(1:4),zeros(1,n_types),-Inf(1,n_types),-Inf(1,n_types)];
        ub = [ub(1:4),zeros(1,n_types),Inf(1,n_types),Inf(1,n_types)];

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
%             [param_num, fval_num] = fmincon(objfun_handle,rand(1,length(param_init)),...
%                 [],[],[],[],-Inf(1,length(param_init)),Inf(1,length(param_init)),[],options);
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
        a_num = param_num(5:(n_types+4)); % not fit
        b_num = param_num((n_types+5):(2*n_types+4));
        p_num = param_num((2*n_types+5):(3*n_types+4));
        ynum = eNDM_general_dir(x0_num,time_stamps,C,U,alpha_num,beta_num,s_num,a_num,b_num,p_num,ipR.solvetype,ipR.volcorrect);
        outputs.endm.(fldname).data = pathology;
        outputs.endm.(fldname).time_stamps = time_stamps;
        outputs.endm.(fldname).predicted = ynum;
        outputs.endm.(fldname).param_fit = param_num;
        outputs.endm.(fldname).fval = fval_num;
        outputs.endm.(fldname).init.C = C;
        outputs.endm.(fldname).init.study = ipR.study;
        if ismember(ipR.study,{'human','mouse'})
            outputs.endm.(fldname).init.study = ['asyn ' ipR.study];
        else
            outputs.endm.(fldname).init.study = ipR.study;
        end
        outputs.endm.(fldname).init.solvetype = ipR.solvetype;
        outputs.endm.(fldname).init.volcorrect = ipR.volcorrect;
        outputs.endm.(fldname).init.normtype = ipR.normtype;
        outputs.endm.(fldname).init.costfun = ipR.costfun;
        outputs.endm.(fldname).init.exclseed_costfun = ipR.exclseed_costfun;
        outputs.endm.(fldname).init.excltpts_costfun = ipR.excltpts_costfun;
        outputs.endm.(fldname).init.w_dir = ipR.w_dir;
        outputs.endm.(fldname).init.param_init = param_init;
        outputs.endm.(fldname).init.ub = ub;
        outputs.endm.(fldname).init.lb = lb;
        outputs.endm.(fldname).init.bootstrapping_endm = ipR.bootstrapping_endm;
        outputs.endm.(fldname).init.resample_rate_endm = ipR.resample_rate_endm;
        outputs.endm.(fldname).init.bounds_type_endm = ipR.bounds_type_endm;
        outputs.endm.(fldname).init.datatype_endm = ipR.datatype_endm;
        outputs.endm.(fldname).init.datalist_endm = ipR.datalist_endm;
        outputs.endm.(fldname).init.datapca_endm = ipR.datapca_endm;
        outputs.endm.(fldname).init.niters_endm = ipR.niters_endm;
        outputs.endm.(fldname).fmincon.optimality_tolerance = ipR.opttol;
        outputs.endm.(fldname).fmincon.function_tolerance = ipR.fxntol;
        outputs.endm.(fldname).fmincon.step_tolerance = ipR.steptol;
        outputs.endm.(fldname).fmincon.algorithm = ipR.algo;
        outputs.endm.(fldname).fmincon.max_evaluations = ipR.maxeval;

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
        outputs.endm.(fldname).results.Corrs = Rvalues;

        P = reshape(pathology, [], 1);
        Y = reshape(ynum, [], 1);
        numObs1 = length(P(~isnan(P)));
        lm_endm = fitlm(Y, P);
        logL = lm_endm.LogLikelihood;
        outputs.endm.(fldname).results.lm_LogL = logL;
        outputs.endm.(fldname).results.lm_AIC = -2*logL + 2*morder;
        outputs.endm.(fldname).results.lm_BIC = -2*logL + log(numObs1)*morder;
        outputs.endm.(fldname).results.lm_intercept = lm_endm.Coefficients.Estimate(1);
        outputs.endm.(fldname).results.lm_pval = lm_endm.Coefficients.pValue(1);
        outputs.endm.(fldname).results.lm_Rsquared_ord = lm_endm.Rsquared.Ordinary;
        outputs.endm.(fldname).results.lm_Rsquared_adj = lm_endm.Rsquared.Adjusted;
        if logical(ipR.verbose_endm)
            disp('--------------------------------------------------')
            disp('General eNDM minimizing quadratic error at all time stamps with fmincon')
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
            disp(['AIC = ' num2str(outputs.endm.(fldname).results.lm_AIC)])
            disp(['BIC = ' num2str(outputs.endm.(fldname).results.lm_BIC)])
            disp(['Intercept = ' num2str(outputs.endm.(fldname).results.lm_intercept)])
            disp(['pValue = ' num2str(outputs.endm.(fldname).results.lm_pval)])
            disp(['Rsqr_ord = ' num2str(outputs.endm.(fldname).results.lm_Rsquared_ord)])
            disp(['Rsqr_adj = ' num2str(outputs.endm.(fldname).results.lm_Rsquared_adj)])
            disp(' ')
        end
        clear ub lb param_init
    end
    fprintf('Creating Optimal eNDM Model\n');
    time_stamps = tpts.(ipR.study);
    pathology = normalizer(data426.(ipR.study),ipR.normtype);  
    seed_location = seed426.(ipR.study);
    fldnames = fieldnames(outputs.endm);
    param_fits = zeros(length(fldnames),length(outputs.endm.(fldnames{1}).param_fit));
    for i = 1:length(fldnames)
        fldname = fldnames{i};
        param_fits(i,:) = outputs.endm.(fldname).param_fit;
    end
    param_opt = mean(param_fits);
    x0_opt = seed_location*param_opt(1);
    alpha_opt = param_opt(2);
    beta_opt = param_opt(3);
    s_opt = param_opt(4); % not fit if directionality is turned off
    a_opt = param_opt(5:(n_types+4)); % not fit
    b_opt = param_opt((n_types+5):(2*n_types+4)); 
    p_opt = param_opt((2*n_types+5):(3*n_types+4)); 
    yopt = eNDM_general_dir(x0_opt,time_stamps,C,U,alpha_opt,beta_opt,s_opt,a_opt,b_opt,p_opt,ipR.solvetype,ipR.volcorrect);
    outputs.endm.Full.data = pathology;
    outputs.endm.Full.time_stamps = time_stamps;
    outputs.endm.Full.predicted = yopt;
    outputs.endm.Full.param_fit = param_opt;
    outputs.endm.Full.fval = objfun_eNDM_general_dir_costopts(param_opt,...
        seed_location,pathology,time_stamps,C,U,ipR.solvetype,ipR.volcorrect,...
        ipR.costfun,ipR.excltpts_costfun,ipR.exclseed_costfun);
    outputs.endm.Full.init = outputs.endm.(fldnames{1}).init;
    outputs.endm.Full.fmincon = [];
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
    outputs.endm.Full.results.Corrs = Rvalues;
    P = reshape(pathology, [], 1);
    Y = reshape(yopt, [], 1);
    numObs1 = length(P(~isnan(P)));
    lm_endm = fitlm(Y, P);
    logL = lm_endm.LogLikelihood;
    outputs.endm.Full.results.lm_LogL = logL;
    outputs.endm.Full.results.lm_AIC = -2*logL + 2*morder;
    outputs.endm.Full.results.lm_BIC = -2*logL + log(numObs1)*morder;
    outputs.endm.Full.results.lm_intercept = lm_endm.Coefficients.Estimate(1);
    outputs.endm.Full.results.lm_pval = lm_endm.Coefficients.pValue(1);
    outputs.endm.Full.results.lm_Rsquared_ord = lm_endm.Rsquared.Ordinary;
    outputs.endm.Full.results.lm_Rsquared_adj = lm_endm.Rsquared.Adjusted;
    
    flow = FlowCalculator(yopt,C,beta_opt,0,U,b_opt);
    for i = 1:size(flow,3)
        flow_ = flow(:,:,i);
        flow_(flow_ < prctile(nonzeros(flow),ipR.flowthresh)) = 0;
        flow(:,:,i) = flow_;
    end
    outputs.endm.Full.flow = flow;
    
    if logical(ipR.verbose_endm)
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

        disp(['AIC = ' num2str(outputs.endm.Full.results.lm_AIC)])
        disp(['BIC = ' num2str(outputs.endm.Full.results.lm_BIC)])
        disp(['Intercept = ' num2str(outputs.endm.Full.results.lm_intercept)])
        disp(['pValue = ' num2str(outputs.endm.Full.results.lm_pval)])
        disp(['Rsqr_ord = ' num2str(outputs.endm.Full.results.lm_Rsquared_ord)])
        disp(['Rsqr_adj = ' num2str(outputs.endm.Full.results.lm_Rsquared_adj)])
        disp(' ')
    end
end

    function indices = NameIndex(names,dattypeendm)
        if strcmp(dattypeendm,'gene')
            load([cd filesep 'raw_data_mouse' filesep 'gene_names_trans.mat'],'gene_names_trans');
            namescell = gene_names_trans;
        elseif strcmp(dattypeendm,'ct_tasic')
            load([cd filesep 'raw_data_mouse' filesep 'classkey_tasic.mat'],'classkey_tasic');
            namescell = classkey_tasic;
        elseif strcmp(dattypeendm,'ct_zeisel')
            load([cd filesep 'raw_data_mouse' filesep 'classkey_zeisel.mat'],'classkey_zeisel');
            namescell = classkey_zeisel;
        end
        
        indices = zeros(1,length(names));
        for n_i = 1:length(names)
            indices(n_i) = find(ismember(namescell,names(n_i)));
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