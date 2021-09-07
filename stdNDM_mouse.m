function outputs = stdNDM_mouse(varargin)
% Parameters to fit
% param(1) = seed rescale factor
% param(2) = alpha
% param(3) = beta
% param(4) = s

% Define defaults and set inputs
study_ = 'IbaHippInj'; 
costfun_ = 'LinR';
solvetype_ = 'analytic';
volcorrect_ = 0;
normtype_ = 'sum';
w_dir_ = 0; 
param_init_ = [NaN,0.5,1,0.5];
ub_ = [Inf,Inf,Inf,1];
lb_ = zeros(1,4);
bootstrapping_ = 0;
resample_rate_ = 0.8;
niters_ = 100;
verbose_ = 1;
fmindisplay_ = 1;
flowthresh_ = 99.93;

ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
% validNonnegative = @(x) isnumeric(x) && all(x(:) >= 0);
validBoolean = @(x) isscalar(x) && (x==0 || x==1);
validChar = @(x) ischar(x);
validStudy = @(x) ismember(x,{'IbaHippInj','IbaStrInj','Clavaguera','Hurtado',...
                                'BolundaDSAD','BolundaCBD','DS4','DS6','DS9',...
                                'asyn_human','asyn_mouse'});
validST = @(x) ismember(x,{'analytic','numeric'});
validParam = @(x) (length(x) == 4);

addParameter(ip, 'study', study_, validStudy);
addParameter(ip, 'costfun', costfun_, validChar);
addParameter(ip, 'solvetype', solvetype_, validST);
addParameter(ip, 'volcorrect', volcorrect_, validBoolean);
addParameter(ip, 'normtype', normtype_, validChar);
addParameter(ip, 'w_dir', w_dir_, validBoolean);
addParameter(ip, 'param_init', param_init_, validParam);
addParameter(ip, 'ub', ub_, validParam);
addParameter(ip, 'lb', lb_, validParam);
addParameter(ip, 'bootstrapping', bootstrapping_, validBoolean);
addParameter(ip, 'resample_rate', resample_rate_, validScalar);
addParameter(ip, 'niters', niters_, validScalar);
addParameter(ip, 'verbose', verbose_, validBoolean);
addParameter(ip, 'fmindisplay', fmindisplay_, validBoolean);
addParameter(ip, 'flowthresh', flowthresh_, validScalar);
parse(ip, varargin{:});
ipR = ip.Results;

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

% Defining LinR (Lin 1989)
% LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

% Define connectome
C = Networks.ret;

% Normalize C (minmax)
cmax = max(max(C));
cmin = min(min(C));
C = (C - cmin)./(cmax-cmin);

% Solve and store results
outputs.ndm = struct;
if ~logical(ipR.bootstrapping)
    fprintf('Creating Optimal NDM Model\n');
    time_stamps = tpts.(ipR.study);
    pathology = normalizer(data426.(ipR.study),ipR.normtype);
    seed_location = seed426.(ipR.study);
    U = zeros(size(C,1),1);
    morder = 1+length(ipR.param_init);
    if isnan(ipR.param_init(1))
        ipR.param_init(1) = nansum(pathology(:,1))/nnz(seed_location); % heuristic default, study-dependent
    end

    if ~logical(ipR.w_dir)
        ipR.param_init(4) = 0.5;
        ipR.ub(4) = 0.5;
        ipR.lb(4) = 0.5;
        morder = morder - 1;
    end
    param_init = [ipR.param_init,zeros(1,3)]; % dummy values for a, b, and p
    lb = [ipR.lb,zeros(1,3)]; % dummy values for a, b, and p
    ub = [ipR.ub,zeros(1,3)]; % dummy values for a, b, and p
    
    opttol = 1e-8; fxntol = 1e-8; steptol = 1e-12; algo = 'sqp'; maxeval = 10000;
    objfun_handle = @(param) objfun_eNDM_general_dir_costopts(param,...
        seed_location,pathology,time_stamps,C,U,ipR.solvetype,ipR.volcorrect,ipR.costfun);
    if logical(ipR.fmindisplay)
        options = optimoptions(@fmincon,'Display','iter','Algorithm',algo,'MaxFunctionEvaluations',maxeval,...
                'OptimalityTolerance',opttol,'FunctionTolerance',fxntol,'StepTolerance',steptol);
    else
        options = optimoptions(@fmincon,'Algorithm',algo,'MaxFunctionEvaluations',maxeval,...
                'OptimalityTolerance',opttol,'FunctionTolerance',fxntol,'StepTolerance',steptol);
    end
    [param_num, fval_num] = fmincon(objfun_handle,param_init,[],[],[],[],lb,ub,[],options);

    % Solve eNDM with the optimal parameters
    x0_num = seed_location*param_num(1);
    alpha_num = param_num(2);
    beta_num = param_num(3);
    s_num = param_num(4); % not fit if directionality is turned off
    a_num = param_num(5); % not fit
    b_num = param_num(6); % not fit
    p_num = param_num(7); % not fit
    ynum = eNDM_general_dir(x0_num,time_stamps,C,U,alpha_num,beta_num,s_num,a_num,b_num,p_num,ipR.solvetype,ipR.volcorrect);
    outputs.ndm.Full.data = pathology;
    outputs.ndm.Full.time_stamps = time_stamps;
    outputs.ndm.Full.predicted = ynum;
    outputs.ndm.Full.param_fit = param_num;
    outputs.ndm.Full.fval = fval_num;
    outputs.ndm.Full.init.C = C;
    if ismember(ipR.study,{'human','mouse'})
        outputs.ndm.Full.init.study = ['asyn ' ipR.study];
    else
        outputs.ndm.Full.init.study = ipR.study;
    end
    outputs.ndm.Full.init.solvetype = ipR.solvetype;
    outputs.ndm.Full.init.volcorrect = ipR.volcorrect;
    outputs.ndm.Full.init.normtype = ipR.normtype;
    outputs.ndm.Full.init.costfun = ipR.costfun;
    outputs.ndm.Full.init.w_dir = ipR.w_dir;
    outputs.ndm.Full.init.param_init = ipR.param_init;
    outputs.ndm.Full.init.ub = ipR.ub;
    outputs.ndm.Full.init.lb = ipR.lb;
    outputs.ndm.Full.init.bootstrapping = ipR.bootstrapping;
    outputs.ndm.Full.init.resample_rate = ipR.resample_rate;
    outputs.ndm.Full.init.niters = ipR.niters;
    outputs.ndm.Full.fmincon.optimality_tolerance = opttol;
    outputs.ndm.Full.fmincon.function_tolerance = fxntol;
    outputs.ndm.Full.fmincon.step_tolerance = steptol;
    outputs.ndm.Full.fmincon.algorithm = algo;
    outputs.ndm.Full.fmincon.max_evaluations = maxeval;

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
    outputs.ndm.Full.results.Corrs = Rvalues;

    P = reshape(pathology, [], 1);
    Y = reshape(ynum, [], 1);
    numObs1 = length(P(~isnan(P)));
    lm_endm = fitlm(Y, P);
    logL = lm_endm.LogLikelihood;
    outputs.ndm.Full.results.lm_AIC = -2*logL + 2*morder;
    outputs.ndm.Full.results.lm_BIC = -2*logL + log(numObs1)*morder;
    outputs.ndm.Full.results.lm_intercept = lm_endm.Coefficients.Estimate(1);
    outputs.ndm.Full.results.lm_pval = lm_endm.Coefficients.pValue(1);
    outputs.ndm.Full.results.lm_Rsquared_ord = lm_endm.Rsquared.Ordinary;
    outputs.ndm.Full.results.lm_Rsquared_adj = lm_endm.Rsquared.Adjusted;
    
    flow = FlowCalculator(ynum,C,beta_num,1,U,b_num);
    for i = 1:size(flow,3)
        flow_ = flow(:,:,i);
        flow_(flow_ < prctile(nonzeros(flow),ipR.flowthresh)) = 0;
        flow(:,:,i) = flow_;
    end
    outputs.ndm.Full.flow = flow;
    
    if logical(ipR.verbose)
        disp('--------------------------------------------------')
        disp('General eNDM minimizing quadratic error at all time stamps with fmincon')
        disp(' ')
        disp(['Optimal seed rescale value = ' num2str(param_num(1))])
        disp(['Optimal alpha = ' num2str(param_num(2))])
        disp(['Optimal beta = ' num2str(param_num(3))])
        if logical(ipR.w_dir)
            disp(['Optimal s = ' num2str(param_num(4))])
        end
        disp(' ')

        % for i = 1:length(type_inds)
            % disp(['Cell Type - ' classkey{type_inds(i)}])
            % disp(['Optimal a = ' num2str(param_num(3 + i))])
            % disp(['Optimal b = ' num2str(param_num(3 + n_types + i))])
            % disp(['Optimal p = ' num2str(param_num(3 + 2*n_types + i))])
            % disp(' ')
        % end

        disp('R values at each time stamp')
        disp(Rvalues)
        disp(' ')
        if strcmp(ipR.costfun,'LinR')
            disp(['Cost Function = ' num2str(length(time_stamps)) ' - sum(LinR)'])
        else
            disp(ipR.costfun)
        end
        disp(fval_num)
        disp(['AIC = ' num2str(outputs.ndm.Full.results.lm_AIC)])
        disp(['BIC = ' num2str(outputs.ndm.Full.results.lm_BIC)])
        disp(['Intercept = ' num2str(outputs.ndm.Full.results.lm_intercept)])
        disp(['pValue = ' num2str(outputs.ndm.Full.results.lm_pval)])
        disp(['Rsqr_ord = ' num2str(outputs.ndm.Full.results.lm_Rsquared_ord)])
        disp(['Rsqr_adj = ' num2str(outputs.ndm.Full.results.lm_Rsquared_adj)])
        disp(' ')
    end
else
    rng(0);
    for i = 1:ipR.niters
        fldname = sprintf('Iter_%d',i);
        time_stamps = tpts.(ipR.study);
        pathology = data426.(ipR.study);    
        seed_location = seed426.(ipR.study);
        U = zeros(size(C,1),1);
        morder = 1+length(ipR.param_init);
        fprintf('NDM Bootstrapping Iteration %d/%d\n',i,ipR.niters);
        notnaninds = find(~isnan(pathology(:,1)));
        settonansize = round((1-ipR.resample_rate)*length(notnaninds));
        settonaninds = randperm(length(notnaninds));
        settonaninds = notnaninds(settonaninds(1:settonansize));
        pathology(settonaninds,:) = NaN;
        pathology = normalizer(pathology,ipR.normtype);
        if isnan(ipR.param_init(1))
            ipR.param_init(1) = nansum(pathology(:,1))/nnz(seed_location); % heuristic default, study-dependent
        end

        if ~logical(ipR.w_dir)
            ipR.param_init(4) = 0.5;
            ipR.ub(4) = 0.5;
            ipR.lb(4) = 0.5;
            morder = morder - 1;
        end
        param_init = [ipR.param_init,zeros(1,3)]; % dummy values for a, b, and p
        lb = [ipR.lb,zeros(1,3)]; % dummy values for a, b, and p
        ub = [ipR.ub,zeros(1,3)]; % dummy values for a, b, and p
        
        opttol = 1e-8; fxntol = 1e-8; steptol = 1e-12; algo = 'sqp'; maxeval = 10000;
        objfun_handle = @(param) objfun_eNDM_general_dir_costopts(param,...
            seed_location,pathology,time_stamps,C,U,ipR.solvetype,ipR.volcorrect,ipR.costfun);
        if logical(ipR.fmindisplay)
            options = optimoptions(@fmincon,'Display','final-detailed','Algorithm',algo,'MaxFunctionEvaluations',maxeval,...
                    'OptimalityTolerance',opttol,'FunctionTolerance',fxntol,'StepTolerance',steptol);
        else
            options = optimoptions(@fmincon,'MaxFunctionEvaluations',maxeval,...
                    'OptimalityTolerance',opttol);
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
        
        outputs.ndm.(fldname).data = pathology;
        outputs.ndm.(fldname).time_stamps = time_stamps;
        outputs.ndm.(fldname).predicted = ynum;
        outputs.ndm.(fldname).param_fit = param_num;
        outputs.ndm.(fldname).fval = fval_num;
        outputs.ndm.(fldname).init.C = C;
        if ismember(ipR.study,{'human','mouse'})
            outputs.ndm.(fldname).init.study = ['asyn ' ipR.study];
        else
            outputs.ndm.(fldname).init.study = ipR.study;
        end
        outputs.ndm.(fldname).init.solvetype = ipR.solvetype;
        outputs.ndm.(fldname).init.volcorrect = ipR.volcorrect;
        outputs.ndm.(fldname).init.normtype = ipR.normtype;
        outputs.ndm.(fldname).init.costfun = ipR.costfun;
        outputs.ndm.(fldname).init.w_dir = ipR.w_dir;
        outputs.ndm.(fldname).init.param_init = ipR.param_init;
        outputs.ndm.(fldname).init.ub = ipR.ub;
        outputs.ndm.(fldname).init.lb = ipR.lb;
        outputs.ndm.(fldname).init.bootstrapping = ipR.bootstrapping;
        outputs.ndm.(fldname).init.resample_rate = ipR.resample_rate;
        outputs.ndm.(fldname).init.niters = ipR.niters;
        outputs.ndm.(fldname).fmincon.optimality_tolerance = opttol;
        outputs.ndm.(fldname).fmincon.function_tolerance = fxntol;
        outputs.ndm.(fldname).fmincon.step_tolerance = steptol;
        outputs.ndm.(fldname).fmincon.algorithm = algo;
        outputs.ndm.(fldname).fmincon.max_evaluations = maxeval;
        
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
        outputs.ndm.(fldname).results.Corrs = Rvalues;
        P = reshape(pathology, [], 1);
        Y = reshape(ynum, [], 1);
        numObs1 = length(P(~isnan(P)));
        lm_endm = fitlm(Y, P);
        logL = lm_endm.LogLikelihood;
        outputs.ndm.(fldname).results.lm_AIC = -2*logL + 2*morder;
        outputs.ndm.(fldname).results.lm_BIC = -2*logL + log(numObs1)*morder;
        outputs.ndm.(fldname).results.lm_intercept = lm_endm.Coefficients.Estimate(1);
        outputs.ndm.(fldname).results.lm_pval = lm_endm.Coefficients.pValue(1);
        outputs.ndm.(fldname).results.lm_Rsquared_ord = lm_endm.Rsquared.Ordinary;
        outputs.ndm.(fldname).results.lm_Rsquared_adj = lm_endm.Rsquared.Adjusted;
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
            disp(['AIC = ' num2str(outputs.ndm.(fldname).results.lm_AIC)])
            disp(['BIC = ' num2str(outputs.ndm.(fldname).results.lm_BIC)])
            disp(['Intercept = ' num2str(outputs.ndm.(fldname).results.lm_intercept)])
            disp(['pValue = ' num2str(outputs.ndm.(fldname).results.lm_pval)])
            disp(['Rsqr_ord = ' num2str(outputs.ndm.(fldname).results.lm_Rsquared_ord)])
            disp(['Rsqr_adj = ' num2str(outputs.ndm.(fldname).results.lm_Rsquared_adj)])
            disp(' ')
        end
    end
    
    fprintf('Creating Optimal NDM Model\n');
    time_stamps = tpts.(ipR.study);
    pathology = normalizer(data426.(ipR.study),ipR.normtype);   
    seed_location = seed426.(ipR.study);
    fldnames = fieldnames(outputs.ndm);
    param_fits = zeros(length(fldnames),length(outputs.ndm.(fldnames{1}).param_fit));
    for i = 1:length(fldnames)
        fldname = fldnames{i};
        param_fits(i,:) = outputs.ndm.(fldname).param_fit;
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
    outputs.ndm.Full.data = pathology;
    outputs.ndm.Full.time_stamps = time_stamps;
    outputs.ndm.Full.predicted = yopt;
    outputs.ndm.Full.param_fit = param_opt;
    outputs.ndm.Full.fval = objfun_eNDM_general_dir_costopts(param_opt,...
        seed_location,pathology,time_stamps,C,U,ipR.solvetype,ipR.volcorrect,ipR.costfun);
    outputs.ndm.Full.init = outputs.ndm.(fldnames{1}).init;
    outputs.ndm.Full.fmincon = [];
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
    outputs.ndm.Full.results.Corrs = Rvalues;
    P = reshape(pathology, [], 1);
    Y = reshape(yopt, [], 1);
    numObs1 = length(P(~isnan(P)));
    lm_endm = fitlm(Y, P);
    logL = lm_endm.LogLikelihood;
    outputs.ndm.Full.results.lm_AIC = -2*logL + 2*morder;
    outputs.ndm.Full.results.lm_BIC = -2*logL + log(numObs1)*morder;
    outputs.ndm.Full.results.lm_intercept = lm_endm.Coefficients.Estimate(1);
    outputs.ndm.Full.results.lm_pval = lm_endm.Coefficients.pValue(1);
    outputs.ndm.Full.results.lm_Rsquared_ord = lm_endm.Rsquared.Ordinary;
    outputs.ndm.Full.results.lm_Rsquared_adj = lm_endm.Rsquared.Adjusted;
    
    flow = FlowCalculator(yopt,C,beta_opt,1,U,b_opt);
    for i = 1:size(flow,3)
        flow_ = flow(:,:,i);
        flow_(flow_ < prctile(nonzeros(flow),ipR.flowthresh)) = 0;
        flow(:,:,i) = flow_;
    end
    outputs.ndm.Full.flow = flow;    
    
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

        disp(['AIC = ' num2str(outputs.ndm.Full.results.lm_AIC)])
        disp(['BIC = ' num2str(outputs.ndm.Full.results.lm_BIC)])
        disp(['Intercept = ' num2str(outputs.ndm.Full.results.lm_intercept)])
        disp(['pValue = ' num2str(outputs.ndm.Full.results.lm_pval)])
        disp(['Rsqr_ord = ' num2str(outputs.ndm.Full.results.lm_Rsquared_ord)])
        disp(['Rsqr_adj = ' num2str(outputs.ndm.Full.results.lm_Rsquared_adj)])
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