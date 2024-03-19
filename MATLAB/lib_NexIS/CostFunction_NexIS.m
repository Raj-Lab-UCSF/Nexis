% Objective function for NexIS_fun
%
% param(1) = gamma (seed rescale factor)
% param(2) = alpha
% param(3) = beta
% param(4) = s
% param(5:(n_types+4)) = b
% param((n_types+5):(2*n_types+4)) = p

function f = CostFunction_NexIS(params_,C_,U_,time_stamps_,seed_,pathol_,...
    solvetype_,volcorrect_,costfun_,excltpts_costfun_,exclseed_costfun_,...
    use_dataspace_,studyname_,logtrans_,lambdaval_,matdir_)

LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

ts = time_stamps_;
pathology = pathol_;
seedregs = seed_;
s_ = params_(4);

% Calculate predictions y with NexIS
predicted = NexIS_fun(C_,U_,ts,seedregs,params_,solvetype_,volcorrect_,matdir_);
if use_dataspace_ && ~strcmp(studyname_,'User-specified') 
    predicted = CCFToData(predicted,studyname_,matdir_);
    pathology = CCFToData(pathology,studyname_,matdir_);
end

% Exclude selected time points from cost function
ts(excltpts_costfun_) = [];
pathology(:,excltpts_costfun_) = [];
predicted(:,excltpts_costfun_) = [];

% Remove seed regions from cost function
if logical(exclseed_costfun_) && ~isequal(seedregs,NaN)
    seedbin = logical(seedregs);
    predicted(seedbin,:) = [];
    pathology(seedbin,:) = [];
elseif logical(exclseed_costfun_) && isequal(seedregs,NaN)
    warning('Seed region undefined - disregarding seed exclusion');
end

% Remove NaN regions, if any - this is written in such a way to account for
% NaNs occuring at different places at different times
if any(isnan(pathology(:)))
    timeinds = [];
    pathology_vec = [];
    predicted_vec = [];
    for i = 1:length(ts)
        timeinds_i = i * ones(size(pathology,1),1);
        pathology_i = pathology(:,i);
        predicted_i = predicted(:,i);
        naninds_i = isnan(pathology_i);
        timeinds_i(naninds_i) = [];
        pathology_i(naninds_i) = [];
        predicted_i(naninds_i) = [];
        timeinds = [timeinds; timeinds_i];
        pathology_vec = [pathology_vec; pathology_i];
        predicted_vec = [predicted_vec; predicted_i];
    end
else
    timeinds = zeros(size(pathology,1),length(ts));
    for i = 1:length(ts)
        timeinds(:,i) = i * ones(size(pathology,1),1);
    end
    timeinds = timeinds(:);
    pathology_vec = pathology(:);
    predicted_vec = predicted(:);
end

% Perform optional log transformation
if strcmp(logtrans_,'log')
    pathology_vec = log(pathology_vec);
    predicted_vec = log(predicted_vec);
    % It is necessary to remove all 0's, which are -Inf after log()
    pathology_vec(~isfinite(pathology_vec)) = [];
    predicted_vec(~isfinite(pathology_vec)) = [];
    timeinds(~isfinite(pathology_vec)) = [];
    pathology_vec(~isfinite(predicted_vec)) = [];
    predicted_vec(~isfinite(predicted_vec)) = [];
    timeinds(~isfinite(predicted_vec)) = [];
elseif strcmp(logtrans_,'log+1')
    pathology_vec = log(pathology_vec + 1);
    predicted_vec = log(predicted_vec + 1);
end

% Define cost function
switch costfun_
    case 'sse_sum'
        f = sum((predicted_vec - pathology_vec).^2);
    case 'sse_end'
        end_ind = max(timeinds);
        f = sum((predicted_vec(timeinds == end_ind) - pathology(timeinds == end_ind)).^2);
    case 'rval'
        f = 1 - corr(predicted_vec,pathology_vec);
    case 'rval_sum'
        Rvalues = zeros(1,length(ts));
        for i = 1:length(ts)
            Rvalues(i) = corr(predicted_vec(timeinds == i),pathology_vec(timeinds == i));
        end
        f = length(ts) - sum(Rvalues);
    case 'rval_end'
        end_ind = max(timeinds);
        f = 1 - corr(predicted_vec(timeinds == end_ind),pathology_vec(timeinds == end_ind));
    case 'linr'
        f = 1 - LinRcalc(predicted_vec,pathology_vec);
    case 'linr_reg_s'
        f = 1 - LinRcalc(predicted_vec,pathology_vec) + lambdaval_*abs(s_ - 0.5);
    case 'linr_sum'
        Rvalues = zeros(1,length(ts));
        for i = 1:length(ts)
            Rvalues(i) = LinRcalc(predicted_vec(timeinds == i),pathology_vec(timeinds == i));
        end
        f = length(ts) - sum(Rvalues);
    case 'linr_end'
        end_ind = max(timeinds);
        f = 1 - LinRcalc(predicted_vec(timeinds == end_ind),pathology_vec(timeinds == end_ind));
end

end

