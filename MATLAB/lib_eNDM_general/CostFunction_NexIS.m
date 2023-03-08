% Objective function for NexIS_fun
%
% param(1) = gamma (seed rescale factor)
% param(2) = alpha
% param(3) = beta
% param(4) = s
% param(5:(n_types+4)) = b
% param((n_types+5):(2*n_types+4)) = p

function [f,newxt,newpath] = CostFunction_NexIS(study_,U_,params_,...
    solvetype_,volcorrect_,costfun_,excltpts_costfun_,exclseed_costfun_,...
    logtrans_,matdir_)

LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

% Load appropriate data struct, load correct connectome
taustudies = {'IbaP301S','IbaHippInj','IbaStrInj','BoludaDSAD','BoludaCBD'...
    'DS1','DS4','DS6','DS6_110','DS7','DS9','DS9_110','DS10'};
asynstudies = {'asyn_mouse','asyn_human','GCI','PFF','Henderson'};

if ismember(study_,taustudies)
    load([matdir_ filesep "Mouse_Tauopathy_Data_HigherQ.mat"], 'mousedata_struct');
elseif ismember(study_,asynstudies)
    load([matdir_ filesep "Mouse_Synuclein_Data.mat"], 'mousedata_struct');
else
    error('Provided study name is invalid.');
end

ts = mousedata_struct.(study_).time_stamps;
pathology = mousedata_struct.(study_).data;
seedregs = mousedata_struct.(study_).seed;

% Calculate predictions y with NexIS
predicted_ccf = NexIS_fun(study_,U_,params_,solvetype_,volcorrect_,matdir_);
predicted = CCFToData(predicted_ccf,study_,matdir_);

% Exclude selected time points from cost function
ts(excltpts_costfun_) = [];
pathology(:,excltpts_costfun_) = [];
predicted(:,exclseed_costfun_) = [];

% Remove seed regions from cost function
if logical(exclseed_costfun_) && ~isequal(seedregs,NaN)
    seedbin = logical(seedregs);
    predicted(seedbin,:) = [];
    pathology(seedbin,:) = [];
elseif logical(exclseed_costfun_) && isequal(seedregs,NaN)
    warning('Seed region undefined - disregarding seed exclusion');
end

% Remove NaN regions, if any
if any(isnan(pathology))
    for i = 1:
    timeinds = 
pathology(isnan(pathology(:,1)),:) = [];
predicted(isnan(pathology(:,1)),:) = [];

% Perform optional log transformation
if strcmp(logtrans_,'log')
    finiteinds_path = isfinite(sum(log(pathology),2));
elseif strcmp(logtrans_,'log+1')
    
end

% Define cost function
switch costfun_
    case 'sse_sum'
        f = sum(sum((predicted - pathology).^2));
    case 'sse_end'
        f = sum(sum((predicted(:,end) - pathology(:,end)).^2));
    case 'rval'
        f = 1 - corr(predicted(:),pathology(:));
    case 'rval_sum'
        Rvalues = zeros(1,length(ts));
        for i = 1:length(ts)
            Rvalues(i) = corr(predicted(:,i),pathology(:,i));
        end
        f = length(ts) - sum(Rvalues);
end

if strcmp(costfun_,'sse_sum')
    
elseif strcmp(costfun_,'rval_sum')
    Rvalues = zeros(1,length(ts));
    for i = 1:length(ts)
        Rvalues(i) = corr(predicted(:,i),pathology(:,i), 'rows','complete');
    end
    f = length(ts) - sum(Rvalues);
elseif strcmp(costfun_,'sse_end')
    f = nansum(nansum((predicted(:,end) - pathology(:,end)).^2));
elseif strcmp(costfun_,'rval_end')
    Rvalues = zeros(1,length(ts));
    for i = 1:length(ts)
        Rvalues(i) = corr(predicted(:,i),pathology(:,i), 'rows','complete');
    end
    f = 1 - Rvalues(end);
elseif strcmp(costfun_,'LinR')
    Rvalues = zeros(1,length(ts));
%     naninds = isnan(pathology(:,1)); % Yuanxi's Comment: Our dataset has some missing data. It could occur at different MPI. 
%                                                             % I recommend that here could be revised to 'naninds = isnan(sum(pathology,2));'
    naninds = isnan(prod(pathology,2));
    newxt = predicted; newxt(naninds,:) = [];
    newpath = pathology; newpath(naninds,:) = [];
    for i = 1:length(ts)
        Rvalues(i) = LinRcalc(newxt(:,i),newpath(:,i));
    end
    f = length(ts) - sum(Rvalues);
%     f = length(ts) - sum(Rvalues) + sum(abs(param))/length(param);
elseif strcmp(costfun_,'LinR_end')
    Rvalues = zeros(1,length(ts));
    naninds = isnan(prod(pathology,2));
    newxt = predicted; newxt(naninds,:) = [];
    newpath = pathology; newpath(naninds,:) = [];
    for i = 1:length(ts)
        Rvalues(i) = LinRcalc(newxt(:,i),newpath(:,i));
    end
    f = 1 - Rvalues(end);
elseif strcmp(costfun_,'log_rval_sum')
    Rvalues = zeros(1,length(ts));
    finiteinds = isfinite(sum(log(pathology),2));
    newxt = predicted;
    for i = 1:length(ts)
        Rvalues(i) = corr(log(predicted(finiteinds,i)),log(pathology(finiteinds,i)), 'rows','complete');
    end
    f = length(ts) - sum(Rvalues);

elseif strcmp(costfun_,'log_rval_end')
    Rvalues = zeros(1,length(ts));
    finiteinds = isfinite(sum(log(pathology),2));
    newxt = predicted;
    for i = 1:length(ts)
        Rvalues(i) = corr(log(predicted(finiteinds,i)),log(pathology(finiteinds,i)), 'rows','complete');
    end
    f = 1 - Rvalues(end);

elseif strcmp(costfun_,'log_LinR_sum')
    Rvalues = zeros(1,length(ts));
    finiteinds = isfinite(sum(log(pathology),2));
    newxt = predicted;
    for i = 1:length(ts)
        Rvalues(i) = LinRcalc(log(predicted(finiteinds,i)),log(pathology(finiteinds,i)));
    end
    f = length(ts) - sum(Rvalues);
end
% fprintf('f = %d\n',f)
% display(seed_location.');
% display(y);
% display(pathology);

end

