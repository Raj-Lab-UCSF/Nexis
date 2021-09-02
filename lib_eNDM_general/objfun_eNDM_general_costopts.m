% Objective function for eNDM_general
%
% param(1) = seed rescale factor
% param(2) = alpha
% param(3) = beta
% param(4:(n_types+3)) = a
% param((n_types+4):(2*n_types+3)) = b
% param((2*n_types+4):(3*n_types+3)) = p

function [f,newxt,newpath] = objfun_eNDM_general_costopts(param,seed_location,pathology,ts,C_,U_,solvetype_,volcorrect_,costfun)

LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);
n_types = size(U_,2);
x0_ = param(1)*seed_location;
alpha_ = param(2);
beta_ = param(3);
a_ = param(4:(n_types+3));
b_ = param((n_types+4):(2*n_types+3));
p_ = param((2*n_types+4):(3*n_types+3));

% Calculate predictions y with eNDM     
    % Solve eNDM; 
    [y] = eNDM_general(x0_,ts,C_,U_,alpha_,beta_,a_,b_,p_,solvetype_,volcorrect_);

% Modify quadratic error objfun to accomodate NaN    
if strcmp(costfun,'sse_sum')
    f = nansum(nansum((y - pathology).^2));
elseif strcmp(costfun,'rval_sum')
    Rvalues = zeros(1,length(ts));
    for jj = 1:length(ts)
        Rvalues(jj) = corr(y(:,jj),pathology(:,jj), 'rows','complete');
    end
    f = length(ts) - sum(Rvalues);
elseif strcmp(costfun,'sse_end')
    f = nansum(nansum((y(:,end) - pathology(:,end)).^2));
elseif strcmp(costfun,'rval_end')
    Rvalues = zeros(1,length(ts));
    for jj = 1:length(ts)
        Rvalues(jj) = corr(y(:,jj),pathology(:,jj), 'rows','complete');
    end
    f = 1 - Rvalues(end);
elseif strcmp(costfun,'LinR')
    Rvalues = zeros(1,length(ts));
    naninds = isnan(pathology(:,1));
    newxt = y; newxt(naninds,:) = [];
    newpath = pathology; newpath(naninds,:) = [];
    for jj = 1:length(ts)
        Rvalues(jj) = LinRcalc(newxt(:,jj),newpath(:,jj));
    end
    f = length(ts) - sum(Rvalues);
%     f = length(ts) - sum(Rvalues) + sum(abs(param))/length(param);
elseif strcmp(costfun,'LinR_end')
    Rvalues = zeros(1,length(ts));
    naninds = isnan(pathology(:,1));
    newxt = y; newxt(naninds,:) = [];
    newpath = pathology; newpath(naninds,:) = [];
    for jj = 1:length(ts)
        Rvalues(jj) = LinRcalc(newxt(:,jj),newpath(:,jj));
    end
    f = 1 - Rvalues(end);
end

end

