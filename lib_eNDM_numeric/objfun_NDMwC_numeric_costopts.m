% Objective function for NDMwC (using numeric integration for ode)
%
% param(1) = beta
% param(2) = x0_value
% param(3) = alpha0 = cte growth/decay term independent of x

function [f] = objfun_NDMwC_numeric_costopts(param,seed_location,pathology,time_stamps,C,alpha0,costfun)

beta = param(1);
x0 = param(2)*seed_location; 
alpha0 = param(3);

% Calculate predictions y with NDMwC     
    % Solve NDM; 
        [y] = NDMwC_numeric(x0,time_stamps,C,beta,alpha0);

% Modify quadratic error objfun to accomodate NaN    
if strcmp(costfun,'sse_sum')
    f = nansum(nansum((y - pathology).^2));
elseif strcmp(costfun,'rval_sum')
    for jj = 1:length(time_stamps)
        Rvalues(jj,:) = corr(y(:,jj),pathology(:,jj), 'rows','complete');
    end
    f = length(time_stamps) - sum(Rvalues);
elseif strcmp(costfun,'sse_end')
    f = nansum(nansum((y(:,end) - pathology(:,end)).^2));
elseif strcmp(costfun,'rval_end')
    for jj = 1:length(time_stamps)
        Rvalues(jj,:) = corr(y(:,jj),pathology(:,jj), 'rows','complete');
    end
    f = 1 - Rvalues(end);
elseif strcmp(costfun,'LinR')
    LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);
    naninds = isnan(pathology(:,1));
    newxt = y; newxt(naninds,:) = [];
    newpath = pathology; newpath(naninds,:) = [];
    for jj = 1:length(time_stamps)
        Rvalues(jj,:) = LinRcalc(newxt(:,jj),newpath(:,jj));
    end
    f = length(time_stamps) - sum(Rvalues);
end

end

