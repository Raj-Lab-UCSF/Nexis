% Objective function for extended NDM  (eNDM) 
%
% param(1) = x0_value
% param(2) = beta1 = diffusivity parameter
% param(3) = alpha1 = linear growth/clearance term in x
% param(4) = alpha2 = cell-type reweighting for alpha1
% param(5) = beta2 =  cell-type reweighting for beta1

function [f] = objfun_eNDM_costopts(param,seed_location,pathology,time_stamps,C,u,costfun)

% Define Diagonal matrix D 
% This is where MICROGLIA u is added!
    D = eye(size(C,1));
    DM = diag(u);    
 
 % Define Laplacian matrix L
    rowdegree = (sum(C, 2)).';
    L = 0*C; %diag(rowdegree) - C;

    for i = 1:size(C,1)
        for j = 1:size(C,1)
            if i ~= j
            L(i,j) = - C(i,j)*(1 + u(i)*param(5));
            else
            L(i,j) = rowdegree(i)*(1 + u(i)*param(5));    
            end
        end        
    end
    
% Calculate predictions y with eNDM     
    for j = 1:length(time_stamps)         
        y(:,j) = expm(  (-param(2)*L+ param(3)*D + param(3)*param(4)*DM )) *time_stamps(j) *seed_location*param(1);
    end    
    
% % Calculate predictions y with eNDM %     for j = 1:length(time_stamps) %
% y(:,j) = expm(  (-param(2)*L +param(2)*param(5)*DM*L+ param(3)*D +
% param(3)*param(4)*DM )) *time_stamps(j) *seed_location*param(1); %
% end

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


