% Objective function for NDM 
%
% param(1) = beta
% param(2) = x0_value

function [f] = objfun_miss(param,pathology,time_stamps,C,beta,T)

stime_stamps = time_stamps - time_stamps(1);

% Define Laplacian matrix L
    rowdegree = (sum(C, 2)).';
    L = diag(rowdegree) - C;

% Calculate predictions y with NDM     
    for j = 1:length(time_stamps)
        y(:,j) = expm(-beta*time_stamps(j)*L)*(T*param );
    end

% Modify quadratic error objfun to accomodate NaN    
f = nansum(nansum((y - pathology).^2));


