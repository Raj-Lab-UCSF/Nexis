% Objective function for NDM 
%
% param(1) = beta
% param(2) = x0_value

function [f] = objfun_NDM(param,seed_location,pathology,time_stamps,C)

%y = param(1)*param(2) - pathology(2,1) + C(1,1) + seed_location(1);

% Define Laplacian matrix L
    rowdegree = (sum(C, 2)).';
    L = diag(rowdegree) - C;

% Calculate predictions y with NDM     
    for j = 1:length(time_stamps)
        y(:,j) = expm(-param(1)*time_stamps(j)*L)*seed_location*param(2);
    end

% Modify quadratic error objfun to accomodate NaN    
f = nansum(nansum((y - pathology).^2));

