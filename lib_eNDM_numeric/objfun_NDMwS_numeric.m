% Objective function for NDM (using numeric integration for ode)
%
% param(1) = beta
% param(2) = x0_value
% param(3) = alpha1

function [f] = objfun_NDMwS_numeric(param,seed_location,pathology,time_stamps,C,alpha1)

beta = param(1);
x0 = param(2)*seed_location; 
alpha1 = param(3);

% Calculate predictions y with NDM     
    % Solve NDM; 
        [y] = NDMwS_numeric(x0,time_stamps,C,beta,alpha1);

% Modify quadratic error objfun to accomodate NaN    
f = nansum(nansum((y - pathology).^2));

