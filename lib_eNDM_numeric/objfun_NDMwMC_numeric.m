% Objective function for NDMwMC (using numeric integration for ode)
%
% param(1) = beta
% param(2) = x0_value
% param(3) = alpha0 = cte growth/decay term on u

function [f] = objfun_NDMwMC_numeric(param,seed_location,pathology,time_stamps,C,u)

beta = param(1);
x0 = param(2)*seed_location; 
alpha0 = param(3);

% Calculate predictions y with NDMwC     
    % Solve NDM; 
        [y] = NDMwMC_numeric(x0,time_stamps,C,u,beta,alpha0);

% Modify quadratic error objfun to accomodate NaN    
f = nansum(nansum((y - pathology).^2));

