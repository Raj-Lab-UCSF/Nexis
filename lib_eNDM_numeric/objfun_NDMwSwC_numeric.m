% Objective function for NDMwSwC (using numeric integration for ode)
%
% param(1) = beta
% param(2) = x0_value
% param(3) = alpha0 = cte growth/decay term independent of x
% param(4) = alpha1

function [f] = objfun_NDMwSwC_numeric(param,seed_location,pathology,time_stamps,C)

beta = param(1);
x0 = param(2)*seed_location; 
alpha0 = param(3);
alpha1 = param(4);

% Calculate predictions y with NDMwC     
    % Solve NDM; 
        [y] = NDMwSwC_numeric(x0,time_stamps,C,beta,alpha0,alpha1);

% Modify quadratic error objfun to accomodate NaN    
f = nansum(nansum((y - pathology).^2));

