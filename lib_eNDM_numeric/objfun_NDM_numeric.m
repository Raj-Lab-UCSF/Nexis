% Objective function for NDM (using numeric integration for ode)
%
% param(1) = beta
% param(2) = x0_value

function [f] = objfun_NDM_numeric(param,seed_location,pathology,time_stamps,C)

beta = param(1);
x0 = param(2)*seed_location; 

% Calculate predictions y with NDM     
    % Solve NDM; 
        [y] = NDM_numeric(x0,time_stamps,C,beta);

% Modify quadratic error objfun to accomodate NaN    
f = nansum(nansum((y - pathology).^2));

