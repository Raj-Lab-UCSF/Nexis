% Objective function for NDM with Source with Microglia (NDMwSwM)
%
% param(1) = beta
% param(2) = x0_value
% param(3) = alpha1 = linear growth/clearance term in x
% param(4) = alpha2 = microglia reweighting for alpha1  

function [f] = objfun_NDMwSwM(param,seed_location,pathology,time_stamps,C,u)

% Define Laplacian matrix L
    rowdegree = (sum(C, 2)).';
    L = diag(rowdegree) - C;

% Define Diagonal matrix D 
% This is where MICROGLIA u is added!
    D = eye(size(C,1));
    DM = diag(u);    
    
% Calculate predictions y with NDM     
    for j = 1:length(time_stamps)         
        y(:,j) = expm(  (-param(1)*L + param(3)*D + param(3)*param(4)*DM )*time_stamps(j)  )*seed_location*param(2);
    end

% Modify quadratic error objfun to accomodate NaN    
f = nansum(nansum((y - pathology).^2));

