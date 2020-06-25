% Objective function for NDM with Source (NDMwS)
%
% param(1) = beta
% param(2) = x0_value
% param(3) = alpha1 = basal (constant) growth/clearance
% param(4) = alpha2 = linear growth/clearance term in x


function [f] = objfun_NDMwS(param,seed_location,pathology,time_stamps,C)

% Define Laplacian matrix L
    rowdegree = (sum(C, 2)).';
    L = diag(rowdegree) - C;

% Define Diagonal matrix D
    D = eye(size(C,2));    
    
% Calculate predictions y with NDM     
    for j = 1:length(time_stamps)
         
        y(:,j) = expm(  (-param(1)*L + param(4)*D )*time_stamps(j)  )*seed_location*param(2) + param(3)*ones(size(C,1),1)*time_stamps(j);

        %y(:,j) = expm(-param(1)*time_stamps(j)*L)*seed_location*param(2);
    end

% Modify quadratic error objfun to accomodate NaN    
f = nansum(nansum((y - pathology).^2));

