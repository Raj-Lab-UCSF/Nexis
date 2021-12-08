% NDM with Source
%
% Evaluate NDM predictions [y] at desired time stamps.
%
% Input: 
%   x0                = Initial Condition
%   time_stamps = provided in units given by experiment
%   C                  = connectivity matrix
%   beta             = diffusivity parameter
%   alpha1          = basal (constant) growth/clearance
%   alpha2          = linear growth/clearance term in x
%
% Output:
% y  =   NDM predicted vectors (columns) at each time stamp

function y = NDMwS(x0,time_stamps,C,beta,alpha1,alpha2)

% Define Laplacian matrix L
    rowdegree = (sum(C, 2)).';
    L = diag(rowdegree) - C;
    
% Define Diagonal matrix D
    D = eye(size(C,2));  

    for j = 1:length(time_stamps)
        y(:,j) = expm(  (-beta*L + alpha2*D )*time_stamps(j)  )*x0 + alpha1*ones(size(C,1),1)*time_stamps(j);
    end

end
