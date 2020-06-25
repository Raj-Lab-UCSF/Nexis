% NDM with Source
%
% Evaluate NDM predictions [y] at desired time stamps.
%
% Input: 
%   x0                = Initial Condition
%   time_stamps = provided in units given by experiment
%   C                  = connectivity matrix
%   beta             = diffusivity parameter
%   alpha1          = linear growth/clearance term in x
%
% Output:
% y  =   NDMwS predicted vectors (columns) at each time stamp

function y = NDMwS_analytic(x0,time_stamps,C,beta,alpha1)

% Define Laplacian matrix L
    rowdegree = (sum(C, 2)).';
    L = diag(rowdegree) - C;
    
% Define Diagonal matrix D
    D = eye(size(C,2));  

    for j = 1:length(time_stamps)
        y(:,j) = expm(  (-beta*L + alpha1*D )*time_stamps(j)  )*x0;
    end

end
