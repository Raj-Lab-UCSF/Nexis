% NDM
%
% Evaluate NDM predictions [y] at desired time stamps.
%
% Input: 
%   x0                = Initial Condition
%   time_stamps = provided in units given by experiment
%   C                  = connectivity matrix
%   beta              = diffusivity parameter
%
% Output:
% y  =   NDM predicted vectors (columns) at each time stamp

function [y] = NDM(x0,time_stamps,C,beta)

% Define Laplacian matrix L
    rowdegree = (sum(C, 2)).';
    L = diag(rowdegree) - C;

    for j = 1:length(time_stamps)
        y(:,j) = expm(-beta*time_stamps(j)*L)*x0;
    end

end
