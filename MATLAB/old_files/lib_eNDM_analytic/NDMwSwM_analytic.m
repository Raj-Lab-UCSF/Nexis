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
%   alpha2          = microglia reweighting for alpha_1
%
% Output:
% y  =   NDM predicted vectors (columns) at each time stamp

function y = NDMwSwM(x0,time_stamps,C,u,beta,alpha1,alpha2)

% Define Laplacian matrix L
    rowdegree = (sum(C, 2)).';
    L = diag(rowdegree) - C;
    
% Define Diagonal matrix D
% This is where MICROGLIA u is added!
  D = eye(size(C,1));
  DM = diag(u);        
  
% Calculate predictions y with NDMwSwM         
    for j = 1:length(time_stamps)
        y(:,j) = expm(  (-beta*L + alpha1*D + alpha1*alpha2*DM)*time_stamps(j)  )*x0;
    end

end
