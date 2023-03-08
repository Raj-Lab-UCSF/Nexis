% eNDM 
%
% Evaluate NDM predictions [y] at desired time stamps.
%
% Input: 
%   x0                = Initial Condition
%   time_stamps = provided in units given by experiment
%   C                  = connectivity matrix
%   beta1           = diffusivity parameter
%   alpha1          = linear growth/clearance term in x
%   alpha2          = microglia reweighting for alpha1
%   beta2            = microglia reweighting for beta1

% Output:
% y  =   NDM predicted vectors (columns) at each time stamp

function y = eNDM_analytic(x0,time_stamps,C,u,beta1,alpha1,alpha2,beta2)

% Define Diagonal matrix D 
% This is where MICROGLIA u is added!
    D = eye(size(C,1));
    DM = diag(u);    
 
 % Define Laplacian matrix L
    rowdegree = (sum(C, 2)).';
    L = 0*C; %diag(rowdegree) - C;

    for i = 1:size(C,1)
        for j = 1:size(C,1)
            if i ~= j
            L(i,j) = - C(i,j)*(1 + u(i)*beta2);
            else
            L(i,j) = rowdegree(i)*(1 + u(i)*beta2);    
            end
        end        
    end
    
% Calculate predictions y with eNDM     
    for j = 1:length(time_stamps)         
        y(:,j) = expm(  (-beta1*L+ alpha1*D + alpha1*alpha2*DM )) *time_stamps(j) *x0;
    end    
 
 

end
