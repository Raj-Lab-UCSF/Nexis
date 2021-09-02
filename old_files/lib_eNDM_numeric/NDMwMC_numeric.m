% NDMwS
%
% Evaluate NDMwMC predictions [y] at desired time stamps.
%
% Input: 
%   x0                = Initial Condition
%   time_stamps = provided in units given by experiment
%   C                  = connectivity matrix
%   beta             = diffusivity parameter
%   alpha0          = cte growth/decay term depending on u
%
% Output:
% y  =   NDMwMC predicted vectors (columns) at each time stamp

function [y] = NDMwMC_numeric(x0,time_stamps,C,u,beta,alpha0)
   
    % Define the RHS of your ODE
    function dydt = ode(t,x,C,u,beta,alpha0)
        
        % Define Laplacian matrix L
        rowdegree = (sum(C, 2)).';
        L = diag(rowdegree) - C;
        
        % Define Diagonal matrix D
        D = eye(size(C,2));  
        
        dydt = -beta*L*x + alpha0*u;        
    end

% WITH OPTIONS
%opts = odeset('RelTol',1e-5,'AbsTol',1e-10);

% NO OPTIONS
[t,y] = ode45(@ode,time_stamps,x0,[],C,u,beta,alpha0);    


y = y'; 
%disp(t)
end
