% NDMwS
%
% Evaluate NDMwS predictions [y] at desired time stamps.
%
% Input: 
%   x0                = Initial Condition
%   time_stamps = provided in units given by experiment
%   C                  = connectivity matrix
%   beta             = diffusivity parameter
%   alpha0          = cte growth/decay term independent of x
%
% Output:
% y  =   NDMwS predicted vectors (columns) at each time stamp

function [y] = NDMwC_numeric(x0,time_stamps,C,beta,alpha0)
   
    % Define the RHS of your ODE
    function dydt = ode(t,x,C,beta,alpha0)
        
        % Define Laplacian matrix L
        rowdegree = (sum(C, 2)).';
        L = diag(rowdegree) - C;
       
        dydt = -beta*L*x + alpha0*ones(size(x));        
    end

% WITH OPTIONS
%opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
%[t,y] = ode45(@ode,time_stamps,x0,opts,C,beta);    

% NO OPTIONS
[t,y] = ode45(@ode,time_stamps,x0,[],C,beta,alpha0);    


y = y'; 
%disp(t)
end
