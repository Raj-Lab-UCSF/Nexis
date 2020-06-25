% NDMwS
%
% Evaluate NDMwS predictions [y] at desired time stamps.
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

function [y] = NDMwS_numeric(x0,time_stamps,C,beta,alpha1)
   
    % Define the RHS of your ODE
    function dydt = ode(t,x,C,beta,alpha1)
        
        % Define Laplacian matrix L
        rowdegree = (sum(C, 2)).';
        L = diag(rowdegree) - C;
        
        % Define Diagonal matrix D
        D = eye(size(C,2));  
        dydt = (-beta*L + alpha1*D)*x;         
        %    dydt = -beta*L*x + alpha0*ones(size(x));        
    end

% WITH OPTIONS
%opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
%[t,y] = ode45(@ode,time_stamps,x0,opts,C,beta);    

% NO OPTIONS
[t,y] = ode45(@ode,time_stamps,x0,[],C,beta,alpha1);    


y = y'; 
%disp(t)
end
