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

function [y] = NDM_numeric(x0,time_stamps,C,beta)
   
    % Define the RHS of your ODE
    function dydt = ode(t,x,C,beta)
        
        % Define Laplacian matrix L
        rowdegree = (sum(C, 2)).';
        L = diag(rowdegree) - C;
    
        dydt = -beta*L*x;        
    end

% WITH OPTIONS
%opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
%[t,y] = ode45(@ode,time_stamps,x0,opts,C,beta);    

% NO OPTIONS
[t,y] = ode15s(@ode,time_stamps,x0,[],C,beta);    


y = y'; 
%disp(t)
end
