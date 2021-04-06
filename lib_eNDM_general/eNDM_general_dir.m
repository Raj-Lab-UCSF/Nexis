% eNDM General w/ Directionality
%
% Evaluate eNDM General predictions [y] at desired time stamps.
%
% Input: 
%   x0              = Initial Condition (n_ROI x 1)
%   time_stamps     = provided in units given by experiment (1 x n_t)
%   C               = connectivity matrix (n_ROI x n_ROI) - SHOULD BE
%                               RETROGRADE-DIRECTED (c(1,2) = connectivity
%                               from ROI_1 to ROI_2)
%   U               = cell types matrix (n_ROI x n_types)
%   alpha           = global amplification/depletion term (scalar)
%   a               = source term as a function of cell types (n_types x 1)
%   beta            = global diffusivity parameter (scalar)
%   s               = global directionality parameter (scalar between 0 and
%                               1, with 1 = fully retrograde spread given
%                               the orientation of C specified above)
%   b               = effect of cell types on transmission (n_types x 1)
%   p               = effect of types on accumulation/depletion (n_types x 1)
%
% Output:
% y  =   eNDM predicted vectors (columns) at each time stamp

function [y] = eNDM_general_dir(x0,time_stamps,C,U,alpha,beta,s,a,b,p,solvetype)
if nargin < 11
    solvetype = 'numeric';
end

% Define source vector s_a (n_ROI x 1)
s_a = U * a;

% Define Diagonal matrix Gamma
s_p = U * p; 
Gamma = diag(alpha + s_p);

% Define Laplacian matrix L
C_dir = (1-s)*C.' + s*C;
coldegree = (sum(C_dir, 1));
L_raw = diag(coldegree) - C_dir;
s_b = U * b;
S_b = repmat(s_b,1,length(s_b)) + ones(length(s_b));
L = L_raw .* S_b.';

% define system dydt = Ax + B
A = Gamma - beta * L;
B = s_a;     

% solve system using analytical or numerical methods
if all(B == 0) && strcmp(solvetype,'analytic')
    y = zeros(size(C,1),length(time_stamps));
    y_analytic = @(t,x0) expm(A*t) * x0;
    for i = 1:length(time_stamps)
        y(:,i) = y_analytic(time_stamps(i),x0);
    end
else
    [~,y] = ode45(@(t,x) (A*x + B),time_stamps,x0);    
    y = y'; 
end
end
