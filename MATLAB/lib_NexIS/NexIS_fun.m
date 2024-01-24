% NexIS Model 
%
% Evaluate NexIS predictions [y] at desired time stamps.
%
% Input: 
%   study_           = study name (char)
%   indiv_           = individual mouse rescale flag (logical)
%   U_               = cell types matrix (n_ROI x n_types)
%   params_          = parameter vector (1 x 6 corresponding to [gamma,
%                               beta, alpha, s, b, p])
%   solvetype_       = 'analytic' or 'numeric' (char)
%   volcorrect_      = volume correction flag (recommended) (logical)
%   matdir_          = directory to load data from (char)
%
% Output:
% y  =  NexIS predicted vectors (columns) at each time stamp

function [y,A] = NexIS_fun(C_,U_,time_stamps_,x0,params_,solvetype_,volcorrect_,matdir_)
if nargin < 8
    matdir_ = [cd filesep 'raw_data_mouse'];
    if nargin < 7
        volcorrect_ = 1;
        if nargin < 6
            solvetype_ = 'numeric';
        end
    end
end

% Define parameters and model-specific conditions
% if ~iterative_
    n_types = (length(params_) - 4)/2;
    gamma = params_(1);
    alpha = params_(2);
    beta = params_(3);
    s = params_(4);
    b = params_(5:(n_types+4));
    p = params_((n_types+5):(2*n_types+4));
    if isequal(unique(x0),[0;1]) % if binary seed location and not using baseline
        x0 = gamma*x0;
    end

% else
% % Fill in later
% end

% % Define source vector s_a (n_ROI x 1)
% if size(a,2) > size(a,1)
%     a = a.';
% end
% s_a = U_ * a; 

% Define Diagonal matrix Gamma
if size(p,2) > size(p,1)
    p = p.';
end
s_p = U_ * p; 
Gamma = diag(alpha + s_p);

% Define Laplacian matrix L
C_dir = (1-s)*C_.' + s*C_;
coldegree = (sum(C_dir, 1));
L_raw = diag(coldegree) - C_dir;
if size(b,2) > size(b,1)
    b = b.';
end
s_b = U_ * b;
S_b = repmat(s_b,1,length(s_b)) + ones(length(s_b));
L = L_raw .* S_b.'; 
if logical(volcorrect_)
    load([matdir_ filesep 'DefaultAtlas.mat'], 'DefaultAtlas');
    voxels_2hem = DefaultAtlas.volumes;
else
    voxels_2hem = ones(size(L,1),1);
end
L = mean(voxels_2hem) * diag(voxels_2hem.^(-1)) * L; % correction proposed by Putra et al. 2021

% define system dydt = Ax
A = Gamma - beta * L;   

% solve system using analytical or numerical methods
if strcmp(solvetype_,'analytic')
    y = zeros(size(C_,1),length(time_stamps_));
    y_analytic = @(t,x0) expm(A*t) * x0;
    for i = 1:length(time_stamps_)
        y(:,i) = y_analytic(time_stamps_(i),x0);
    end
else
    [~,y] = ode45(@(t,x) (A*x + B),time_stamps_,x0);    
    y = y'; 
end
end
