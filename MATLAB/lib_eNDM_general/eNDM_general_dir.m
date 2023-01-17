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
%   alpha           = global amplification/depletion term  (scalar)
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

function [y,A] = eNDM_general_dir(x0,time_stamps,C,U,alpha,beta,s,a,b,p,solvetype,volcorrect)
if nargin < 12
    volcorrect = 0;
    if nargin < 11
        solvetype = 'numeric';
    end
end

% Define source vector s_a (n_ROI x 1)
if size(a,2) > size(a,1)
    a = a.';
end
s_a = U * a; 

% Define Diagonal matrix Gamma
if size(p,2) > size(p,1)
    p = p.';
end
s_p = U * p; 
Gamma = diag(alpha + s_p); %% Yuanxi's Comment: if s_p is an effect matrix, why not alpha*s_p?

% Define Laplacian matrix L
C_dir = (1-s)*C.' + s*C;
coldegree = (sum(C_dir, 1));
L_raw = diag(coldegree) - C_dir;
if size(b,2) > size(b,1)
    b = b.';
end
s_b = U * b;
S_b = repmat(s_b,1,length(s_b)) + ones(length(s_b));
L = L_raw;
if logical(volcorrect)
    load([cd filesep 'raw_data_mouse' filesep 'regionvoxels.mat'],'voxels');
    voxels_2hem = cat(1,voxels,voxels)/2; % approximation, splitting over the two hemispheres
    L = mean(voxels_2hem) * diag(voxels_2hem.^(-1)) * L; % correction proposed by Putra et al. 2021
end

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
