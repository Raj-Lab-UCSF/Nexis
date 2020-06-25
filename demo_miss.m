%% Demo: MISS  

%%  0.  Setup   

% Initialization
    clear all; close all; clc; 

% Add folder with raw data to path  
    addpath([pwd '/raw_data_mouse'])

% Add library with eNDM functions to path
    addpath([pwd '/lib_eNDM_numeric'])
    addpath([pwd '/lib_eNDM_analytic'])

% Load dataset of interest
   load celltypedata_mrx3_540.mat
   
%   load MouseDataForPedro.mat
   load Hurtado426.mat

% Select connectome C (426x426)
    %C = Networks.nd;        % symmmetric  
    C = Networks.ret;      % non-symmetric  
    %C = Networks.ant;     % non-symmetric      

% Normalize C
    cmax = max(max(C));
    cmin = min(min(C)); 
    C = (C - cmin)./(cmax-cmin);
        
% Define number of regions of interest (nroi) 
    nroi = size(C,1);

% Load time stamps, pathology measurements, and the seed_location 
   
        % Hurtado   
           time_stamps = [2,4,6,8];
           pathology = data426.Hurtado  ;
           
        % Hippocampus Injection   
%            time_stamps = [1,3,6];
%            pathology = data426.IbaHippInj;

        % Striatal Injection 
%            time_stamps = [1,3,6];
%            pathology = data426.IbaStrInj;
           

%% 

% Normalize every column of cell type matrix T 
    for j = 1:25    
        u = celltypedata(:,j); 
        u = (u - min(u))./(max(u) - min(u));
        T(:,j) = u;
    end

   % Init Guess
   
   init_guess_params = zeros(size(T,2),1); 
   lb = 0*ones(size(T,2),1);  
   ub = 10*ones(size(T,2),1);
    
    result_obj = [];
    result_param = [];
% Span beta values    
    
    for beta = 0:.5:3
        
    % Apply fmincon
    [param2, fval] = fmincon(@(param)objfun_miss(param,pathology,time_stamps,C,beta,T),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
     
                       
     result_obj = [result_obj;fval];
     result_param = [result_param,param2];
    end
    
    figure
    plot(0:.5:3,result_obj)

corr(T*result_param(:,1),pathology(:,1),'Rows','complete')

corr(T*result_param(:,1),pathology(:,2),'Rows','complete')

corr(T*result_param(:,1),pathology(:,3),'Rows','complete')

corr(T*result_param(:,1),pathology(:,4),'Rows','complete')
    
