%% Demo [Analytics vs Numeric]
%
% The goal of this demo is to apply eNDM to mouse data. 
% Uses numerical solutions for models using ode solvers in objective functions. 
% Several consistency checks are provided. 
%
% Written by: Pedro D. Maia
%
% Created: Feb/6/2020
%
% Last Modified: Feb/19/2020

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
   load MouseDataForPedro.mat

% Select connectome C (426x426)
    C = Networks.nd;        % symmmetric  
    %C = Networks.ret;      % non-symmetric  
    %C = Networks.ant;     % non-symmetric      

% Normalize C
    cmax = max(max(C));
    cmin = min(min(C)); 
    C = (C - cmin)./(cmax-cmin);
        
% Define number of regions of interest (nroi) 
    nroi = size(C,1);

% Load time stamps, pathology measurements, and the seed_location 
   
        % Hippocampus Injection   
           time_stamps = [1,3,6];
           pathology = data426.IbaHippInj;
           
           pathology = pathology./norm(pathology,2);
           seed_location = seed426.IbaHippInj;  

        % Striatal Injection 
%            time_stamps = [1,3,6];
%            pathology = data426.IbaStrInj;
%            seed_location = seed426.IbaStrInj;
           
%% 1.  [Analytics vs Numeric] NDM, fmincon minimizing MSE  

% Set ndm_comp = 1 to run
ndm_comp = 0;
if ndm_comp==1

% Parameters to fit
    % param(1) = beta
    % param(2) = x0_value

% Extra input to required by objective function
    % seed_location 
    % pathology
    % time_stamps
    % C
    
% Set initial guess for beta
    init_guess_params(1) = 1;

% Set initial guess for x0_value
    init_guess_params(2) = nansum(pathology(:,1));

% Set lower bounds (lb) for parameters
    lb = [0,0];

% Set upper bounds (ub) for parameters
%    ub = [10,nansum(pathology(:,1))];
ub = [3,3];
% Apply fmincon w/ numeric 
    [param_num, fval_num] = fmincon(@(param)objfun_NDM_numeric(param,seed_location,pathology,time_stamps,C),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDM with the optimal parameters
    beta_num = param_num(1);
    x0_num = seed_location*param_num(2);
    ynum = NDM_numeric(x0_num,time_stamps,C,beta_num);
    
    % Apply fmincon w/ analytic
    [param_ana, fval_ana] = fmincon(@(param)objfun_NDM_analytic(param,seed_location,pathology,time_stamps,C),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDM with the optimal parameters
    beta_ana = param_ana(1);
    x0_ana = seed_location*param_ana(2);
    yana = NDM_analytic(x0_ana,time_stamps,C,beta_ana);
    
    %DEBUG
    %beta_ana = param_num(1);
    %x0_ana = seed_location*param_num(2);
    %yana = NDM_analytic(x0_ana,time_stamps,C,beta_ana);

 %% Display Numeric Results
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(ynum(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDM minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['NUMERIC Beta = ' num2str(param_num(1)) ', x0 value = ' num2str(param_num(2))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval_num)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(ynum,pathology,time_stamps)

%% Display Analytic Results
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(yana(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDM minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['ANALYTIC Beta = ' num2str(param_ana(1)) ', x0 value = ' num2str(param_ana(2))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval_ana)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(yana,pathology,time_stamps)
clearvars y 

%% Compare values
%table(ynum,yana)
figure
subplot(3,1,1)
plot(ynum(:,1),yana(:,1),'o')
hold on
line([0 3], [0 3],'Color','k')
%
subplot(3,1,2)
plot(ynum(:,2),yana(:,2),'o')
hold on
line([0 3], [0 3],'Color','k')
%
subplot(3,1,3)
plot(ynum(:,3),yana(:,3),'o')
hold on
line([0 3], [0 3],'Color','k')
%
suptitle('Comparison: Numeric vs Analytic') 

disp('rel error at t1') 
disp(    norm(ynum(:,1)-yana(:,1),2) /norm(ynum(:,1),2)      )

disp('rel error at t2') 
disp(    norm(ynum(:,2)-yana(:,2),2) /norm(ynum(:,2),2)      )


disp('rel error at t3') 
disp(    norm(ynum(:,3)-yana(:,3),2) /norm(ynum(:,3),2)      )

end

%% 2.  [Analytics vs Numeric] NDMwS, fmincon minimizing MSE 

% Set ndmws_comp = 1 to run
ndmws_comp = 1;

if ndmws_comp ==1

% Parameters to fit
    % param(1) = beta
    % param(2) = x0_value
    % param(3) = alpha1 = linear growth/clearance term in x

% Extra input to required by objective function
    % seed_location 
    % pathology
    % time_stamps
    % C
    
% Set initial guess for beta
    init_guess_params(1) = 1;

% Set initial guess for x0_value
    init_guess_params(2) = nansum(pathology(:,1));

% Set initial guess for alpha1 and alpha2
    init_guess_params(3) = .5;
    
% Set lower bounds (lb) for parameters
    lb = [0,0,-5];

% Set upper bounds (ub) for parameters
%    ub = [10,nansum(pathology(:,end)),5];
    ub = [3,3,3];

% Apply fmincon w/ numeric 
    [param_num, fval_num] = fmincon(@(param)objfun_NDMwS_numeric(param,seed_location,pathology,time_stamps,C),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDMwS with the optimal parameters
    beta_num = param_num(1);
    x0_num = seed_location*param_num(2);
    alpha1_num = param_num(3);
    ynum = NDMwS_numeric(x0_num,time_stamps,C,beta_num,alpha1_num);
    
    % Apply fmincon w/ analytic
    [param_ana, fval_ana] = fmincon(@(param)objfun_NDMwS_analytic(param,seed_location,pathology,time_stamps,C),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDMwS with the optimal parameters
    beta_ana = param_ana(1);
    x0_ana = seed_location*param_ana(2);
    alpha1_ana = param_ana(3);
    yana = NDMwS_analytic(x0_ana,time_stamps,C,beta_ana,alpha1_ana);
    
    %DEBUG
    %beta_ana = param_num(1);
    %x0_ana = seed_location*param_num(2);
    %yana = NDM_analytic(x0_ana,time_stamps,C,beta_ana);

    
 %% Display Numeric Results  
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(ynum(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDMwS minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['NUMERIC Beta = ' num2str(param_num(1)) ', x0 value = ' num2str(param_num(2)) ...
        ', alpha1 = ' num2str(param_num(3))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval_num)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(ynum,pathology,time_stamps)

%% Display Analytic Results   
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(yana(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDMwS minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['ANALYTIC Beta = ' num2str(param_ana(1)) ', x0 value = ' num2str(param_ana(2)) ...
         ', alpha1 = ' num2str(param_ana(3))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval_ana)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(yana,pathology,time_stamps)
clearvars y 

%% Compare values   
%table(ynum,yana)
figure
subplot(3,1,1)
plot(ynum(:,1),yana(:,1),'o')
hold on
line([0 3], [0 3],'Color','k')
%
subplot(3,1,2)
plot(ynum(:,2),yana(:,2),'o')
hold on
line([0 3], [0 3],'Color','k')
%
subplot(3,1,3)
plot(ynum(:,3),yana(:,3),'o')
hold on
line([0 3], [0 3],'Color','k')
%
suptitle('Comparison: Numeric vs Analytic') 

disp('rel error at t1') 
disp(    norm(ynum(:,1)-yana(:,1),2) /norm(ynum(:,1),2)      )

disp('rel error at t2') 
disp(    norm(ynum(:,2)-yana(:,2),2) /norm(ynum(:,2),2)      )

disp('rel error at t3') 
disp(    norm(ynum(:,3)-yana(:,3),2) /norm(ynum(:,3),2)      )


table(ynum,yana)
end

%% A.  [Numeric] NDMwC, fmincon minimizing MSE, no analytic counterpart 

% Set ndmwc_comp = 1 to run 
ndmwc_comp = 0;

if ndmwc_comp ==1

% Parameters to fit
    % param(1) = beta
    % param(2) = x0_value
    % param(3) = alpha0  cte growth/decay term independent of x
    
% Extra input to required by objective function
    % seed_location 
    % pathology
    % time_stamps
    % C
    
% Set initial guess for beta
    init_guess_params(1) = 1;

% Set initial guess for x0_value
    init_guess_params(2) = nansum(pathology(:,1));

% Set initial guess for alpha1 and alpha2
    init_guess_params(3) = .5;
    
% Set lower bounds (lb) for parameters
    lb = [0,0,-3];

% Set upper bounds (ub) for parameters
%    ub = [10,nansum(pathology(:,end)),5];
    ub = [3,3,3];

% Apply fmincon w/ numeric 
    [param_num, fval_num] = fmincon(@(param)objfun_NDMwC_numeric(param,seed_location,pathology,time_stamps,C),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDMwS with the optimal parameters
    beta_num = param_num(1);
    x0_num = seed_location*param_num(2);
    alpha0_num = param_num(3);
    ynum = NDMwC_numeric(x0_num,time_stamps,C,beta_num,alpha0_num);
      
 %% Display Numeric Results  
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(ynum(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDMwC minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['NUMERIC Beta = ' num2str(param_num(1)) ', x0 value = ' num2str(param_num(2)) ...
        ', alpha0 = ' num2str(param_num(3))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval_num)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(ynum,pathology,time_stamps)
clearvars y 
end

%% B.  [Numeric] NDMwSwC, fmincon minimizing MSE, no analytic counterpart 
%
% NDM with Source with Cte term (NDMwSwC)

% Set ndmwc_comp = 1 to run 
ndmwswc_comp = 0;

if ndmwswc_comp ==1

% Parameters to fit
    % param(1) = beta
    % param(2) = x0_value
    % param(3) = alpha0  cte growth/decay term independent of x
    % param(4) = alpha1 = linear growth/clearance term in x
    
% Extra input to required by objective function
    % seed_location 
    % pathology
    % time_stamps
    % C
    
% Set initial guess for beta
    init_guess_params(1) = 1;

% Set initial guess for x0_value
    init_guess_params(2) = nansum(pathology(:,1));

% Set initial guess for alpha1 and alpha2
    init_guess_params(3) = .5;
    init_guess_params(4) = .5;

% Set lower bounds (lb) for parameters
    lb = [0,0,-3,-3];

% Set upper bounds (ub) for parameters
%    ub = [10,nansum(pathology(:,end)),5];
    ub = [3,3,3,3];

% Apply fmincon w/ numeric 
    [param_num, fval_num] = fmincon(@(param)objfun_NDMwSwC_numeric(param,seed_location,pathology,time_stamps,C),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDMwS with the optimal parameters
    beta_num = param_num(1);
    x0_num = seed_location*param_num(2);
    alpha0_num = param_num(3);
    alpha1_num = param_num(4);
    ynum = NDMwSwC_numeric(x0_num,time_stamps,C,beta_num,alpha0_num,alpha1_num);
      
 %% Display Numeric Results  
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(ynum(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDMwSwC minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['NUMERIC Beta = ' num2str(param_num(1)) ', x0 value = ' num2str(param_num(2)) ...
        ', alpha0 = ' num2str(param_num(3))  ', alpha1 = ' num2str(param_num(4))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval_num)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(ynum,pathology,time_stamps)
clearvars y 
end

%% C.  [Numeric] NDMwMC, fmincon minimizing MSE, no analytic counterpart 
%
% NDM with Microglia Cte term (NDMwMC)

% Set ndmwmc_comp = 1 to run 
ndmwmc_comp = 0;

if ndmwmc_comp ==1

    
for j = 1:25
    
disp([' Results for Cell Type #' num2str(j)])
u = celltypedata(:,j); 
u = (u - min(u))./(max(u) - min(u));
    
% Parameters to fit
    % param(1) = beta
    % param(2) = x0_value
    % param(3) = alpha0  cte growth/decay term dependent of u
    
% Extra input to required by objective function
    % seed_location 
    % pathology
    % time_stamps
    % C
    
% Set initial guess for beta
    init_guess_params(1) = 1;

% Set initial guess for x0_value
    init_guess_params(2) = nansum(pathology(:,1));

% Set initial guess for alpha0 
    init_guess_params(3) = .5;

% Set lower bounds (lb) for parameters
    lb = [0,0,-3];

% Set upper bounds (ub) for parameters
%    ub = [10,nansum(pathology(:,end)),5];
    ub = [3,3,3];

% Apply fmincon w/ numeric 
    [param_num, fval_num] = fmincon(@(param)objfun_NDMwMC_numeric(param,seed_location,pathology,time_stamps,C,u),...
                           init_guess_params,[],[],[],[],lb,ub,[]);
 
% Solve NDMwS with the optimal parameters
    beta_num = param_num(1);
    x0_num = seed_location*param_num(2);
    alpha0_num = param_num(3);
    ynum = NDMwMC_numeric(x0_num,time_stamps,C,u,beta_num,alpha0_num);
      
 %% Display Numeric Results  
 % Save Rvalues in a matrix 
            for jj = 1:length(time_stamps)
                 Rvalues(jj,:) = corr(ynum(:,jj),pathology(:,jj), 'rows','complete');
            end

 % Display results
    disp('--------------------------------------------------')
    disp('NDMwMC minimizing quadratic error at all time stamps with fmincon');   
    disp(' ')
    disp(['NUMERIC Beta = ' num2str(param_num(1)) ', x0 value = ' num2str(param_num(2)) ...
        ', alpha0 = ' num2str(param_num(3))])
    disp(' ')
    disp(['R values at each time stamp'])
    disp(Rvalues)
    disp(' ')
    disp('Square error') 
    disp(fval_num)
    
% Plot prediction vs data using optimal parameters
plot_pred_vs_data(ynum,pathology,time_stamps)
clearvars y 
end
end

