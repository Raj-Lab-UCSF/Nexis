%% Regression w/ MISS + Hurtado 

%%  Setup   

% Initialization
    clear all; close all; clc; 
    
% Load datasets of interest
   load celltypedata_mrx3_540.mat
   load vulnerable_genes_expression.mat
   load Hurtado426.mat

%   load MouseDataForPedro.mat
   load Hurtado426.mat

% Read time stamps, pathology measurements, and the seed_location 
   
% Hurtado   
   time_stamps = [2,4,6,8];
   pathology = data426.Hurtado  ;

%% Min-max normalization    
   
% Normalize every column of cell type matrix T 
    for j = 1:size(celltypedata,2)    
        u = celltypedata(:,j); 
        u = (u - min(u))./(max(u) - min(u));
        T(:,j) = u;
    end

% Normalize every column of pathology matrix   
for j = 1:size(pathology,2)
        u = pathology(:,j);
        u = (u - min(u))./(max(u) - min(u));
        pathology(:,j) = u;
end

% Normalize every column of vulnerability expression genes matrix G
for j = 1:size(vuln_expr_reg,2)
     u = vuln_expr_reg(:,j);
     u = (u - min(u))./(max(u) - min(u));
     G(:,j) = u;
end

%% Stepwise LM 

% For vulnerability genes in G matrix
for j = 1:length(time_stamps)
    
    disp('-----------------------------------------')
    disp(['Stepwise regression Vuln. genes vs pathology at time  ' num2str(time_stamps(j))])
    disp(' ')
    mdl = stepwiselm(G,pathology(:,j),'PEnter',0.0001)
    disp('-----------------------------------------')
    disp(' ')
    disp(' ')

    % SAVE PREDICTION VECTORS FOR CHRIS
    % Predicted value by regression at every time stamp (column)
    % Pls do glassbrain of the desired column
    pred_with_genes(:,j) = mdl.Fitted; 
    
end


% For cell types in T matrix
for j = 1:length(time_stamps)
    
    disp('-----------------------------------------')
    disp(['Stepwise regression cell types vs pathology at time  ' num2str(time_stamps(j))])
    disp(' ')
    mdl = stepwiselm(T,pathology(:,j),'PEnter',0.00001)
    disp('-----------------------------------------')
    disp(' ')
    disp(' ')

    % SAVE PREDICTION VECTORS FOR CHRIS
    % Predicted value by regression at every time stamp (column)
    % Pls do glassbrain of the desired column
    pred_with_ctypes(:,j) = mdl.Fitted; 
    
end

%% Stepwise using only {'App'} and {'Mapt'}

sG = G(:,[4,13]);

% For vulnerability genes in G matrix
for j = 1:length(time_stamps)
    
    disp('-----------------------------------------')
    disp(['Stepwise regression Vuln. genes vs pathology at time  ' num2str(time_stamps(j))])
    disp(' ')
    mdl = stepwiselm(sG,pathology(:,j),'PEnter',0.01)
    disp('-----------------------------------------')
    disp(' ')
    disp(' ')

    % SAVE PREDICTION VECTORS FOR CHRIS
    % Predicted value by regression at every time stamp (column)
    % Pls do glassbrain of the desired column
    pred_only_app_mapt(:,j) = mdl.Fitted; 
end
    
%% Lasso (gave bad results...)

% Lasso for genes vs pathology at time stamp 4

    % [B_genes,FitInfo_genes] = lasso(G,pathology(:,4),'CV',10)
    % figure
    % lassoPlot(B_genes,FitInfo_genes,'PlotType','CV');
    % legend('show') % Show legend


% Lasso for cell types vs pathology at time stamp 4

    % [B_ctypes,FitInfo_ctypes] = lasso(T,pathology(:,4),'CV',10)
    % figure
    % lassoPlot(B_genes,FitInfo_genes,'PlotType','CV');
    % legend('show') % Show legend

