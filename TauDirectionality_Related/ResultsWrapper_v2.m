%% 0. Loading data
rng(0); clear; clc; close all;
codedir = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis'; cd(codedir);
matdir = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
figdir = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Figures/TauDirectionality';
output_dir = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Results_Tables_TauDir';
load([matdir filesep 'Connectomes.mat'],'Connectomes');
load([matdir filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct')
studynames = fieldnames(mousedata_struct);
studynames(ismember(studynames,'IbaP301S')) = []; % Remove this study, too few time points
C = Connectomes.default; 

%% 1. Model-free analysis
%% 1.1 Connectome heatmap
savenclose = 0;
AMBCAHeatmap(C,savenclose,figdir);

%% 1.2 Graph metric analyses
studynames_plot = {'IbaStrInj'};
savenclose = 0;
whichplots = {'C_seed','All'};
tpt_flag = 'All';
for i = 1:length(studynames_plot)
    for j = 1:length(whichplots)
        DegreeEigenvectorSeedPlots(mousedata_struct,studynames_plot{i},...
            C,whichplots{j},matdir,savenclose,figdir);
    end
end
[R_vals, ttest_results] = CompareGraphMetricPlot(mousedata_struct,C,tpt_flag,...
    matdir,savenclose,figdir);

%% 1.3 Pathology brainframes, all timepoints
% datasets_bf = studynames; % note that xfac weighting may need to be 
% adjusted for different datasets; currently set up for IbaStrInj, which is
% featured in the manuscript
datasets_bf = {'IbaStrInj'};
tptsplot = 'All'; 
savenclose = 0;
for i = 1:length(datasets_bf)
    BrainframePathologyPlot(mousedata_struct,datasets_bf{i},tptsplot,matdir,...
        savenclose,figdir)
end

%% 1.4 Connectivity to/from seed brainframes
dataset_bf = 'IbaStrInj';
tptsplot = 3;
seedconntypes = {'In','Out'};
savenclose = 0;
for i = 1:length(seedconntypes)
    BrainframeSeedConnectivityPlot(mousedata_struct,dataset_bf,tptsplot,...
        C,seedconntypes{i},matdir,savenclose,figdir);
end

%% 2. NexIS:global w/directionality modeling
% fit longitudinally alpha/beta/s (if fit_s)
% fix gamma and alpha for per-timepoint, use LinR

%% 2.1 Longitudinal models
% Input parameters
saveoutputs = 1;
filename_out = 'outputs_all';
outputs_all = struct;
modelnames = {'fit_s','ret','ant','nd'};
use_dataspace = 1;
w_dir = 1;
volcorrect = 1;
bootstrapping = 0;
exclseed_outputs = 0;
param_init = [NaN,0.5,1,0.5];

% Run NexIS_global
for i = 1:length(studynames)
    tablename = [filename_out '_' studynames{i}];
    sumtable = [];
    ub = [Inf,Inf,Inf,1]; ubs = repmat(ub,4,1); ubs(:,end) = [1,1,0,0.5].';
    lb = [0,0,0,0]; lbs = repmat(lb,4,1); lbs(:,end) = [0,1,0,0.5].';
    for j = 1:length(modelnames)
        fprintf('Study %d of %d, Model %s\n',i,length(studynames),modelnames{j})
        outputs = NexIS_global('study',studynames{i},...
                                'w_dir',w_dir,...
                                'volcorrect',volcorrect,...
                                'param_init',param_init,...
                                'ub',ubs(j,:),...
                                'lb',lbs(j,:),...
                                'use_dataspace',use_dataspace,...
                                'costfun','linr',...
                                'bootstrapping',bootstrapping,...
                                'exclseed_outputs',exclseed_outputs);
        outputs_all.(studynames{i}).(modelnames{j}) = outputs;
        sumtable_i = Output2Table(outputs,0,'null','null');
        sumtable_i.Properties.RowNames{1} = modelnames{j};
        sumtable = [sumtable; sumtable_i];
    end
    if saveoutputs
        writetable(sumtable,[output_dir filesep tablename '.csv'],'WriteRowNames',true)
    end
end
if saveoutputs
    save([output_dir filesep filename_out '.mat'],'outputs_all');
end

%% 2.2 Figures per 2.1
preload = 1;
filename_out = 'outputs_all';
if preload
    load([output_dir filesep filename_out '.mat'],'outputs_all');
end
savenclose = 0;
pertpt = 0;
tpt_plot = 3; % vs. last time point
datset_bf = 'IbaStrInj';

RvstPlots(outputs_all,tpt_plot,1,matdir,savenclose,figdir);
[R,s,tstatstruct] = CompareDirPlots_deltaR_s(outputs_all,pertpt,savenclose,figdir);
% save([output_dir filesep 'CompareDirLong.mat'],'R','s','tstatstruct');
BrainframeModelPredPlot(outputs_all,datset_bf,matdir,savenclose,figdir);

%% 2.3 Per-timepoint models, Lin R cost function, fix gamma and alpha
output_dir = '/Users/justintorok/Documents/MATLAB/Nexis_Project'; % DON'T USE 
% Input parameters
saveoutputs = 1;
outputs_all_tpt = struct;
filename_out = 'outputs_all_tpt_fixgammaalpha';
modelnames = {'fit_s','ret','ant','nd'};
costfun = 'linr';
use_dataspace = 1;
w_dir = 1;
volcorrect = 1;
bootstrapping = 0;
exclseed_outputs = 0;
studynames = {'IbaStrInj','Hurtado'};

% Run NexIS_global
for i = 1:length(studynames)
    tablename = [filename_out '_' studynames{i}];
    sumtable = [];
    for j = 1:length(modelnames)
        fprintf('Study %d of %d, Model %s\n',i,length(studynames),modelnames{j})
        params_opt = outputs_all.(studynames{i}).(modelnames{j}).nexis_global.Full.param_fit;
        gammaval = params_opt(1); alphaval = params_opt(2); % Fix gamma/alpha to longitudinal vals
        ub = [gammaval,alphaval,Inf,1]; ubs = repmat(ub,4,1); ubs(:,end) = [1,1,0,0.5].';
        lb = [gammaval,alphaval,0,0]; lbs = repmat(lb,4,1); lbs(:,end) = [0,1,0,0.5].';
        excl_tpts = [[2,3];[1,3];[1,2]];
        for k = 1:size(excl_tpts,1)
            fprintf('Timepoint %d of %d\n',k,size(excl_tpts,1))
            excl_tpt = excl_tpts(k,:);
            tpt_str = ['t_' num2str(setdiff(1:3,excl_tpt))];
            outputs = NexIS_global('study',studynames{i},...
                                    'w_dir',w_dir,...
                                    'volcorrect',volcorrect,...
                                    'param_init',param_init,...
                                    'ub',ubs(j,:),...
                                    'lb',lbs(j,:),...
                                    'use_dataspace',use_dataspace,...
                                    'bootstrapping',bootstrapping,...
                                    'costfun',costfun,...
                                    'excltpts_costfun',excl_tpt,...
                                    'exclseed_outputs',exclseed_outputs);
            outputs_all_tpt.(studynames{i}).(modelnames{j}).(tpt_str) = outputs;
            sumtable_i = Output2Table(outputs,0,'null','null');
            sumtable_i.Properties.RowNames{1} = [modelnames{j} ', ' tpt_str];
            sumtable_i.Properties.VariableNames{16} = 'R';
            sumtable = [sumtable; sumtable_i];
        end
    end
    if saveoutputs
        writetable(sumtable,[output_dir filesep tablename '.csv'],'WriteRowNames',true)
    end
end
if saveoutputs
    save([output_dir filesep filename_out '.mat'],'outputs_all_tpt');
end

%% 2.4 Figures per 2.3
filename_out = 'outputs_all_tpt_fixgammaalpha';
preload = 1;
if preload
    load([output_dir filesep filename_out '.mat'],'outputs_all_tpt');
end
savenclose = 0;
pertpt = 1;

[R_fix,s_fix,tstatstruct_fix] = CompareDirPlots_deltaR_s(outputs_all_tpt,pertpt,savenclose,figdir);
% save([output_dir filesep 'CompareDirPerTpt.mat'],'R_fix','s_fix','tstatstruct_fix');
plottypes = {'alpha_s','beta_s','beta_alpha'};
for i = 1:length(plottypes)
    [amat,bmat,smat] = salphabetaPlot(outputs_all_tpt,plottypes{i},savenclose,figdir);
end

for i = 1:2
    PerTimepointPlot_sbeta(outputs_all_tpt,i-1,savenclose,figdir);
    PerTimepointRegressionPlot_sbeta(outputs_all_tpt,i-1,savenclose,figdir);
    % CorrComparePlot_Combined(outputs_all,outputs_all_tpt,i-1,savenclose,figdir);
end

%% 2.5.1 Longitudinal models, bootstrapping
% Input parameters
saveoutputs = 1;
filename_out = 'outputs_all_bs';
outputs_all = struct;
modelnames = {'fit_s'};
use_dataspace = 1;
w_dir = 1;
volcorrect = 1;
bootstrapping = 1;
niters = 100;
exclseed_outputs = 0;
param_init = [NaN,0.5,1,0.5];

% Run NexIS_global
for i = 1:length(studynames)
    tablename = [filename_out '_' studynames{i}];
    sumtable = [];
    ub = [Inf,Inf,Inf,1]; ubs = repmat(ub,4,1); ubs(:,end) = [1,1,0,0.5].';
    lb = [0,0,0,0]; lbs = repmat(lb,4,1); lbs(:,end) = [0,1,0,0.5].';
    for j = 1:length(modelnames)
        fprintf('Study %d of %d, Model %s\n',i,length(studynames),modelnames{j})
        outputs = NexIS_global('study',studynames{i},...
                                'w_dir',w_dir,...
                                'volcorrect',volcorrect,...
                                'param_init',param_init,...
                                'ub',ubs(j,:),...
                                'lb',lbs(j,:),...
                                'use_dataspace',use_dataspace,...
                                'costfun','linr',...
                                'bootstrapping',bootstrapping,...
                                'niters',niters,...
                                'exclseed_outputs',exclseed_outputs);
        outputs_all.(studynames{i}).(modelnames{j}) = outputs;
        sumtable_i = Output2Table(outputs,0,'null','null');
        sumtable_i.Properties.RowNames{1} = modelnames{j};
        sumtable = [sumtable; sumtable_i];
    end
    if saveoutputs
        writetable(sumtable,[output_dir filesep tablename '.csv'],'WriteRowNames',true)
    end
end
if saveoutputs
    save([output_dir filesep filename_out '.mat'],'outputs_all');
end

%% 2.6 Figure per 2.3
preload = 1;
filename_out_bs = 'outputs_all_bs';
filename_out = 'outputs_all';
if preload
    load([output_dir filesep filename_out_bs '.mat'],'outputs_all');
    outputs_all_bs = outputs_all;
    load([output_dir filesep filename_out '.mat'],'outputs_all');
end
paramnames = {'alpha','beta','s','R'};
savenclose = 0;
for i = 1:length(paramnames)
    BootstrappingPlotter_NexIS_1param(outputs_all_bs,outputs_all,paramnames{i},...
        'fit_s',savenclose,figdir);
end
