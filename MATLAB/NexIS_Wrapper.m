% Demo and debugging script for all major NexIS functions. Make sure you
% are running this code in the top-level 'Nexis' folder of the repository!

%% 1. NexIS:global
%% 1.1 No bootstrapping of parameters, preloaded tau studies
% Running NexIS
rng(0); clear; clc;
studylist = {'IbaHippInj','Hurtado'}; % Cell array of test datasets
wdir = [0,1]; % Toggle directionality fitting (s) on and off
volcorrect = 1; % ***Always keep set to 1***
usedataspace = [0,1]; % Toggle data space fitting on and off
param_init = [NaN,0,1,0.5]; % Initial fmincon parameter guesses; {gamma, alpha, beta, s}
ub = [Inf,Inf,Inf,1]; % Upper bounds for fmincon
lb = zeros(1,4); % Lower bounds for fmincon
excltpts_costfun = {[],1}; % Exclude selected time points from cost function
bootstrapping_glob = 0; % Flag for bootstrapping
niters_glob = 3; % Number of bootstrapped iterations (does nothing if bootstrap flag is 0)

% Table output parameters
writetofile = 1; % Create .csv from MATLAB table
filename_out = 'NexIS_Wrapper_Global_1-1-1_NoBootstrap'; % Name of output file
filepath_out = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Save path

% Run model and create output tables for each dataset, if writetofile = 1
outputs_all = struct;
numsims = length(studylist)*length(wdir)*length(usedataspace)*length(excltpts_costfun);
for i = 1:length(studylist)
    study_i = studylist{i};
    tablename = [filename_out '_' study_i]; % Create one output table per study
    sumtable_i = [];
    for j = 1:length(wdir)
        for m = 1:length(excltpts_costfun)
            for k = 1:length(usedataspace)
                tablerowname = ['global, ' study_i ', '];
                wdir_j = wdir(j);
                if wdir_j
                    tablerowname = [tablerowname 'with dir., '];
                else
                    tablerowname = [tablerowname 'no dir., '];
                end
                if isempty(excltpts_costfun{m})
                    tablerowname = [tablerowname 'all tpts, '];
                else
                    tablerowname = [tablerowname 'excl tpts, '];
                end
                useds_k = usedataspace(k);
                simno = length(wdir)*length(usedataspace)*length(excltpts_costfun)*(i-1)...
                    + length(excltpts_costfun)*(m-1) + length(excltpts_costfun)*length(wdir)*(j-1) + k;
                fprintf('NexIS Wrapper Test %d/%d\n',simno,numsims)
                if useds_k
                    tablerowname = [tablerowname 'data space'];
                else
                    tablerowname = [tablerowname 'CCF space'];
                end
                outputs_ng_ijkm = NexIS_global('study',study_i,...
                                              'w_dir',wdir_j,...
                                              'use_dataspace',useds_k,...
                                              'bootstrapping',bootstrapping_glob,...
                                              'niters',niters_glob,...
                                              'volcorrect',volcorrect,...
                                              'param_init',param_init,...
                                              'ub',ub,...
                                              'lb',lb,...
                                              'excltpts_costfun',excltpts_costfun{m});
                fieldname_ijkm = [study_i '_wdir_' num2str(wdir_j) '_useds_' num2str(useds_k) '_excl_tpts_' num2str(isempty(excltpts_costfun{m}))];
                outputs_all.(fieldname_ijkm) = outputs_ng_ijkm;
                sumtable_ijkm = Output2Table(outputs_ng_ijkm,0,'null','null'); % create table row
                sumtable_ijkm.Properties.RowNames{1} = tablerowname; % label table row
                sumtable_i = [sumtable_i; sumtable_ijkm]; % add row to table
            end
        end
    end
    if writetofile
        writetable(sumtable_i,[filepath_out filesep tablename '.csv'],'WriteRowNames',true)
    end
end

if writetofile
    save([filepath_out filesep filename_out '.mat'],'outputs_all');
end 

%% 1.2 With bootstrapping of parameters, preloaded tau studies
% Running NexIS
rng(0); clc;
studylist = {'IbaHippInj','Hurtado'}; % Cell array of test datasets
wdir = [0,1]; % Toggle directionality fitting (s) on and off
volcorrect = 1; % ***Always keep set to 1***
usedataspace = [0,1]; % Toggle data space fitting on and off
param_init = [NaN,0,1,0.5]; % Initial fmincon parameter guesses; {gamma, alpha, beta, s}
ub = [Inf,Inf,Inf,1]; % Upper bounds for fmincon
lb = zeros(1,4); % Lower bounds for fmincon
excltpts_costfun = {[],1}; % Exclude selected time points from cost function
bootstrapping_glob = 1; % Flag for bootstrapping
niters_glob = 3; % Number of bootstrapped iterations (does nothing if bootstrap flag is 0)

% Table output parameters
writetofile = 1; % Create .csv from MATLAB table
filename_out = 'NexIS_Wrapper_Global_1-1-1_WithBootstrap'; % Name of output file
filepath_out = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Save path

% Run model and create output tables for each dataset, if writetofile = 1
outputs_all = struct;
numsims = length(studylist)*length(wdir)*length(usedataspace)*length(excltpts_costfun);
for i = 1:length(studylist)
    study_i = studylist{i};
    tablename = [filename_out '_' study_i]; % Create one output table per study
    sumtable_i = [];
    for j = 1:length(wdir)
        for m = 1:length(excltpts_costfun)
            for k = 1:length(usedataspace)
                tablerowname = ['global, ' study_i ', '];
                wdir_j = wdir(j);
                if wdir_j
                    tablerowname = [tablerowname 'with dir., '];
                else
                    tablerowname = [tablerowname 'no dir., '];
                end
                if isempty(excltpts_costfun{m})
                    tablerowname = [tablerowname 'all tpts, '];
                else
                    tablerowname = [tablerowname 'excl tpts, '];
                end
                useds_k = usedataspace(k);
                simno = length(wdir)*length(usedataspace)*length(excltpts_costfun)*(i-1)...
                    + length(excltpts_costfun)*(m-1) + length(excltpts_costfun)*length(wdir)*(j-1) + k;
                fprintf('NexIS Wrapper Test %d/%d\n',simno,numsims)
                if useds_k
                    tablerowname = [tablerowname 'data space'];
                else
                    tablerowname = [tablerowname 'CCF space'];
                end
                outputs_ng_ijkm = NexIS_global('study',study_i,...
                                              'w_dir',wdir_j,...
                                              'use_dataspace',useds_k,...
                                              'bootstrapping',bootstrapping_glob,...
                                              'niters',niters_glob,...
                                              'volcorrect',volcorrect,...
                                              'param_init',param_init,...
                                              'ub',ub,...
                                              'lb',lb,...
                                              'excltpts_costfun',excltpts_costfun{m});
                fieldname_ijkm = [study_i '_wdir_' num2str(wdir_j) '_useds_' num2str(useds_k) '_excl_tpts_' num2str(isempty(excltpts_costfun{m}))];
                outputs_all.(fieldname_ijkm) = outputs_ng_ijkm;
                sumtable_ijkm = Output2Table(outputs_ng_ijkm,0,'null','null'); % create table row
                sumtable_ijkm.Properties.RowNames{1} = tablerowname; % label table row
                sumtable_i = [sumtable_i; sumtable_ijkm]; % add row to table
            end
        end
    end
    if writetofile
        writetable(sumtable_i,[filepath_out filesep tablename '.csv'],'WriteRowNames',true)
    end
end

if writetofile
    save([filepath_out filesep filename_out '.mat'],'outputs_all');
end 

%% 1.3 User-specified connectome & pathology

%% 1.4 A-syn stuff
% Deal with later

%% 1.5 Relevant plotting of outputs
% Deal with later

%% 2. NexIS:SV
%% 2.1 No bootstrapping of parameters, preloaded tau studies

%% 2.2 With bootstrapping of parameters, preloaded tau studies

%% 2.3 User-specified connectome & pathology

%% 2.4 A-syn stuff
% Deal with later

%% 2.5 Relevant plotting of outputs
% Deal with later

%%
%     outputs_ct_ds_nobs = NexIS_SV('study',studylist{i},...
%                                   'bootstrapping',0,...
%                                   'niters',5,...
%                                   'w_dir',1,...
%                                   'volcorrect',1,...
%                                   'param_init',[NaN,0,1,0.5],...
%                                   'ub',[Inf,Inf,Inf,1],...
%                                   'lb',zeros(1,4),...
%                                   'use_dataspace',1,...
%                                   'datatype_nexis_sv','gene',...
%                                   'datalist_nexis_sv',{'Pvalb','Vip','Sst'});

%%
% C = Connectomes.default;
% C = C/max(C(:));
% seed = seed426.DS4;
% x0 = seed * 0.8;
% time_stamps = [4,8,12];
% y = eNDM_general_dir(x0,time_stamps,C,zeros(426,1),0.2,3.5,0.65,0,0,0,'analytic',0);
% corr(y(:), data426.DS4(:),'rows','complete')^2