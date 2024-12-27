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
filename_out = 'NexIS_Wrapper_Global_1-1_NoBootstrap'; % Name of output file
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
                fieldname_ijkm = [study_i '_wdir_' num2str(wdir_j) '_useds_' num2str(useds_k) '_excl_tpts_' num2str(~isempty(excltpts_costfun{m}))];
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
wdir = 1; % Toggle directionality fitting (s) on and off
volcorrect = 0; % ***Always keep set to 1***
usedataspace = 0; % Toggle data space fitting on and off
param_init = [NaN,0,1,0.5]; % Initial fmincon parameter guesses; {gamma, alpha, beta, s}
ub = [Inf,Inf,Inf,1]; % Upper bounds for fmincon
lb = zeros(1,4); % Lower bounds for fmincon
excltpts_costfun = {[]}; % Exclude selected time points from cost function
bootstrapping_glob = 1; % Flag for bootstrapping
niters_glob = 3; % Number of bootstrapped iterations (does nothing if bootstrap flag is 0)

% Table output parameters
writetofile = 1; % Create .csv from MATLAB table
filename_out = 'NexIS_Wrapper_Global_1-2_WithBootstrap_ccf_test'; % Name of output file
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
                fieldname_ijkm = [study_i '_wdir_' num2str(wdir_j) '_useds_' num2str(useds_k) '_excl_tpts_' num2str(~isempty(excltpts_costfun{m}))];
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

%% 1.3 Preloaded a-syn studies (simple test)
% Running NexIS
rng(0); clear; clc;
studylist = {'asyn_mouse','asyn_human','Henderson','PFF','GCI'}; % Cell array of test datasets
wdir = 1; % Toggle directionality fitting (s) on and off
volcorrect = [1,1,0,0,0]; % ***Can only be set to 1 for asyn_mouse/human studies***
usedataspace = 0; % ***Required to be 0 for a-syn, all data in native connectome space***
param_init = [NaN,0,1,0.5]; % Initial fmincon parameter guesses; {gamma, alpha, beta, s}
ub = [Inf,Inf,Inf,1]; % Upper bounds for fmincon
lb = zeros(1,4); % Lower bounds for fmincon
excltpts_costfun = []; % Exclude selected time points from cost function
bootstrapping_glob = 0; % Flag for bootstrapping

% Table output parameters
writetofile = 1; % Create .csv from MATLAB table
filename_out = 'NexIS_Wrapper_Global_1-3_NoBootstrap_a-syn'; % Name of output file
filepath_out = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Save path

% Run model and create output tables for each dataset, if writetofile = 1
outputs_all = struct;
numsims = length(studylist);
for i = 1:length(studylist)
    study_i = studylist{i};
    tablename = [filename_out '_' study_i]; % Create one output table per study
    tablerowname = ['global, ' study_i ', '];
    if wdir
        tablerowname = [tablerowname 'with dir., '];
    else
        tablerowname = [tablerowname 'no dir., ']; %#ok<UNRCH>
    end
    if isempty(excltpts_costfun)
        tablerowname = [tablerowname 'all tpts, '];
    else
        tablerowname = [tablerowname 'excl tpts, '];
    end
    simno = i;
    fprintf('NexIS Wrapper Test %d/%d\n',simno,numsims)
    outputs_ng_i = NexIS_global('study',study_i,...
                                  'w_dir',wdir,...
                                  'use_dataspace',usedataspace,...
                                  'bootstrapping',bootstrapping_glob,...
                                  'niters',niters_glob,...
                                  'volcorrect',volcorrect(i),...
                                  'param_init',param_init,...
                                  'ub',ub,...
                                  'lb',lb,...
                                  'excltpts_costfun',excltpts_costfun);
    fieldname_i = [study_i '_wdir_' num2str(wdir) '_excl_tpts_' num2str(~isempty(excltpts_costfun))];
    outputs_all.(fieldname_i) = outputs_ng_i;
    sumtable_i = Output2Table(outputs_ng_i,0,'null','null'); % create table row
    sumtable_i.Properties.RowNames{1} = tablerowname; % label table row
    if writetofile
        writetable(sumtable_i,[filepath_out filesep tablename '.csv'],'WriteRowNames',true)
    end
end

if writetofile
    save([filepath_out filesep filename_out '.mat'],'outputs_all');
end

%% 1.4 No bootstrapping, user-specified connectome & pathology
% Running NexIS
rng(0); clear; clc;
matdir = '~/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
studylist = {'IbaHippInj','Hurtado'}; % Cell array of test dataset names
wdir = 1; % Toggle directionality fitting (s) on and off
volcorrect = 0; % ***Can only be set to 1 for preloaded tau studies and asyn_mouse/human***
usedataspace = 0; % ***Leave off when user-specified dataset***
param_init = [NaN,0,1,0.5]; % Initial fmincon parameter guesses; {gamma, alpha, beta, s}
ub = [Inf,Inf,Inf,1]; % Upper bounds for fmincon
lb = zeros(1,4); % Lower bounds for fmincon
excltpts_costfun = []; % Exclude selected time points from cost function
bootstrapping_glob = 0; % Flag for bootstrapping

% Table output parameters
writetofile = 1; % Create .csv from MATLAB table
filename_out = 'NexIS_Wrapper_Global_1-4_NoBootstrap_user-specified'; % Name of output file
filepath_out = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Save path

% Load relevant data files for test (theoretically, could use anything)
load([matdir filesep 'Mouse_Tauopathy_Data_HigherQ_CCF.mat'],'mousedata_struct_ccf'); % CCF versions of data
load([matdir filesep 'Connectomes.mat'],'Connectomes');

% Connectome definition and min-max normalization
C = Connectomes.default;
cmax = max(max(C));
cmin = min(min(C));
C = (C - cmin)./(cmax-cmin);

% Run model and create output tables for each dataset, if writetofile = 1
outputs_all = struct;
numsims = length(studylist);
for i = 1:length(studylist)
    study_i = studylist{i};
    data_i = mousedata_struct_ccf.(study_i).data; % Extract data
    ts_i = mousedata_struct_ccf.(study_i).time_stamps; % Extract time points
    seed_i = mousedata_struct_ccf.(study_i).seed; % Extract seed
    tablename = [filename_out '_' study_i]; % Create one output table per study
    tablerowname = ['global, ' study_i ' CCF, '];
    if wdir
        tablerowname = [tablerowname 'with dir., '];
    else
        tablerowname = [tablerowname 'no dir., ']; %#ok<UNRCH>
    end
    if isempty(excltpts_costfun)
        tablerowname = [tablerowname 'all tpts'];
    else
        tablerowname = [tablerowname 'excl tpts'];
    end
    simno = i;
    fprintf('NexIS Wrapper Test %d/%d\n',simno,numsims)
    outputs_ng_i = NexIS_global('study','User_specified',...
                                  'C',C,...
                                  'data',data_i,...
                                  'seed',seed_i,...
                                  'tpts',ts_i,...
                                  'w_dir',wdir,...
                                  'use_dataspace',usedataspace,...
                                  'bootstrapping',bootstrapping_glob,...
                                  'volcorrect',volcorrect,...
                                  'param_init',param_init,...
                                  'ub',ub,...
                                  'lb',lb,...
                                  'excltpts_costfun',excltpts_costfun);
    fieldname_i = [study_i '_wdir_' num2str(wdir) '_excl_tpts_' num2str(~isempty(excltpts_costfun))];
    outputs_all.(fieldname_i) = outputs_ng_i;
    sumtable_i = Output2Table(outputs_ng_i,0,'null','null'); % create table row
    sumtable_i.Properties.RowNames{1} = tablerowname; % label table row
    if writetofile
        writetable(sumtable_i,[filepath_out filesep tablename '.csv'],'WriteRowNames',true)
    end
end

if writetofile
    save([filepath_out filesep filename_out '.mat'],'outputs_all');
end

%% 1.5 With bootstrapping, user-specified connectome & pathology
% Running NexIS
rng(0); clear; clc;
matdir = '~/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
studylist = {'IbaHippInj','Hurtado'}; % Cell array of test dataset names
wdir = 1; % Toggle directionality fitting (s) on and off
volcorrect = 0; % ***Can only be set to 1 for preloaded tau studies and asyn_mouse/human***
usedataspace = 0; % ***Leave off when user-specified dataset***
param_init = [NaN,0,1,0.5]; % Initial fmincon parameter guesses; {gamma, alpha, beta, s}
ub = [Inf,Inf,Inf,1]; % Upper bounds for fmincon
lb = zeros(1,4); % Lower bounds for fmincon
excltpts_costfun = []; % Exclude selected time points from cost function
bootstrapping_glob = 1; % Flag for bootstrapping
niters_glob = 3; % Number of bootstrapped iterations (does nothing if bootstrap flag is 0)

% Table output parameters
writetofile = 1; % Create .csv from MATLAB table
filename_out = 'NexIS_Wrapper_Global_1-5_WithBootstrap_user-specified'; % Name of output file
filepath_out = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Save path

% Load relevant data files for test (theoretically, could use anything)
load([matdir filesep 'Mouse_Tauopathy_Data_HigherQ_CCF.mat'],'mousedata_struct_ccf'); % CCF versions of data
load([matdir filesep 'Connectomes.mat'],'Connectomes');

% Connectome definition and min-max normalization
C = Connectomes.default;
cmax = max(max(C));
cmin = min(min(C));
C = (C - cmin)./(cmax-cmin);

% Run model and create output tables for each dataset, if writetofile = 1
outputs_all = struct;
numsims = length(studylist);
for i = 1:length(studylist)
    study_i = studylist{i};
    data_i = mousedata_struct_ccf.(study_i).data; % Extract data
    ts_i = mousedata_struct_ccf.(study_i).time_stamps; % Extract time points
    seed_i = mousedata_struct_ccf.(study_i).seed; % Extract seed
    tablename = [filename_out '_' study_i]; % Create one output table per study
    tablerowname = ['global, ' study_i ' CCF, '];
    if wdir
        tablerowname = [tablerowname 'with dir., '];
    else
        tablerowname = [tablerowname 'no dir., ']; %#ok<UNRCH>
    end
    if isempty(excltpts_costfun)
        tablerowname = [tablerowname 'all tpts'];
    else
        tablerowname = [tablerowname 'excl tpts'];
    end
    simno = i;
    fprintf('NexIS Wrapper Test %d/%d\n',simno,numsims)
    outputs_ng_i = NexIS_global('study','User_specified',...
                                  'C',C,...
                                  'data',data_i,...
                                  'seed',seed_i,...
                                  'tpts',ts_i,...
                                  'w_dir',wdir,...
                                  'use_dataspace',usedataspace,...
                                  'bootstrapping',bootstrapping_glob,...
                                  'niters',niters_glob,...
                                  'volcorrect',volcorrect,...
                                  'param_init',param_init,...
                                  'ub',ub,...
                                  'lb',lb,...
                                  'excltpts_costfun',excltpts_costfun);
    fieldname_i = [study_i '_wdir_' num2str(wdir) '_excl_tpts_' num2str(~isempty(excltpts_costfun))];
    outputs_all.(fieldname_i) = outputs_ng_i;
    sumtable_i = Output2Table(outputs_ng_i,0,'null','null'); % create table row
    sumtable_i.Properties.RowNames{1} = tablerowname; % label table row
    if writetofile
        writetable(sumtable_i,[filepath_out filesep tablename '.csv'],'WriteRowNames',true)
    end
end

if writetofile
    save([filepath_out filesep filename_out '.mat'],'outputs_all');
end

%% 1.6 Relevant plotting of outputs
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