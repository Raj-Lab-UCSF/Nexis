% Demo and debugging script for all major NexIS functions. Make sure you
% are running this code in the top-level 'Nexis' folder of the repository!
%
%
%
%% 1. NexIS:global
%
%
%
%% 1.1 No bootstrapping of parameters, preloaded tau studies
%
%
%
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
                fprintf('NexIS Wrapper Test 1.1, %d/%d\n',simno,numsims)
                if useds_k
                    tablerowname = [tablerowname 'data space'];
                else
                    tablerowname = [tablerowname 'CCF space'];
                end
                fprintf('Simulation: %s\n',tablerowname)
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
%
%
%
%% 1.2 With bootstrapping of parameters, preloaded tau studies
%
%
%
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
                fprintf('NexIS Wrapper Test 1.2, %d/%d\n',simno,numsims)
                if useds_k
                    tablerowname = [tablerowname 'data space'];
                else
                    tablerowname = [tablerowname 'CCF space'];
                end
                fprintf('Simulation: %s\n',tablerowname)
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
%
%
%
%% 1.3 Preloaded a-syn studies (simple test)
%
%
%
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
    fprintf('NexIS Wrapper Test 1.3, %d/%d\n',simno,numsims)
    fprintf('Simulation: %s\n',tablerowname)
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
%
%
%
%% 1.4 No bootstrapping, user-specified connectome & pathology
%
%
%
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
    fprintf('NexIS Wrapper Test 1.4, %d/%d\n',simno,numsims)
    fprintf('Simulation: %s\n',tablerowname)
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
%
%
%
%% 1.5 With bootstrapping, user-specified connectome & pathology
%
%
%
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
    fprintf('NexIS Wrapper Test 1.5, %d/%d\n',simno,numsims)
    fprintf('Simulation: %s\n',tablerowname)
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
%
%
%
%% 2. NexIS:SV
%
%
%
%% 2.1 No bootstrapping of parameters, preloaded tau studies
%
%
%
% Loading previously run NexIS_global struct. Not required to do this to run
% NexIS:SV, but recommended for efficiency reasons, particularly if running
% through multiple factors (i.e., genes, cell types) for the same pathology
% dataset
rng(0); clear; clc;
filename_out = 'NexIS_Wrapper_SV_2-1_NoBootstrap'; % Name of output file
filename_in_glob = 'NexIS_Wrapper_Global_1-1_NoBootstrap'; % NexIS:global input file
filepath_in = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Load path, sims
filepath_out = filepath_in; % Save path
writetofile = 1;

% Pull names of (single) gene & cell types 
nexsvnames = cell(2,1);
nexsvnames{1} = 'Yao';
nexsvnames{2} = 'Oligo'; % Single SV factor for debugging purposes

% Note: lines of code below were used to test to make sure all kinds of
% 'datatype_nexis_sv' inputs worked, and this was successful even though
% the outputs were not all run and saved (12/31/24)
%
% fp_in_sv = '~/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
% load([fp_in_sv filesep 'CellTypeMaps.mat'],'CellTypeMaps');
% load([fp_in_sv filesep 'GeneExpressionMaps.mat'],'GeneExpressionMaps');
% nexsvnames = cell(2,6);
% nexsvnames(1,1:size(nexsvnames,2)) = {'gene','Yao','Tasic','Zeisel',...
%     'Zhuang_Class','Zhuang_Subclass'};
% for j = 1:length(nexsvnames)
%     if strcmp(nexsvnames{1,j},'gene')
%         allnames = GeneExpressionMaps.All.gene_names;
%     else
%         allnames = CellTypeMaps.(nexsvnames{1,j}).classkey;
%     end
%     testind_j = datasample(1:length(allnames),1);
%     nexsvnames{2,j} = allnames{testind_j};
% end

% Running NexIS:SV (Exhaustive test of no bootstrap)
nexglob_outputs = load([filepath_in filesep filename_in_glob '.mat'],'outputs_all');
mdlinstances = fieldnames(nexglob_outputs.outputs_all);
bootstrapping_sv = 0; % Flag for bootstrapping for SV
bounds_type_nexis_sv = 'unconstrained'; % Bounds type for NexIS:global parameters
niters_sv = 3; % Number of bootstrapped iterations for SV
outputs_all = struct;
numsims = length(mdlinstances)*size(nexsvnames,2);
sumtable_IbaHippInj = []; % One large table per study; predetermined to be 2
sumtable_Hurtado = []; % One large table per study; predetermined to be 2
for i = 1:length(mdlinstances)
    % Grab each previously run NexIS:global instance
    mdlinstance_i = mdlinstances{i};
    outputs_i = nexglob_outputs.outputs_all.(mdlinstance_i);
    outputs_i_inputstruct = outputs_i.nexis_global.Full;
    % Pull inputs from NexIS:global instance
    study_i = outputs_i_inputstruct.init.study;
    wdir_i = outputs_i_inputstruct.init.w_dir;
    volcorrect_i = outputs_i_inputstruct.init.volcorrect;
    % Note for below: use_dataspace was erroneously not saved earlier in 'init', is now
    if size(outputs_i_inputstruct.data,1) == size(outputs_i_inputstruct.init.C)
        usedataspace_i = 0;
    else
        usedataspace_i = 1;
    end
    excltpts_costfun_i = outputs_i_inputstruct.init.excltpts_costfun;
    % Output table row definitions
    sumtable_i = [];
    tablerowname = ['SV, ' study_i ', '];
    if wdir_i
        tablerowname = [tablerowname 'with dir., '];
    else
        tablerowname = [tablerowname 'no dir., '];
    end
    if isempty(excltpts_costfun_i)
        tablerowname = [tablerowname 'all tpts, '];
    else
        tablerowname = [tablerowname 'excl tpts, '];
    end
    if usedataspace_i
        tablerowname = [tablerowname 'data space'];
    else
        tablerowname = [tablerowname 'CCF space'];
    end
    for j = 1:size(nexsvnames,2)
        simno = size(nexsvnames,2)*(i-1) + j;
        fieldname_ij = [mdlinstance_i '_' nexsvnames{2,j}];
        fprintf('NexIS Wrapper Test 2.1, %d/%d\n',simno,numsims)
        fprintf('Simulation: %s\n',[tablerowname ', ' nexsvnames{1,j} ' ' nexsvnames{2,j}])
        outputs_nsv_ij = NexIS_SV('outputs_nexisglobal',outputs_i,...
                                  'study',study_i,...
                                  'w_dir',wdir_i,...
                                  'use_dataspace',usedataspace_i,...
                                  'volcorrect',volcorrect_i,...
                                  'excltpts_costfun',excltpts_costfun_i,...
                                  'datatype_nexis_sv',nexsvnames{1,j},...
                                  'datalist_nexis_sv',nexsvnames(2,j),...                
                                  'bootstrapping_nexis_sv',bootstrapping_sv,...              
                                  'bounds_type_nexis_sv',bounds_type_nexis_sv,...
                                  'niters_nexis_sv',niters_sv);
        outputs_all.(fieldname_ij) = outputs_nsv_ij;
        sumtable_ij = Output2Table(outputs_nsv_ij,0,'null','null'); % create table row
        sumtable_ij.Properties.RowNames{1} = ['Global' tablerowname(3:end)]; % label table row
        sumtable_ij.Properties.RowNames{2} = [tablerowname ', Factor ' num2str(j)]; % label table row
        if j == 1
            sumtable_i = [sumtable_i; sumtable_ij]; % add row to table
        else
            sumtable_i = [sumtable_i; sumtable_ij(2,:)]; % add row to table
        end
    end
    if strcmp(study_i,'IbaHippInj')
        sumtable_IbaHippInj = [sumtable_IbaHippInj; sumtable_i];
    elseif strcmp(study_i,'Hurtado')
        sumtable_Hurtado = [sumtable_Hurtado; sumtable_i];
    end
end

if writetofile
    % Non-programmatically splitting and saving tables by study
    tablename_IbaHippInj = [filename_out '_' 'IbaHippInj'];
    writetable(sumtable_IbaHippInj,[filepath_out filesep tablename_IbaHippInj '.csv'],'WriteRowNames',true)
    tablename_Hurtado = [filename_out '_' 'Hurtado'];
    writetable(sumtable_Hurtado,[filepath_out filesep tablename_Hurtado '.csv'],'WriteRowNames',true)
    % Save .mat file
    save([filepath_out filesep filename_out '.mat'],'outputs_all');
end 
%
%
%
%% 2.2 With bootstrapping of parameters, preloaded tau studies
%
%
%
% Loading previously run NexIS_global struct
rng(0); clear; clc;
filename_out = 'NexIS_Wrapper_SV_2-2_WithBootstrap'; % Name of output file
filename_in_glob = 'NexIS_Wrapper_Global_1-2_WithBootstrap'; % NexIS:global input file
filepath_in = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Load path, sims
filepath_out = filepath_in; % Save path
writetofile = 1;

% Pull name of one cell type
nexsvnames = cell(2,1);
nexsvnames{1} = 'Yao';
nexsvnames{2} = 'Oligo'; % Single SV factor for debugging purposes

% Running NexIS:SV (Exhaustive test of bootstrap)
nexglob_outputs = load([filepath_in filesep filename_in_glob '.mat'],'outputs_all');
mdlinstances = fieldnames(nexglob_outputs.outputs_all);
bootstrapping_sv = 1; % Flag for bootstrapping for SV
bounds_type_nexis_sv = 'CI_95'; % Bounds type for NexIS:global parameters
niters_sv = 3; % Number of bootstrapped iterations for SV
outputs_all = struct;
numsims = length(mdlinstances)*size(nexsvnames,2);
sumtable_IbaHippInj = []; % One large table per study; predetermined to be 2
sumtable_Hurtado = []; % One large table per study; predetermined to be 2
for i = 1:length(mdlinstances)
    % Grab each previously run NexIS:global instance
    mdlinstance_i = mdlinstances{i};
    outputs_i = nexglob_outputs.outputs_all.(mdlinstance_i);
    outputs_i_inputstruct = outputs_i.nexis_global.Full;
    % Pull inputs from NexIS:global instance
    study_i = outputs_i_inputstruct.init.study;
    wdir_i = outputs_i_inputstruct.init.w_dir;
    volcorrect_i = outputs_i_inputstruct.init.volcorrect;
    % Note for below: use_dataspace was erroneously not saved earlier in 'init', is now
    if size(outputs_i_inputstruct.data,1) == size(outputs_i_inputstruct.init.C)
        usedataspace_i = 0;
    else
        usedataspace_i = 1;
    end
    excltpts_costfun_i = outputs_i_inputstruct.init.excltpts_costfun;
    % Output table row definitions
    sumtable_i = [];
    tablerowname = ['SV, ' study_i ', '];
    if wdir_i
        tablerowname = [tablerowname 'with dir., '];
    else
        tablerowname = [tablerowname 'no dir., '];
    end
    if isempty(excltpts_costfun_i)
        tablerowname = [tablerowname 'all tpts, '];
    else
        tablerowname = [tablerowname 'excl tpts, '];
    end
    if usedataspace_i
        tablerowname = [tablerowname 'data space'];
    else
        tablerowname = [tablerowname 'CCF space'];
    end
    for j = 1:size(nexsvnames,2)
        simno = size(nexsvnames,2)*(i-1) + j;
        fieldname_ij = [mdlinstance_i '_' nexsvnames{2,j}];
        fprintf('NexIS Wrapper Test 2.2, %d/%d\n',simno,numsims)
        fprintf('Simulation: %s\n',[tablerowname ', ' nexsvnames{1,j} ' ' nexsvnames{2,j}])
        outputs_nsv_ij = NexIS_SV('outputs_nexisglobal',outputs_i,...
                                  'study',study_i,...
                                  'w_dir',wdir_i,...
                                  'use_dataspace',usedataspace_i,...
                                  'volcorrect',volcorrect_i,...
                                  'excltpts_costfun',excltpts_costfun_i,...
                                  'datatype_nexis_sv',nexsvnames{1,j},...
                                  'datalist_nexis_sv',nexsvnames(2,j),...                
                                  'bootstrapping_nexis_sv',bootstrapping_sv,...              
                                  'bounds_type_nexis_sv',bounds_type_nexis_sv,...
                                  'niters_nexis_sv',niters_sv);
        outputs_all.(fieldname_ij) = outputs_nsv_ij;
        sumtable_ij = Output2Table(outputs_nsv_ij,0,'null','null'); % create table row
        sumtable_ij.Properties.RowNames{1} = ['Global' tablerowname(3:end)]; % label table row
        sumtable_ij.Properties.RowNames{2} = [tablerowname ', Factor ' num2str(j)]; % label table row
        if j == 1
            sumtable_i = [sumtable_i; sumtable_ij]; % add row to table
        else
            sumtable_i = [sumtable_i; sumtable_ij(2,:)]; % add row to table
        end
    end
    if strcmp(study_i,'IbaHippInj')
        sumtable_IbaHippInj = [sumtable_IbaHippInj; sumtable_i];
    elseif strcmp(study_i,'Hurtado')
        sumtable_Hurtado = [sumtable_Hurtado; sumtable_i];
    end
end

if writetofile
    % Non-programmatically splitting and saving tables by study
    tablename_IbaHippInj = [filename_out '_' 'IbaHippInj'];
    writetable(sumtable_IbaHippInj,[filepath_out filesep tablename_IbaHippInj '.csv'],'WriteRowNames',true)
    tablename_Hurtado = [filename_out '_' 'Hurtado'];
    writetable(sumtable_Hurtado,[filepath_out filesep tablename_Hurtado '.csv'],'WriteRowNames',true)
    % Save .mat file
    save([filepath_out filesep filename_out '.mat'],'outputs_all');
end
%
%
%
%% 2.3 De novo NexIS:global, with and without bootstrap, Brundin a-syn studies only
%
%
%
% Running NexIS, testing out running NexIS:global within the NexIS:SV call
% as well as all combinations of bootstrapping (for completeness)
rng(0); clear; clc;
studylist = {'asyn_mouse','asyn_human'}; % ***Of preloaded, only these can be run with SV***
wdir = 1;
volcorrect = 1;
usedataspace = 0; % ***Required to be 0 for a-syn, all data in native connectome space***
param_init = [NaN,0,1,0.5]; % Initial fmincon parameter guesses; {gamma, alpha, beta, s}
ub = [Inf,Inf,Inf,1]; % Upper bounds for fmincon
lb = zeros(1,4); % Lower bounds for fmincon
excltpts_costfun = []; % Exclude selected time points from cost function
bootstrapping_glob = [0,1]; % Flag for bootstrapping
niters_glob = 3;
bootstrapping_sv = [0,1]; % Flag for bootstrapping for SV
bounds_type_nexis_sv = {'old','CI_95'}; % Bounds type for NexIS:global parameters
niters_sv = 3; % Number of bootstrapped iterations for SV

% Pull name of one cell type
nexsvnames = cell(2,1);
nexsvnames{1} = 'Zeisel';
nexsvnames{2} = 'HBNOR'; % Single SV factor for debugging purposes

% Table output parameters
writetofile = 1; % Create .csv from MATLAB table
filename_out = 'NexIS_Wrapper_SV_2-3_WithBootstrap_a-syn'; % Name of output file
filepath_out = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Save path

% Run model and create output tables for each dataset, if writetofile = 1
outputs_all = struct;
numsims = length(studylist)*length(bootstrapping_sv)*length(bootstrapping_glob)*size(nexsvnames,2);
for i = 1:length(studylist)
    study_i = studylist{i};
    tablename = [filename_out '_' study_i]; % Create one output table per study
    sumtable_i = [];
    for j = 1:length(bootstrapping_glob)
        for k = 1:length(bootstrapping_sv)
            for m = 1:size(nexsvnames,2)
                tablerowname = ['SV, ' study_i ', '];
                if bootstrapping_glob(j)
                    tablerowname = [tablerowname 'with b.s. glob, '];
                else
                    tablerowname = [tablerowname 'no b.s. glob, '];
                end
                if bootstrapping_sv(k)
                    tablerowname = [tablerowname 'with b.s. sv'];
                else
                    tablerowname = [tablerowname 'no b.s. sv'];
                end
                simno = size(nexsvnames,2)*length(bootstrapping_sv)*length(bootstrapping_glob)*(i-1)...
                    + size(nexsvnames,2)*(k-1) + size(nexsvnames,2)*length(bootstrapping_glob)*(j-1) + m;
                fprintf('NexIS Wrapper Test 2.3, %d/%d\n',simno,numsims)
                fprintf('Simulation: %s\n',[tablerowname ', ' nexsvnames{1,m} ' ' nexsvnames{2,m}])
                outputs_sv_ijkm = NexIS_SV('study',study_i,...
                                              'w_dir',wdir,...
                                              'use_dataspace',usedataspace,...
                                              'bootstrapping',bootstrapping_glob(j),...
                                              'niters',niters_glob,...
                                              'volcorrect',volcorrect,...
                                              'param_init',param_init,...
                                              'ub',ub,...
                                              'lb',lb,...
                                              'excltpts_costfun',excltpts_costfun,...
                                              'datatype_nexis_sv',nexsvnames{1,m},...
                                              'datalist_nexis_sv',nexsvnames(2,m),...                
                                              'bootstrapping_nexis_sv',bootstrapping_sv(k),...              
                                              'bounds_type_nexis_sv',bounds_type_nexis_sv{k},...
                                              'niters_nexis_sv',niters_sv);
                fieldname_ijkm = [study_i '_glob_bs_' num2str(bootstrapping_glob(j))...
                    '_sv_bs_' num2str(bootstrapping_sv(k)) '_' nexsvnames{2,m}];
                outputs_all.(fieldname_ijkm) = outputs_sv_ijkm;
                sumtable_ijkm = Output2Table(outputs_sv_ijkm,0,'null','null'); % create table row
                sumtable_ijkm.Properties.RowNames{1} = ['Global' tablerowname(3:end)]; % label table row
                sumtable_ijkm.Properties.RowNames{2} = [tablerowname ', Factor ' num2str(m)]; % label table row
                if m == 1
                    sumtable_i = [sumtable_i; sumtable_ijkm]; % add row to table
                else
                    sumtable_i = [sumtable_i; sumtable_ijkm(2,:)]; % add row to table
                end
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
%
%
%
%% 2.4 No bootstrapping, user-specified connectome & pathology
%
%
%
% Loading previously run NexIS_global struct. Not required to do this to run
% NexIS:SV, but recommended for efficiency reasons, particularly if running
% through multiple factors (i.e., genes, cell types) for the same pathology
% dataset
rng(0); clear; clc;
filename_out = 'NexIS_Wrapper_SV_2-4_NoBootstrap_user-specified'; % Name of output file
filename_in_glob = 'NexIS_Wrapper_Global_1-4_NoBootstrap_user-specified'; % NexIS:global input file
filepath_in = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Load path, sims
filepath_out = filepath_in; % Save path
writetofile = 1;

% Define U vectors a priori instead of calling one
nexsvnames = cell(2,2);
nexsvnames{1,1} = 'User_specified';
nexsvnames{2,1} = rand(426,1); % Know n_ROI ahead of time
nexsvnames{1,2} = 'User_specified';
nexsvnames{2,2} = randn(426,1);

% Running NexIS:SV (Exhaustive test of no bootstrap)
nexglob_outputs = load([filepath_in filesep filename_in_glob '.mat'],'outputs_all');
mdlinstances = fieldnames(nexglob_outputs.outputs_all);
bootstrapping_sv = 0; % Flag for bootstrapping for SV
bounds_type_nexis_sv = 'unconstrained'; % Bounds type for NexIS:global parameters
niters_sv = 3; % Number of bootstrapped iterations for SV
outputs_all = struct;
numsims = length(mdlinstances)*size(nexsvnames,2);
for i = 1:length(mdlinstances)
    % Grab each previously run NexIS:global instance
    mdlinstance_i = mdlinstances{i};
    tablename_i = [filename_out '_' mdlinstance_i]; % Create one output table per output instance
    outputs_i = nexglob_outputs.outputs_all.(mdlinstance_i);
    outputs_i_inputstruct = outputs_i.nexis_global.Full;
    % Pull inputs from NexIS:global instance
    study_i = outputs_i_inputstruct.init.study;
    C_i = outputs_i_inputstruct.init.C;
    data_i = outputs_i_inputstruct.data;
    seed_i = outputs_i_inputstruct.init.seed;
    if isnan(seed_i) % Check if simulation was run from baseline
        data_i = [outputs_i_inputstruct.baseline, data_i];
    end
    ts_i = outputs_i_inputstruct.time_stamps;
    wdir_i = outputs_i_inputstruct.init.w_dir;
    volcorrect_i = outputs_i_inputstruct.init.volcorrect;
    usedataspace_i = outputs_i_inputstruct.init.use_dataspace;
    excltpts_costfun_i = outputs_i_inputstruct.init.excltpts_costfun;
    bootstrapping_glob_i = outputs_i_inputstruct.init.bootstrapping;
    % Output table row definitions
    sumtable_i = [];
    tablerowname = ['SV, ' study_i '_Data'];
    for j = 1:size(nexsvnames,2)
        simno = size(nexsvnames,2)*(i-1) + j;
        fieldname_ij = [mdlinstance_i '_' nexsvnames{1,j} '_Factor_' num2str(j)];
        fprintf('NexIS Wrapper Test 2.4, %d/%d\n',simno,numsims)
        fprintf('Simulation: %s\n',[tablerowname ', ' nexsvnames{1,j} '_Factor_' num2str(j)])
        outputs_nsv_ij = NexIS_SV('outputs_nexisglobal',outputs_i,...                                  
                                  'C',C_i,...
                                  'data',data_i,...
                                  'seed',seed_i,...
                                  'tpts',ts_i,...
                                  'w_dir',wdir_i,...
                                  'study', study_i,...
                                  'use_dataspace',usedataspace_i,...
                                  'bootstrapping',bootstrapping_glob_i,...
                                  'volcorrect',volcorrect_i,...
                                  'datatype_nexis_sv',nexsvnames{1,j},...
                                  'datalist_nexis_sv',nexsvnames{2,j},...                
                                  'bootstrapping_nexis_sv',bootstrapping_sv,...              
                                  'bounds_type_nexis_sv',bounds_type_nexis_sv,...
                                  'niters_nexis_sv',niters_sv);
        outputs_all.(fieldname_ij) = outputs_nsv_ij;
        sumtable_ij = Output2Table(outputs_nsv_ij,0,'null','null'); % create table row
        sumtable_ij.Properties.RowNames{1} = ['Global' tablerowname(3:end)]; % label table row
        sumtable_ij.Properties.RowNames{2} = [tablerowname ', Factor ' num2str(j)]; % label table row
        if j == 1
            sumtable_i = [sumtable_i; sumtable_ij]; % add row to table
        else
            sumtable_i = [sumtable_i; sumtable_ij(2,:)]; % add row to table
        end
    end
    if writetofile
        writetable(sumtable_i,[filepath_out filesep tablename_i '.csv'],'WriteRowNames',true)
    end
end

if writetofile
    % Save .mat file
    save([filepath_out filesep filename_out '.mat'],'outputs_all');
end 
%
%
%
%% 2.5 With bootstrapping, user-specified connectome & pathology
%
%
%
% Loading previously run NexIS_global struct. Not required to do this to run
% NexIS:SV, but recommended for efficiency reasons, particularly if running
% through multiple factors (i.e., genes, cell types) for the same pathology
% dataset
rng(0); clear; clc;
filename_out = 'NexIS_Wrapper_SV_2-5_WithBootstrap_user-specified'; % Name of output file
filename_in_glob = 'NexIS_Wrapper_Global_1-5_WithBootstrap_user-specified'; % NexIS:global input file
filepath_in = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Load path, sims
filepath_out = filepath_in; % Save path
writetofile = 1;

% Define U vectors a priori instead of calling one
nexsvnames = cell(2,2);
nexsvnames{1,1} = 'User_specified';
nexsvnames{2,1} = rand(426,1); % Know n_ROI ahead of time
nexsvnames{1,2} = 'User_specified';
nexsvnames{2,2} = randn(426,1);

% Running NexIS:SV (Exhaustive test of no bootstrap)
nexglob_outputs = load([filepath_in filesep filename_in_glob '.mat'],'outputs_all');
mdlinstances = fieldnames(nexglob_outputs.outputs_all);
bootstrapping_sv = 1; % Flag for bootstrapping for SV
bounds_type_nexis_sv = 'unconstrained'; % Bounds type for NexIS:global parameters
niters_sv = 3; % Number of bootstrapped iterations for SV
outputs_all = struct;
numsims = length(mdlinstances)*size(nexsvnames,2);
sumtable_IbaHippInj = []; % One large table per study; predetermined to be 2
sumtable_Hurtado = []; % One large table per study; predetermined to be 2
for i = 1:length(mdlinstances)
    % Grab each previously run NexIS:global instance
    mdlinstance_i = mdlinstances{i};
    tablename_i = [filename_out '_' mdlinstance_i]; % Create one output table per output instance
    outputs_i = nexglob_outputs.outputs_all.(mdlinstance_i);
    outputs_i_inputstruct = outputs_i.nexis_global.Full;
    % Pull inputs from NexIS:global instance
    study_i = outputs_i_inputstruct.init.study;
    C_i = outputs_i_inputstruct.init.C;
    data_i = outputs_i_inputstruct.data;
    seed_i = outputs_i_inputstruct.init.seed;
    if isnan(seed_i) % Check if simulation was run from baseline
        data_i = [outputs_i_inputstruct.baseline, data_i];
    end
    ts_i = outputs_i_inputstruct.time_stamps;
    wdir_i = outputs_i_inputstruct.init.w_dir;
    volcorrect_i = outputs_i_inputstruct.init.volcorrect;
    usedataspace_i = outputs_i_inputstruct.init.use_dataspace;
    excltpts_costfun_i = outputs_i_inputstruct.init.excltpts_costfun;
    bootstrapping_glob_i = outputs_i_inputstruct.init.bootstrapping;
    % Output table row definitions
    sumtable_i = [];
    tablerowname = ['SV, ' study_i '_Data'];
    for j = 1:size(nexsvnames,2)
        simno = size(nexsvnames,2)*(i-1) + j;
        fieldname_ij = [mdlinstance_i '_' nexsvnames{1,j} '_Factor_' num2str(j)];
        fprintf('NexIS Wrapper Test 2.5, %d/%d\n',simno,numsims)
        fprintf('Simulation: %s\n',[tablerowname ', ' nexsvnames{1,j} '_Factor_' num2str(j)])
        outputs_nsv_ij = NexIS_SV('outputs_nexisglobal',outputs_i,...                                  
                                  'C',C_i,...
                                  'data',data_i,...
                                  'seed',seed_i,...
                                  'tpts',ts_i,...
                                  'w_dir',wdir_i,...
                                  'study', study_i,...
                                  'use_dataspace',usedataspace_i,...
                                  'bootstrapping',bootstrapping_glob_i,...
                                  'volcorrect',volcorrect_i,...
                                  'datatype_nexis_sv',nexsvnames{1,j},...
                                  'datalist_nexis_sv',nexsvnames{2,j},...                
                                  'bootstrapping_nexis_sv',bootstrapping_sv,...              
                                  'bounds_type_nexis_sv',bounds_type_nexis_sv,...
                                  'niters_nexis_sv',niters_sv);
        outputs_all.(fieldname_ij) = outputs_nsv_ij;
        sumtable_ij = Output2Table(outputs_nsv_ij,0,'null','null'); % create table row
        sumtable_ij.Properties.RowNames{1} = ['Global' tablerowname(3:end)]; % label table row
        sumtable_ij.Properties.RowNames{2} = [tablerowname ', Factor ' num2str(j)]; % label table row
        if j == 1
            sumtable_i = [sumtable_i; sumtable_ij]; % add row to table
        else
            sumtable_i = [sumtable_i; sumtable_ij(2,:)]; % add row to table
        end
    end
    if writetofile
        writetable(sumtable_i,[filepath_out filesep tablename_i '.csv'],'WriteRowNames',true)
    end
end

if writetofile
    % Save .mat file
    save([filepath_out filesep filename_out '.mat'],'outputs_all');
end 
%
%
%
%% 2.6 Multiple factors, with and without PCA
%
%
%
% Running NexIS, testing out running NexIS:global within the NexIS:SV call
% as well as all combinations of bootstrapping (for completeness)
rng(0); clear; clc;
studylist = {'IbaHippInj'}; % ***Code won't work with more than one study at a time***
wdir = 1;
volcorrect = 1;
usedataspace = 1;
param_init = [NaN,0,1,0.5]; % Initial fmincon parameter guesses; {gamma, alpha, beta, s}
ub = [Inf,Inf,Inf,1]; % Upper bounds for fmincon
lb = zeros(1,4); % Lower bounds for fmincon
bootstrapping_glob = [0,1]; % Flag for bootstrapping
niters_glob = 3;
bootstrapping_sv = [0,1]; % Flag for bootstrapping for SV
bounds_type_nexis_sv = {'old','CI_95'}; % Bounds type for NexIS:global parameters
niters_sv = 3; % Number of bootstrapped iterations for SV
usepca = [0,1]; % Flag for using PC 1 of cell type distributions

% Pull name of one cell type
fp_in_sv = '~/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
load([fp_in_sv filesep 'CellTypeMaps.mat'],'CellTypeMaps');
nexsvnames = cell(2,2);
nexsvnames(1,:) = {'Yao','User_specified'};
nexsvnames{2,1} = {'Pvalb','Sst','Vip'};
interneuron_bool = ismember(CellTypeMaps.(nexsvnames{1,1}).classkey,nexsvnames{2,1});
nexsvnames{2,2} = CellTypeMaps.(nexsvnames{1,1}).maps(:,interneuron_bool);

% Table output parameters
writetofile = 1; % Create .csv from MATLAB table
filename_out = 'NexIS_Wrapper_SV_2-6_WithPCA'; % Name of output file
filepath_out = '~/Documents/MATLAB/Nexis_Project/Results_Files_NexISWrapper'; % Save path

% Run model and create output tables for each dataset, if writetofile = 1
outputs_all = struct;
numsims = length(studylist)*length(bootstrapping_glob)*length(usepca)*size(nexsvnames,2);
sumtable_noPCA = []; % w/ and w/o PCA have to be separate tables
sumtable_withPCA = [];
for k = 1:length(usepca)
    sumtable_k = [];
    for i = 1:length(studylist) % ***Code won't work with more than one study***
        study_i = studylist{i};
        tablename = [filename_out '_' study_i '_PCA_' num2str(k)]; % Create one output table per study
        for j = 1:length(bootstrapping_glob) % Use bootstrapping parameter for global and SV together  
            for m = 1:size(nexsvnames,2)
                tablerowname = ['SV, ' study_i ', '];
                if bootstrapping_glob(j)
                    tablerowname = [tablerowname 'with b.s., '];
                else
                    tablerowname = [tablerowname 'no b.s., '];
                end
                if usepca(k)
                    tablerowname = [tablerowname 'with PCA'];
                else
                    tablerowname = [tablerowname 'no PCA'];
                end
                simno = size(nexsvnames,2)*length(studylist)*length(bootstrapping_glob)*(k-1)...
                    + size(nexsvnames,2)*(j-1) + size(nexsvnames,2)*length(bootstrapping_glob)*(i-1) + m;
                fprintf('NexIS Wrapper Test 2.6, %d/%d\n',simno,numsims)
                fprintf('Simulation: %s\n',[tablerowname ', ' nexsvnames{1,m} ' ' 'Interneuron Test ' num2str(m)])
                outputs_sv_ijkm = NexIS_SV('study',study_i,...
                                              'w_dir',wdir,...
                                              'use_dataspace',usedataspace,...
                                              'bootstrapping',bootstrapping_glob(j),...
                                              'niters',niters_glob,...
                                              'volcorrect',volcorrect,...
                                              'param_init',param_init,...
                                              'ub',ub,...
                                              'lb',lb,...
                                              'datatype_nexis_sv',nexsvnames{1,m},...
                                              'datalist_nexis_sv',nexsvnames{2,m},...                
                                              'bootstrapping_nexis_sv',bootstrapping_sv(j),...        
                                              'bounds_type_nexis_sv',bounds_type_nexis_sv{j},...
                                              'niters_nexis_sv',niters_sv,...
                                              'datapca_nexis_sv',usepca(k));
                fieldname_ijkm = [study_i '_bs_' num2str(bootstrapping_glob(j))...
                    '_Interneuron_' num2str(m) '_PCA_' num2str(usepca(k))];
                outputs_all.(fieldname_ijkm) = outputs_sv_ijkm;
                sumtable_ijkm = Output2Table(outputs_sv_ijkm,0,'null','null'); % create table row
                sumtable_ijkm.Properties.RowNames{1} = ['Global' tablerowname(3:end)]; % label table row
                sumtable_ijkm.Properties.RowNames{2} = [tablerowname ', Factor ' num2str(m)]; % label table row
                if m == 1
                    sumtable_k = [sumtable_k; sumtable_ijkm]; % add row to table
                else
                    sumtable_k = [sumtable_k; sumtable_ijkm(2,:)]; % add row to table
                end
            end
        end
    end
    if usepca(k)
        sumtable_withPCA = [sumtable_withPCA; sumtable_k];
    else
        sumtable_noPCA = [sumtable_noPCA; sumtable_k];
    end
end
if writetofile
    % Non-programmatically splitting and saving tables by study
    tablename_withPCA = [filename_out '_' 'withPCA'];
    writetable(sumtable_withPCA,[filepath_out filesep tablename_withPCA '.csv'],'WriteRowNames',true)
    tablename_noPCA = [filename_out '_' 'noPCA'];
    writetable(sumtable_noPCA,[filepath_out filesep tablename_noPCA '.csv'],'WriteRowNames',true)
    % Save .mat file
    save([filepath_out filesep filename_out '.mat'],'outputs_all');
end
% Consistency check appeared to work - both interneuron tests yielded same
% results (1/2/25)
%
%
%
%% 3 Relevant plotting of outputs
%
%
%
% Deal with later