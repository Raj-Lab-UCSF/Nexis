rng(0); clear; clc;
studylist = {'DS4'};
matdir = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
% x = [0.536129657164085,0.25,1.89,0.500000000000000];
% niters = 10;
% R2s = zeros(1,niters);
for i = 1:length(studylist)
    
    % outputs_ng_ccf = NexIS_global('study',studylist{i},'bootstrapping',0,...
    %     'w_dir',0,'volcorrect',0,'param_init',x,'ub',x,...
    %     'lb',x,'use_dataspace',0);    
    outputs_ng_ccf = NexIS_global('study',studylist{i},'bootstrapping',0,...
        'w_dir',1,'volcorrect',1,'param_init',[NaN,0,1,0.5],'ub',[Inf,Inf,Inf,1],...
        'lb',zeros(1,4),'use_dataspace',0);
    
    % outputs_ng_ds = NexIS_global('study',studylist{i},'bootstrapping',0,...
    %     'w_dir',1,'volcorrect',1,'param_init',[NaN,0,1,0.5],'ub',[Inf,Inf,Inf,1],...
    %     'lb',zeros(1,4),'use_dataspace',1);
end


% C = Connectomes.default;
% C = C/max(C(:));
% seed = seed426.DS4;
% x0 = seed * 0.8;
% time_stamps = [4,8,12];
% y = eNDM_general_dir(x0,time_stamps,C,zeros(426,1),0.2,3.5,0.65,0,0,0,'analytic',0);
% corr(y(:), data426.DS4(:),'rows','complete')^2