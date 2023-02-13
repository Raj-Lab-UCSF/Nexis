% homeo = [2344, 801, 1162, 2305, 1540, 3062, 3240, 1717];
rng(0);
% studylist = {'asyn_mouse', 'asyn_human'};
studylist = {'Henderson'};
niters = 1;
R2s = zeros(1,niters);
for i = 1:length(studylist)
    outputs_ndm = stdNDM_mouse('study',studylist{i},'bootstrapping',0,'niters',niters,...
        'w_dir',1,'volcorrect',0,'param_init',[NaN,0,1,0.5],'ub',[Inf,Inf,Inf,1],...
        'lb',zeros(1,4));
%     for j = 1:niters
%         outputs_endm = eNDM_mouse('outputs_ndm',outputs_ndm,'study',studylist{i},...
%             'w_dir',1,'bootstrapping_endm',0,'volcorrect',1,...
%             'bounds_type_endm','CI_95','datatype_endm','ct_zeisel',...
%             'datalist_endm',{'random'});
%         R2s(j) = outputs_endm.endm.Full.results.lm_Rsquared_adj;
%         Output2Table(outputs_endm,1,[studylist{i} '_output_table'])
%         BootstrappingPlotter(outputs_endm);
%         CorrelationPlotter(outputs_endm);
%     end
end
pred = outputs_ndm.ndm.Full.predicted;
load('/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse/regionvoxels.mat','voxels');
voxels_ = [voxels; voxels];
concentration_sum = sum(pred)
mass_sum = voxels_.' * pred


%%

outputs6 = eNDM_mouse('outputs_ndm',outputs4,'study','IbaHippInj','bootstrapping_endm',...
    0,'niters_endm',1,'w_dir',0,'exclseed_costfun',1);
z = Output2Table(outputs6); 
studies = {'asyn_mouse','asyn_human'};
genesets = {{'Trem2'},{'P2ry12','Cx3cr1','Fcrls','Olfml3','Hexb','Siglech','Sox5','Jun'},...
    {'P2ry12','Cx3cr1','Fcrls','Olfml3','Hexb','Siglech','Sox5','Jun','Trem2'}};
datapca = [0 1 1];
setnames = {'Trem2_only','Homeostatic','Homeostatic+Trem2'};
for i = 1:length(studies)
    outputs1 = stdNDM_mouse('study',studies{i},'bootstrapping',1,'niters',100,'verbose',0,'fmindisplay',0);
    for j = 1:length(genesets)
        outputs2 = eNDM_mouse('outputs_ndm',outputs1,'study',studies{i},...
            'bootstrapping_endm',1,'niters_endm',100,'verbose',0,'fmindisplay',0,...
            'verbose_endm',0,'fmindisplay_endm',0,'datalist_endm',genesets{j},'datapca_endm',datapca(j));
        Output2Table(outputs2,1,[studies{i} '_' setnames{j} '_output_table'])
        BootstrappingPlotter(outputs2);
        CorrelationPlotter(outputs2);
        save([studies{i} '_' setnames{j} '_outputs.mat'],'outputs2','-v7.3');
    end
end



% 
% outputs3 = eNDM_mouse('study','DS9','bootstrapping',1,'niters',2,...
%     'bootstrapping_endm',1,'niters_endm',2,'verbose',0,'fmindisplay',0,...
%     'verbose_endm',0,'fmindisplay_endm',0,'outputs_ndm',outputs1);
% 
% homeo = {'P2ry12','Cx3cr1','Fcrls','Olfml3','Hexb','Siglech','Sox5','Jun'};
% % BootstrappingPlotter(outputs1);
% BootstrappingPlotter(outputs2);
% 
% rng(3);
% outputs3 = eNDM_mouse('study','DS9','bootstrapping',1,'niters',3,'w_dir',1,...
%     'bootstrapping_endm',1,'niters_endm',3,'verbose',0,'fmindisplay',0,...
%     'verbose_endm',0,'fmindisplay_endm',0,'datalist_endm',randperm(3855,3),'datapca_endm',1);
% Output2Table(outputs3);
% BootstrappingPlotter(outputs3);
% CorrelationPlotter(outputs_endm);

% For Sam - This is a test for the alpha synuclein data, running a PCA on
% the homeostatic markers of microglia and using the first component as a
% predictor. It runs both the standard NDM and the eNDM with or without
% bootstrapping.
% homeo = {'P2ry12','Cx3cr1','Fcrls','Olfml3','Hexb','Siglech','Sox5','Jun'};
% rescalevals = [1,2];
% alphavals = [0.25 0.5];
% betavals = [];
% studies = {'asyn_mouse','asyn_human'};
% outputcell = cell(8,2); % 
% index = 1;
% for k = 1:length(studies)
%     for i = 1:length(rescalevals)
%         for j = 1:length(alphavals)
%             outputs_asynmouse = stdNDM_mouse('study',studies{k},'bootstrapping',0,...
%                 'verbose',0,'fmindisplay',0,'normtype','none','param_init',...
%                 [rescalevals(i),alphavals(j),8,0.5],'algo','interior-point');
% %             flname = sprintf('param_init_test_%s_rescaleval_%d_alphaval_%d',studies{k},rescalevals(i),alphavals(j));
%             outputstable = Output2Table(outputs_asynmouse);            
%             % add entries to cell array from outputstable here
%             outputcell{index,1} = outputstable.('alpha (Mean)')(1);
%             outputcell{index,2} = outputstable.('beta (Mean)')(1);
%             index = index + 1;
%             clear outputs_asynmouse outputstable
%         end
%     end
% end
% % Output2Table(outputs_asynmouse)
% % BootstrappingPlotter(outputs_asynmouse);
% % CorrelationPlotter(outputs_asynmouse);
