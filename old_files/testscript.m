% outputs1 = stdNDM_mouse('study','DS9','bootstrapping',1,'niters',5,'verbose',0,'fmindisplay',0);
% 
% % clc;
% homeo = {'P2ry12','Cx3cr1','Fcrls','Olfml3','Hexb','Siglech','Sox5','Jun'};
% homeo = [2344, 801, 1162, 2305, 1540, 3062, 3240, 1717];
% outputs2 = eNDM_mouse('study','DS9','bootstrapping',0,'niters',5,...
%     'bootstrapping_endm',0,'niters_endm',5,'verbose',0,'fmindisplay',0,...
%     'verbose_endm',0,'fmindisplay_endm',0);
% Output2Table(outputs2)
% BootstrappingPlotter(outputs2);
% CorrelationPlotter(outputs2);
% 
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
% CorrelationPlotter(outputs3);

% For Sam - This is a test for the alpha synuclein data, running a PCA on
% the homeostatic markers of microglia and using the first component as a
% predictor. It runs both the standard NDM and the eNDM with or without
% bootstrapping.
homeo = {'P2ry12','Cx3cr1','Fcrls','Olfml3','Hexb','Siglech','Sox5','Jun'};
rescalevals = [1,2,4,8,12];
alphavals = [0.25 0.5 1];
betavals = [];
studies = {'asyn_mouse','asyn_human'};
% outputcell = cell((length(studies)*length..., however many outputs you want)
index = 1;
for k = 1:length(studies)
    for i = 1:length(rescalevals)
        for j = 1:length(alphavals)
            outputs_asynmouse = eNDM_mouse('study',studies{k},'bootstrapping',0,...
                'verbose',0,'fmindisplay',0,'normtype','none','param_init',...
                [rescalevals(i),alphavals(j),8,0.5],'algo','interior-point');
%             flname = sprintf('param_init_test_%s_rescaleval_%d_alphaval_%d',studies{k},rescalevals(i),alphavals(j));
            outputstable = Output2Table(outputs_asynmouse)            
            % add entries to cell array from outputstable here
            %outputcell{index, column of interest)  = value of interest
            index = index + 1;
        end
    end
end
% Output2Table(outputs_asynmouse)
% BootstrappingPlotter(outputs_asynmouse);
% CorrelationPlotter(outputs_asynmouse);
