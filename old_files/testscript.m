outputs1 = stdNDM_mouse('study','DS9','bootstrapping',1,'niters',5,'verbose',0,'fmindisplay',0);

% clc;
% homeo = {'P2ry12','Cx3cr1','Fcrls','Olfml3','Hexb','Siglech','Sox5','Jun'};
outputs2 = eNDM_mouse('study','DS9','bootstrapping',0,'niters',5,...
    'bootstrapping_endm',0,'niters_endm',5,'verbose',0,'fmindisplay',0,...
    'verbose_endm',0,'fmindisplay_endm',0,'datalist_endm',{'Trem2'});
Output2Table(outputs2,1)
BootstrappingPlotter(outputs2);
CorrelationPlotter(outputs2);
outputs3 = eNDM_mouse('study','DS9','bootstrapping',1,'niters',2,...
    'bootstrapping_endm',1,'niters_endm',2,'verbose',0,'fmindisplay',0,...
    'verbose_endm',0,'fmindisplay_endm',0,'outputs_ndm',outputs1);

homeo = {'P2ry12','Cx3cr1','Fcrls','Olfml3','Hexb','Siglech','Sox5','Jun'};
% BootstrappingPlotter(outputs1);
BootstrappingPlotter(outputs2);

rng(3);
outputs3 = eNDM_mouse('study','DS9','bootstrapping',1,'niters',3,'w_dir',1,...
    'bootstrapping_endm',1,'niters_endm',3,'verbose',0,'fmindisplay',0,...
    'verbose_endm',0,'fmindisplay_endm',0,'datalist_endm',randperm(3855,3),'datapca_endm',1);
Output2Table(outputs3);
BootstrappingPlotter(outputs3);
CorrelationPlotter(outputs3);