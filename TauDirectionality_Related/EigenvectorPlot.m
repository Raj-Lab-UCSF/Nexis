function EigenvectorPlot(studyname_,C_,matdir)
%% 2.1 Figure 1, scatterplots
% Hurtado & IbaHippInj
hurtado_end = mousedata_struct.Hurtado.data(:,end);
ibahippinj_end = mousedata_struct.IbaStrInj.data(:,end);

outdeg = sum(C_ret,2); 
outdeg_hurtado = CCFToData(outdeg,'Hurtado',matdir);
outdeg_ibahippinj = CCFToData(outdeg,'IbaStrInj',matdir);

indeg = sum(C_ret,1).'; 
indeg_hurtado = CCFToData(indeg,'Hurtado',matdir);
indeg_ibahippinj = CCFToData(indeg,'IbaStrInj',matdir);

u1_ret = v_ret(:,1); 
u1_ret_hurtado = CCFToData(u1_ret,'Hurtado',matdir);
u1_ret_ibahippinj = CCFToData(u1_ret,'IbaStrInj',matdir);

u1_ant = v_ant(:,1); 
u1_ant_hurtado = CCFToData(u1_ant,'Hurtado',matdir);
u1_ant_ibahippinj = CCFToData(u1_ant,'IbaStrInj',matdir);

ylim1 = [0 max(ibahippinj_end)];
ylim2 = [0 max(hurtado_end)];
xlim1 = [0 max(outdeg_ibahippinj)];
xlim2 = [0 max(indeg_ibahippinj)];
xlim3 = [0 max(u1_ret_ibahippinj)];
xlim4 = [0 max(u1_ant_ibahippinj)];

figure('Units','inches','Position',[0 0 11 6]); 
tiledlayout(2,4,'TileSpacing','compact','Padding','tight');
nexttile;
scatter(outdeg_ibahippinj,ibahippinj_end,'bo','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim1); xlim(xlim1);
xticks([xlim1(1) (xlim1(1)+xlim1(2))/2 xlim1(2)]);
yticks([ylim1(1) (ylim1(1)+ylim1(2))/2 ylim1(2)]);
xticklabels({''});
ylabel('IbaHippInj Pathology')
text(0.6*xlim1(2),0.1*ylim1(2),sprintf('R = %.2f',corr(outdeg_ibahippinj,ibahippinj_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(indeg_ibahippinj,ibahippinj_end,'ro','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim1); xlim(xlim2); xticklabels({''}); yticklabels({''}); 
xticks([xlim2(1) (xlim2(1)+xlim2(2))/2 xlim2(2)]);
yticks([ylim1(1) (ylim1(1)+ylim1(2))/2 ylim1(2)]);
text(0.6*xlim2(2),0.1*ylim1(2),sprintf('R = %.2f',corr(indeg_ibahippinj,ibahippinj_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(u1_ret_ibahippinj,ibahippinj_end,'bd','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim1); xlim(xlim3); xticklabels({''}); yticklabels({''});
xticks([xlim3(1) (xlim3(1)+xlim3(2))/2 xlim3(2)]);
yticks([ylim1(1) (ylim1(1)+ylim1(2))/2 ylim1(2)]);
text(0.6*xlim3(2),0.1*ylim1(2),sprintf('R = %.2f',corr(u1_ret_ibahippinj,ibahippinj_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(u1_ant_ibahippinj,ibahippinj_end,'rd','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim1); xlim(xlim4); xticklabels({''}); yticklabels({''});
xticks([xlim4(1) (xlim4(1)+xlim4(2))/2 xlim4(2)]);
yticks([ylim1(1) (ylim1(1)+ylim1(2))/2 ylim1(2)]);
text(0.6*xlim4(2),0.1*ylim1(2),sprintf('R = %.2f',corr(u1_ant_ibahippinj,ibahippinj_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(outdeg_hurtado,hurtado_end,'bo','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim2); xlim(xlim1);
xticks([xlim1(1) (xlim1(1)+xlim1(2))/2 xlim1(2)]);
yticks([ylim2(1) (ylim2(1)+ylim2(2))/2 ylim2(2)]);
xtickformat('%.1d');
xlabel('Out Degree'); ylabel('Hurtado Pathology');
text(0.6*xlim1(2),0.1*ylim2(2),sprintf('R = %.2f',corr(outdeg_hurtado,hurtado_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(indeg_hurtado,hurtado_end,'ro','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim2); xlim(xlim2); yticklabels({''});
xticks([xlim2(1) (xlim2(1)+xlim2(2))/2 xlim2(2)]);
yticks([ylim2(1) (ylim2(1)+ylim2(2))/2 ylim2(2)]);
xtickformat('%.1d');
xlabel('In Degree'); 
text(0.6*xlim2(2),0.1*ylim2(2),sprintf('R = %.2f',corr(indeg_hurtado,hurtado_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(u1_ret_hurtado,hurtado_end,'bd','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim2); xlim(xlim3); yticklabels({''}); 
xticklabels({'0',num2str((xlim3(1)+xlim3(2))/2,'%.2f'),num2str(xlim3(2),'%.2f')})
xticks([xlim3(1) (xlim3(1)+xlim3(2))/2 xlim3(2)]);
yticks([ylim2(1) (ylim2(1)+ylim2(2))/2 ylim2(2)]);
xlabel('u1_r_e_t');
text(0.6*xlim3(2),0.1*ylim2(2),sprintf('R = %.2f',corr(u1_ret_hurtado,hurtado_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(u1_ant_hurtado,hurtado_end,'rd','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim2); xlim(xlim4); yticklabels({''});
xticklabels({'0',num2str((xlim4(1)+xlim4(2))/2,'%.2f'),num2str(xlim4(2),'%.2f')})
xticks([xlim4(1) (xlim4(1)+xlim4(2))/2 xlim4(2)]);
yticks([ylim2(1) (ylim2(1)+ylim2(2))/2 ylim2(2)]);
xlabel('u1_a_n_t'); 
text(0.6*xlim4(2),0.1*ylim2(2),sprintf('R = %.2f',corr(u1_ant_hurtado,hurtado_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

% print([figdir filesep 'NoModelScatterplots'],'-dtiffn','-r300'); close;

end