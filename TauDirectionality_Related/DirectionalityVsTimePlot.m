function DirectionalityVsTimePlot(outstruct,use_s,titlestr)
     
studynames = fieldnames(outstruct);
studynames(ismember(studynames,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
studylabels = cellfun(@(x)strrep(x,'_',' '),studynames,'UniformOutput',0);
modelnames = fieldnames(outstruct.(studynames{1}));

tptnames = fieldnames(outstruct.(studynames{1}).(modelnames{1}));
s_deltaRmat = NaN(length(studynames),length(tptnames)); % 4 models
tptmat = s_deltaRmat;
for i = 1:size(s_deltaRmat,1)
    if use_s
        resstruct = outstruct.(studynames{i}).('fit_s');
        ylab = 's';
        for k = 1:length(tptnames)
            resstruct_tpt = resstruct.(tptnames{k}).nexis_global.Full;
            tptmat(i,k) = resstruct_tpt.time_stamps;
            s_deltaRmat(i,k) = resstruct_tpt.param_fit(4);
        end
    else
        resstruct_ret = outstruct.(studynames{i}).('ret');
        resstruct_ant = outstruct.(studynames{i}).('ant');
        ylab = '\DeltaR_d_i_r';
        for k = 1:length(tptnames)
            resstruct_ret_tpt = resstruct_ret.(tptnames{k}).nexis_global.Full;
            tptmat(i,k) = resstruct_ret_tpt.time_stamps;
            datavec_ret = resstruct_ret_tpt.data(:);
            predvec_ret = resstruct_ret_tpt.predicted(:);
            corr_ret = corr(datavec_ret,predvec_ret,'rows','complete');
            resstruct_ant_tpt = resstruct_ant.(tptnames{k}).nexis_global.Full;
            datavec_ant = resstruct_ant_tpt.data(:);
            predvec_ant = resstruct_ant_tpt.predicted(:);
            corr_ant = corr(datavec_ant,predvec_ant,'rows','complete');
            s_deltaRmat(i,k) = corr_ret - corr_ant;
        end
    end
end

cmap = hsv(length(studynames));
shapes = {'o','s','d','^','v','<','>','p','h','+','x'};
figure('Units','inches','Position',[0 0 11 10]); hold on;
plothands = {};
for i = 1:length(studynames)
    s = scatter(tptmat(i,:),s_deltaRmat(i,:),75,shapes{i},...
        'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),...
        'MarkerFaceAlpha',0.3);
    plothands = [plothands,s];
end
lm = fitlm(tptmat(:),s_deltaRmat(:));
x_lm = linspace(0,10,100).'; 
[y_lm, y_ci] = predict(lm, x_lm);
plot(x_lm,y_ci(:,1),'k:'); plot(x_lm,y_ci(:,2),'k:'); 
fill([x_lm; flipud(x_lm)],[y_ci(:,1); flipud(y_ci(:,2))],[1 0 0.25],...
    'EdgeColor','none','FaceAlpha',0.15);
h = plot(x_lm,y_lm,'k','LineWidth',2);

xticks([0,3,6,9]); xlim([0,10]); xlabel('Time (Months)')
Rptext = sprintf(['R^2 = %.2f,' newline 'p = %.1d'],...
            lm.Rsquared.Adjusted,lm.ModelFitVsNullModel.Pvalue);
if ~use_s
    yplotmax = max(s_deltaRmat(:)); yplotmin = min(s_deltaRmat(:));
    ytickvec = [yplotmin,(yplotmax+yplotmin)/2,yplotmax];
    ylim([1.1*yplotmin,1.25*yplotmax]); yticks(ytickvec); 
    yticklabels({num2str(ytickvec(1),'%.2f'),num2str(ytickvec(2),'%.2f'),num2str(ytickvec(3),'%.2f')}); 
    loc = 'northeast';    
    text(0.75,0.1,Rptext,'FontSize',20,'FontName','Times','Units','normalized');
else
    yplotmax = 1.1; yplotmin = -0.1; ylim([yplotmin,yplotmax]); 
    yticks([0,0.5,1]); yticklabels({'0','0.5','1'}); 
    loc = 'southeast';
    text(0.75,0.95,Rptext,'FontSize',20,'FontName','Times','Units','normalized');
end
ylabel(ylab); title(titlestr);
legend(plothands,studylabels,'Location',loc,'NumColumns',3,'FontSize',20);
set(gca,'FontSize',24,'FontName','Times');
% 
% gplot = gscatter(tptvec,s_deltaRvec,gposvec,...
%     cmap,shapes,7,'doleg','off');

% xticks(xpos); xlim([min(xposvec)-0.5,max(xposvec)+0.5]); xticklabels(xlabs);
% xlabel([]);
% yplotmax = max(s_deltaRmat(:)); yplotmin = 0;
% ylim([yplotmin,1.1*yplotmax]); yticks([yplotmin,(yplotmin+yplotmax)/2,yplotmax]);
% yticklabels({'0',num2str((yplotmin+yplotmax)/2,'%.2f'),num2str(yplotmax,'%.2f')})
% ylabel('R'); title('Per Timepoint');
% legend(legstr,{'fit s','ret','ant','nd'},'Location','northeast','NumColumns',4,'FontSize',22,'box','off')
% set(gca,'FontSize',26,'FontName','Times');
end