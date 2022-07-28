%% Run Nexis:global on all tau models and filter out low-performing datasets
clear; clc;
matdir = '/Users/justintorok/Documents/MATLAB/Nexis/raw_data_mouse';
load([matdir filesep 'eNDM_mousedata.mat'],'data426');
studylist = fieldnames(data426).';
load([matdir filesep 'KaufmanDiamond_datasets_dat&seed.mat'],'data426');
studylist = [studylist fieldnames(data426).'];
R2s = zeros(1,length(studylist));
outputs_R2studies = struct;
rng(0);

for i = 1:length(studylist)
    outputs_ndm = stdNDM_mouse('study',studylist{i},'bootstrapping',0,...
        'w_dir',1,'volcorrect',1);
    R2s(i) = outputs_ndm.ndm.Full.results.lm_Rsquared_adj;
    outputs_R2studies.(studylist{i}) = outputs_ndm.ndm.Full;
end
thresh = 0.25; % R across all timepoints > 0.5
datsetnames = studylist(R2s > thresh);
for i = 1:length(studylist)
    if ~ismember(studylist{i},datsetnames)
        outputs_R2studies = rmfield(outputs_R2studies,studylist{i});
    end
end

%% Create correlation plot with the Tasic, et al. dataset
addpath('/Users/justintorok/Documents/MATLAB/Nexis/MATLAB/lib_eNDM_general')
datapath = '/Users/justintorok/Documents/MATLAB/Nexis/raw_data_mouse';

% Load datasets of interest
load([datapath filesep 'Tasic_CTMaps.mat']);
load([datapath filesep 'classkey_tasic.mat']);

% Create x-axis labels
subclasses = ['NDM', classkey_tasic];
strepper1 = @(x) strrep(x,'_','-');
subclasses = cellfun(strepper1,subclasses,'UniformOutput',0);
strepper2 = @(x) strrep(x,'Macro','Micro');
subclasses = cellfun(strepper2,subclasses,'UniformOutput',0);
strepper3 = @(x) strrep(x,'L2-3','L2/3');
subclasses = cellfun(strepper3,subclasses,'UniformOutput',0);

% Calculate and store correlations
resultsstruct = struct;
corrmat_all = zeros(length(datsetnames),length(classkey_tasic)+1);

for j = 1:length(datsetnames)
    datsetname = datsetnames{j};
    Y = outputs_R2studies.(datsetname).data;
    Y2 = outputs_R2studies.(datsetname).predicted;
    naninds = isnan(Y(:,1));
    Y_end = Y(~naninds,end);
    Y2_end = Y2(~naninds,end);
    Xct = Tasic_ng606(~naninds,:);
    testmat = [Y_end Y2_end Xct];
    [corrmat,corrp] = corrcoef(testmat);
    resultsstruct.(datsetname).corrmat = corrmat;
    resultsstruct.(datsetname).logpmat = -log10(corrp);
    corrmat_all(j,:) = corrmat(1,2:end);
end

plotting = 1;
if logical(plotting)
    rng('default');
    cmap = hsv(length(classkey_tasic)+1);
    shapes = {'o','s','d','^','v','<','>','p','h','+','.','x'};
    g = 1:(length(classkey_tasic)+1);
    xpos = zeros(length(datsetnames),length(g));
    xposscatter = @(y) 0.2 * (2*rand - 1) + y;
    for i = 1:length(datsetnames)
        xpos(i,:) = xposscatter(g);
    end
    f3 = figure('Position',[0 0 1500 400]); hold on;
    b = boxplot(corrmat_all,g,'Colors',cmap,'Symbol','');
    set(b,{'linew'},{2})
    plothands_leg = [];
    for i = 1:length(datsetnames)
        v = gscatter(xpos(i,:),corrmat_all(i,:),g,cmap,shapes{i},7,'off');
        if strcmp('Hurtado',datsetnames{i})
            scatter(xpos(i,:),corrmat_all(i,:),75,['k' shapes{i}],'filled')
        end
        w = gscatter(NaN(1,length(g)),NaN(1,length(g)),g,'k',shapes{i},7,'off');
        for n = 1:length(v)
            if ~strcmp('Hurtado',datsetnames{i})
                set(v(n), 'MarkerFaceColor', cmap(n,:));
            end
            set(w(n), 'MarkerFaceColor', 'k');
        end
        plothands_leg = [plothands_leg, w];
    end
    plothands_leg = plothands_leg(1,:);
%     plot([0 length(classkey_tasic)+1],[median(corrmat_all) median(corrmat_all)],'k--','LineWidth',2);
    datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames,'UniformOutput',0);
    hleg = legend(plothands_leg,datsetlabels,'Location','northeast','NumColumns',4);
    set(gca,'XTick',1:length(subclasses),'XTickLabel',subclasses,'XTickLabelRotation',45);
    set(gca,'YTick',[-1 -0.5 0 0.5 1])
    ylabel("Pearson's R");
    xlabel('');
    set(gca,'TickLength',[0 0])
    title('Correlations with End Tau Pathology','FontSize',28);
    set(gca, 'FontSize', 24, 'FontName', 'Times');
end
print('AAIC2022_CorrelationsBarChart','-dpng','-r300'); close;

%% LASSO-based cell-type rankings w.r.t. t(end) tau pathology 
% 
% Bs = cell(1,length(datsetnames));
% FitInfos = Bs;
% top3_types = Bs;
% sortindsmat = zeros(length(datsetnames),length(classkey_tasic));
% for i = 1:length(datsetnames)
%     datsetname = datsetnames{i};
%     Y = outputs_R2studies.(datsetname).data;
%     naninds = isnan(Y(:,1));
%     Y_end = Y(~naninds,end);
%     Xct = Tasic_ng606(~naninds,:);
%     Xct_sums = repmat(max(Xct),size(Xct,1),1);
%     Xct_norm = Xct./Xct_sums;
%     [B, FitInfo] = lasso(Xct_norm,Y_end,'CV',10);
%     nonzeroinds = zeros(1,size(Xct,2));
%     for j = 1:size(Xct,2)
%         nonzeroinds(j) = find(B(j,:),1,'last');
%     end
%     [~,sortinds] = sort(nonzeroinds,'descend');
%     sortindsmat(i,:) = sortinds;
%     top3_types{i} = classkey_tasic(sortinds(1:3));
%     Bs{i} = B; FitInfos{i} = FitInfo;
% end
% ranksmat = zeros(length(datsetnames),length(subclasses));
% for i = 1:length(datsetnames)
%     sortinds = sortindsmat(i,:);
%     for j = 1:length(subclasses)
%         ranksmat(i,j) = find(sortinds == j);
%     end
% end
% ranksmat_mean = mean(ranksmat);
% [~,ctranks] = sort(ranksmat_mean);
% ranksmat_rankorder = ranksmat(:,ctranks);
% subclasses_rankorder = subclasses(ctranks);
% 
% plotting = 1;
% if logical(plotting)
%     figure('units','inch','position',[0 0 5 9]); 
%     cmap = flipud(pink(25));
%     xs = [1 length(datsetnames)]; ys = [1 length(subclasses)];
%     imagesc(xs,ys,ranksmat_rankorder.'); 
%     h = colorbar; set(h,'YDir','reverse');
%     colormap(cmap);
%     set(gca,'ytick',1:length(subclasses));
%     figlabs = subclasses_rankorder;
%     set(gca,'YTickLabel',figlabs);
%     set(gca,'xtick',1:length(datsetnames));
%     set(gca,'TickLength',[0 0]);
%     set(gca,'XAxisLocation','top');
%     set(gca,'XTickLabelRotation',45);
%     datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames,'UniformOutput',0);
%     set(gca,'XTickLabel',datsetlabels);
% %     ylabel('Cell Types','FontSize',20);
%     set(gca,'FontSize',16,'FontName','Times');
% end
% print('LASSORankings','-dpng','-r300'); close;
% 
% %% Linear models for top three cell types for each study (LASSO)
% mdls = cell(1,length(datsetnames));
% Y_ends = mdls;
% Y_preds = mdls;
% subclasses_top3 = mdls;
% for i = 1:length(datsetnames)
%     datsetname = datsetnames{i};
%     B = Bs{i}; FitInfo = FitInfos{i};
%     ind_top3 = find(FitInfo.DF==3,1,'last');
%     lambda = FitInfo.Lambda(ind_top3);
%     ct_top3 = top3_types{i};
%     ctinds_top3 = zeros(1,length(ct_top3));
%     for j = 1:length(ctinds_top3)
%         ctinds_top3(j) = find(ismember(classkey_tasic,ct_top3{j}));
%     end
%     subclasses_top3{i} = subclasses(ctinds_top3);
%     Y = outputs_R2studies.(datsetname).data;
%     naninds = isnan(Y(:,1));
%     Y_end = Y(~naninds,end);
%     Xct = Tasic_ng606(~naninds,:);
%     Xct_sums = repmat(max(Xct),size(Xct,1),1);
%     Xct_norm = Xct./Xct_sums;
%     Xct_norm_top3 = Xct_norm(:,ctinds_top3);
%     mdls{i} = fitlm(Xct_norm_top3,Y_end);
%     Y_preds{i} = [ones(length(Y_end),1), Xct_norm_top3] * mdls{i}.Coefficients.Estimate;
%     Y_ends{i} = Y_end;
% end
% 
% if plotting
%     f = figure('units','inch','position',[0 0 15 7]);
%     cmap = hsv(length(datsetnames));
%     datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames,'UniformOutput',0);
%     for i = 1:length(datsetnames)
%         subplot(2,4,i); % fix later for not 8 datasets
%         scatter(Y_preds{i},Y_ends{i},50,cmap(i,:),shapes{i},'filled');
%         h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
%         max_x = max(Y_preds{i}); max_y = max(Y_ends{i});
%         min_x = min(Y_preds{i}); min_y = min(Y_ends{i});
%         xlim([min_x max_x]); ylim([min_y max_y*1.05]);
%         xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
%         % xtickformat('%.1d'); ytickformat('%.1d'); 
%         if round(min_x,2) ~= 0 
%             xticklabels({num2str(min_x,'%.2f'),num2str((max_x+min_x)/2,'%.2f'),...
%                 num2str(max_x,'%.2f')});
%         else
%             xticklabels({'0',num2str((max_x+min_x)/2,'%.2f'),...
%                 num2str(max_x,'%.2f')});
%         end
%         if strcmp(datsetnames{i},'Clavaguera')
%             yticklabels({num2str(min_y,'%.3f'),num2str((max_y+min_y)/2,'%.2f'),...
%                     num2str(max_y,'%.2f')});
%         else
%             if round(min_y,2) ~= 0 
%                 yticklabels({num2str(min_y,'%.2f'),num2str((max_y+min_y)/2,'%.2f'),...
%                     num2str(max_y,'%.2f')});
%             else
%                 yticklabels({'0',num2str((max_y+min_y)/2,'%.2f'),...
%                     num2str(max_y,'%.2f')});
%             end
%         end
%         % ax = gca; ax.YAxis.Exponent = -1; ax.XAxis.Exponent = -1;
%         subs_top3 = subclasses_top3{i};
%         text(0.05*(max_x-min_x)+min_x, 0.925*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls{i}.Rsquared.Adjusted),...
%             'FontName','Times','FontSize',16);
%         title(datsetlabels{i},'FontSize',20);
%         set(gca,'FontSize',16,'FontName','Times');
%     end
%     han=axes(f,'visible','off'); 
%     han.XLabel.Visible='on';
%     han.YLabel.Visible='on';
%     ylh = ylabel(han,'Observed Pathology at End Time Point','FontSize',20,...
%         'FontName','Times','FontWeight','bold');
%     ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.75);
%     xlh = xlabel(han,'Predicted Pathology at End Time Point','FontSize',20,...
%         'FontName','Times','FontWeight','bold');
%     xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.5);
%     print('Top3Scatters_ByLasso','-dpng','-r300'); close;
% end

%% Stepwise regression w.r.t. t(end) tau pathology 

% mdls_stepwise = cell(1,length(datsetnames));
% top_celltypes_stepwise = mdls_stepwise;
% sortindsmat = zeros(length(datsetnames),length(classkey_tasic));
% for i = 1:length(datsetnames)
%     datsetname = datsetnames{i};
%     Y = outputs_R2studies.(datsetname).data;
%     naninds = isnan(Y(:,1));
%     Y_end = Y(~naninds,end);
%     Xct = Tasic_ng606(~naninds,:);
%     Xct_sums = repmat(max(Xct),size(Xct,1),1);
%     Xct_norm = Xct./Xct_sums;
%     mdls_stepwise{i} = stepwiselm(Xct_norm,Y_end,'constant','Upper','linear','VarNames',...
%         [subclasses 'Y_end'],'PEnter',0.01/length(subclasses),...
%         'PRemove',0.05/length(subclasses));
%     nonzeroinds = zeros(1,size(Xct,2));
% end

%% Linear models for top three cell types for each study (Correlations)
mdls_cts = cell(1,length(datsetnames));
mdls_ndm = mdls_cts;
mdls_comb = mdls_cts;
Y_ends = mdls_cts;
Y_preds_cts = mdls_cts;
Y_preds_ndm = mdls_cts;
Y_preds_comb = mdls_cts;
subclasses_top3_corrs = mdls_cts;
subclasses_cts = subclasses(2:end);
for i = 1:length(datsetnames)
    datsetname = datsetnames{i};
    corrs_i = abs(corrmat_all(i,2:end));
    [~, inds_i] = sort(corrs_i,'descend');
    subclasses_top3_corrs{i} = subclasses_cts(inds_i(1:3));
    Y = outputs_R2studies.(datsetname).data;
    Y2 = outputs_R2studies.(datsetname).predicted;
    naninds = isnan(Y(:,1));
    Y_end = Y(~naninds,end);
    Y2_end = Y2(~naninds,end);
    Xct = Tasic_ng606(~naninds,:);
    Xct_sums = repmat(max(Xct),size(Xct,1),1);
    Xct_norm = Xct./Xct_sums;
    Xct_norm_top3 = Xct_norm(:,inds_i(1:3));
    mdls_cts{i} = fitlm(Xct_norm_top3,Y_end);
    Y_preds_cts{i} = [ones(length(Y_end),1), Xct_norm_top3] * mdls_cts{i}.Coefficients.Estimate;
    mdls_ndm{i} = fitlm(Y2_end,Y_end);
    Y_preds_ndm{i} = [ones(length(Y_end),1), Y2_end] * mdls_ndm{i}.Coefficients.Estimate;
    mdls_comb{i} = fitlm([Y2_end, Xct_norm_top3],Y_end);
    Y_preds_comb{i} = [ones(length(Y_end),1), Y2_end, Xct_norm_top3] * mdls_comb{i}.Coefficients.Estimate;
    Y_ends{i} = Y_end;
end

if plotting
    cmap = hsv(length(datsetnames));
    datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames,'UniformOutput',0);
    for i = 1:length(datsetnames)
        f = figure('units','inch','position',[0 0 12 2.5]);
        subplot(1,3,1);
        scatter(Y_preds_ndm{i},Y_ends{i},50,cmap(i,:),'d','filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds_ndm{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds_ndm{i}); min_y = min(Y_ends{i});
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
        if round(min_x,3) ~= 0 
            xticklabels({num2str(min_x,'%.3f'),num2str((max_x+min_x)/2,'%.3f'),...
                num2str(max_x,'%.3f')});
        else
            xticklabels({'0',num2str((max_x+min_x)/2,'%.3f'),...
                num2str(max_x,'%.3f')});
        end

        if round(min_y,3) ~= 0 
            yticklabels({num2str(min_y,'%.3f'),num2str((max_y+min_y)/2,'%.3f'),...
                num2str(max_y,'%.3f')});
        else
            yticklabels({'0',num2str((max_y+min_y)/2,'%.3f'),...
                num2str(max_y,'%.3f')});
        end
        % ax = gca; ax.YAxis.Exponent = -1; ax.XAxis.Exponent = -1;
%         subs_top3 = subclasses_top3{i};
        text(0.5*(max_x-min_x)+min_x, 0.25*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls_ndm{i}.Rsquared.Adjusted),...
            'FontName','Times','FontSize',20,'EdgeColor','k');
        ylabel([datsetlabels{i} ' Tau'],'FontWeight','bold','FontSize',24)
%         if i == 1
%             title('NDM','FontSize',20);
%         end
        set(gca,'FontSize',20,'FontName','Times');

        subplot(1,3,2);
        scatter(Y_preds_cts{i},Y_ends{i},50,cmap(i,:),'o','filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds_cts{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds_cts{i}); min_y = min(Y_ends{i});
        % xtickformat('%.1d'); ytickformat('%.1d'); 
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
        if round(min_x,3) ~= 0 
            xticklabels({num2str(min_x,'%.3f'),num2str((max_x+min_x)/2,'%.3f'),...
                num2str(max_x,'%.3f')});
        else
            xticklabels({'0',num2str((max_x+min_x)/2,'%.3f'),...
                num2str(max_x,'%.3f')});
        end
        yticklabels([]);
%         if round(min_y,3) ~= 0 
%             yticklabels({num2str(min_y,'%.3f'),num2str((max_y+min_y)/2,'%.3f'),...
%                 num2str(max_y,'%.3f')});
%         else
%             yticklabels({'0',num2str((max_y+min_y)/2,'%.3f'),...
%                 num2str(max_y,'%.3f')});
%         end
     
        % ax = gca; ax.YAxis.Exponent = -1; ax.XAxis.Exponent = -1;
%         subs_top3 = subclasses_top3{i};
        text(0.5*(max_x-min_x)+min_x, 0.25*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls_cts{i}.Rsquared.Adjusted),...
            'FontName','Times','FontSize',20,'EdgeColor','k');
%         if i == 1
%             title('Cell Types','FontSize',20);
%         end
        set(gca,'FontSize',20,'FontName','Times');
        
        subplot(1,3,3);
        scatter(Y_preds_comb{i},Y_ends{i},50,cmap(i,:),'s','filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds_comb{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds_comb{i}); min_y = min(Y_ends{i});
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);

        if round(min_x,3) ~= 0 
            xticklabels({num2str(min_x,'%.3f'),num2str((max_x+min_x)/2,'%.3f'),...
                num2str(max_x,'%.3f')});
        else
            xticklabels({'0',num2str((max_x+min_x)/2,'%.3f'),...
                num2str(max_x,'%.3f')});
        end
        yticklabels([]);
%         if round(min_y,3) ~= 0 
%             yticklabels({num2str(min_y,'%.3f'),num2str((max_y+min_y)/2,'%.3f'),...
%                 num2str(max_y,'%.3f')});
%         else
%             yticklabels({'0',num2str((max_y+min_y)/2,'%.3f'),...
%                 num2str(max_y,'%.3f')});
%         end

        % ax = gca; ax.YAxis.Exponent = -1; ax.XAxis.Exponent = -1;
%         subs_top3 = subclasses_top3{i};
        text(0.5*(max_x-min_x)+min_x, 0.25*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls_comb{i}.Rsquared.Adjusted),...
            'FontName','Times','FontSize',20,'EdgeColor','k');
        set(gca,'FontSize',20,'FontName','Times');
        print(['AAIC2022_Scatters_' datsetnames{i}],'-dpng','-r300'); close;
    end
%     han=axes(f,'visible','off'); 
%     han.XLabel.Visible='on';
%     han.YLabel.Visible='on';
%     ylh = ylabel(han,'Observed Pathology at End Time Point','FontSize',20,...
%         'FontName','Times','FontWeight','bold');
%     ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.75);
%     xlh = xlabel(han,'Predicted Pathology at End Time Point','FontSize',20,...
%         'FontName','Times','FontWeight','bold');
%     xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.5);
%     print('Top3Scatters_ByCorr','-dpng','-r300'); close;
end

%% Cell Type Cross-Correlation Matrix
