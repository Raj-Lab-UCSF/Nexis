%% 1. Run Nexis:global on all tau models and filter out low-performing datasets
clear; clc; rng(0);
matdir = '/Users/justintorok/Documents/MATLAB/Nexis/raw_data_mouse';
load([matdir filesep 'OutputStruct_CTCorrsReg.mat']); % code below already run
% load([matdir filesep 'eNDM_mousedata.mat'],'data426');
% studylist = fieldnames(data426).';
% load([matdir filesep 'KaufmanDiamond_datasets_dat&seed.mat'],'data426');
% studylist = [studylist fieldnames(data426).'];
% R2s = zeros(1,length(studylist));
% outputs_R2studies = struct;
% rng(0);
% 
% for i = 1:length(studylist)
%     outputs_ndm = stdNDM_mouse('study',studylist{i},'bootstrapping',0,...
%         'w_dir',0,'volcorrect',1);
%     R2s(i) = outputs_ndm.ndm.Full.results.lm_Rsquared_adj;
%     outputs_R2studies.(studylist{i}) = outputs_ndm.ndm.Full;
% end
% thresh = 0.25; % R across all timepoints > 0.5
% datsetnames = studylist(R2s > thresh);
% for i = 1:length(studylist)
%     if ~ismember(studylist{i},datsetnames)
%         outputs_R2studies = rmfield(outputs_R2studies,studylist{i});
%     end
% end

%% 2. Create correlation plots with the Tasic, et al. dataset
plotting = 1;
savenclose = 0;
addpath('/Users/justintorok/Documents/MATLAB/Nexis/MATLAB/lib_eNDM_general')
datapath = '/Users/justintorok/Documents/MATLAB/Nexis/raw_data_mouse';
datapath2 = '/Users/justintorok/Documents/MATLAB/MISS/MISS-MatFiles';

% Load datasets of interest
load([datapath filesep 'classkey_tasic.mat']);
load([datapath filesep 'Tasic_CTMaps.mat']);
load([datapath2 filesep 'Tasic_AllGenes_Inputs.mat'],'genevct_allgenes')
load([datapath2 filesep 'Tasic_Inputs.mat'],'genevct')
load([datapath2 filesep 'MRx3_L90_inds.mat'],'geneinds');

% Create x-axis labels
subclasses = classkey_tasic;
strepper1 = @(x) strrep(x,'_','-');
subclasses = cellfun(strepper1,subclasses,'UniformOutput',0);
strepper2 = @(x) strrep(x,'Macro','Micro');
subclasses = cellfun(strepper2,subclasses,'UniformOutput',0);
strepper3 = @(x) strrep(x,'L2-3','L2/3');
subclasses = cellfun(strepper3,subclasses,'UniformOutput',0);
datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames,'UniformOutput',0);

% Calculate and store correlations
resultsstruct = struct;
corrmat_all = zeros(length(datsetnames),length(classkey_tasic));
parcorrmat_all = corrmat_all;
logpmat_all = corrmat_all; 
parlogpmat_all = corrmat_all;

for j = 1:length(datsetnames)
    datsetname = datsetnames{j};
    Y = outputs_R2studies.(datsetname).data;
    naninds = isnan(Y(:,1));
    Y_end = Y(~naninds,end);
    Xct = Tasic_ng606(~naninds,:);
    testmat = [Y_end Xct];
    [corrmat,corrp] = corrcoef(testmat);
    resultsstruct.(datsetname).corrmat = corrmat;
    resultsstruct.(datsetname).logpmat = -log10(corrp);
    corrmat_all(j,:) = corrmat(1,2:end);
    logpmat_all(j,:) = -log10(corrp(1,2:end));
    [pcorrmat,pcorrp] = partialcorr(testmat);
    resultsstruct.(datsetname).parcorrmat = pcorrmat;
    resultsstruct.(datsetname).parlogpmat = -log10(pcorrp);
    parcorrmat_all(j,:) = pcorrmat(1,2:end);
    parlogpmat_all(j,:) = -log10(pcorrp(1,2:end));
end

if plotting
    cmap = hsv(length(classkey_tasic)); 
    shapes = {'o','s','d','^','v','<','>','p','h','+','.','x'};
    g = 1:length(classkey_tasic);
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
    plot([0,length(classkey_tasic)+1],[0,0],'LineStyle','--','Color',[0.25 0.25 0.25],'LineWidth',1.5);
    xlim([0,length(classkey_tasic)+1])
    plothands_leg = plothands_leg(1,:);
    hleg = legend(plothands_leg,datsetlabels,'Location','northeast','NumColumns',4,'FontSize',20);
    set(gca,'XTick',1:length(subclasses),'XTickLabel',subclasses,'XTickLabelRotation',45);
    set(gca,'YTick',[-1 -0.5 0 0.5 1])
    ylabel("Pearson's R");
    xlabel('');
    set(gca,'TickLength',[0 0])
    title('Pearson Correlations with End Tau Pathology','FontSize',32);
    set(gca, 'FontSize', 28, 'FontName', 'Times');
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('CorrelationsBarChart_Tasic','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    f3 = figure('Position',[0 0 1500 400]); hold on;
    b = boxplot(parcorrmat_all,g,'Colors',cmap,'Symbol','');
    set(b,{'linew'},{2})
    plothands_leg = [];
    for i = 1:length(datsetnames)
        v = gscatter(xpos(i,:),parcorrmat_all(i,:),g,cmap,shapes{i},7,'off');
        if strcmp('Hurtado',datsetnames{i})
            scatter(xpos(i,:),parcorrmat_all(i,:),75,['k' shapes{i}],'filled')
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
    plot([0,length(classkey_tasic)+1],[0,0],'LineStyle','--','Color',[0.25 0.25 0.25],'LineWidth',1.5);    
    xlim([0,length(classkey_tasic)+1]);
    plothands_leg = plothands_leg(1,:);
    hleg = legend(plothands_leg,datsetlabels,'Location','northeast','NumColumns',4,'FontSize',20);
    set(gca,'XTick',1:length(subclasses),'XTickLabel',subclasses,'XTickLabelRotation',45);
    set(gca,'YTick',[-1 -0.5 0 0.5 1])
    ylabel("Partial Correlation \rho");
    xlabel('');
    set(gca,'TickLength',[0 0])
    title('Partial Correlations with End Tau Pathology','FontSize',32);
    set(gca, 'FontSize', 28, 'FontName', 'Times');
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('PartialCorrelationsBarChart_Tasic','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    figure('Units','inches','Position',[0 0 14 6])
    I = imagesc(logpmat_all,[0,10]); colormap hot; colorbar;
    set(gca,'TickLength',[0 0],'DataAspectRatio',[1 1 1],'XTick',1:25,...
        'XTickLabelRotation',90,'XTickLabel',subclasses,...
        'YTick',1:8,'YTickLabel',datsetlabels,'FontName','Times','FontSize',20)
    title("-log_1_0(p-value), Pearson's R")
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('PValuePlot_RawCorrs','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    figure('Units','inches','Position',[0 0 14 6])
    I = imagesc(parlogpmat_all,[0,10]); colormap hot; colorbar;
    set(gca,'TickLength',[0 0],'DataAspectRatio',[1 1 1],'XTick',1:25,...
        'XTickLabelRotation',90,'XTickLabel',subclasses,...
        'YTick',1:8,'YTickLabel',datsetlabels,'FontName','Times','FontSize',20)
    title("-log_1_0(p-value), Partial Correlation \rho")
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('PValuePlot_PartialCorrs','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end
end

sigtypes_partialp05 = cell(1,length(datsetlabels)); sigthreshp05 = -log10(0.05);
sigtypes_partialp01 = cell(1,length(datsetlabels)); sigthreshp01 = -log10(0.01);
sigtypes_partialp001 = cell(1,length(datsetlabels)); sigthreshp001 = -log10(0.001);
for i = 1:length(datsetnames)
    sigtypes_partialp05{i} = find(parlogpmat_all(i,:) > sigthreshp05);
    sigtypes_partialp01{i} = find(parlogpmat_all(i,:) > sigthreshp01);
    sigtypes_partialp001{i} = find(parlogpmat_all(i,:) > sigthreshp001);
end

%% 3. Cross-correlation matrices with hierarchical clustering
plotting = 1;
savenclose = 0;

% Hierarchical clustering of cell-type distributions for gene-expression
% cross-correlation (all genes)
Xct_corr_g = genevct_allgenes;
ct_dist_g = pdist(Xct_corr_g.','cityblock');
ct_linkage_g = linkage(ct_dist_g,'average');
ct_reorder_g = optimalleaforder(ct_linkage_g,ct_dist_g);
subclasses_crosscorr_g = subclasses(ct_reorder_g);
ct_corrs_g = corrcoef(Xct_corr_g);

% Hierarchical clustering of cell-type distributions for gene-expression
% cross-correlation (MRx3 genes)
Xct_corr_mrx3 = genevct(sort(geneinds(1:606)),:); % First 606 MRx3 genes, Mezias, et al. 2022
ct_dist_mrx3 = pdist(Xct_corr_mrx3.','cityblock');
ct_linkage_mrx3 = linkage(ct_dist_mrx3,'average');
ct_reorder_mrx3 = optimalleaforder(ct_linkage_mrx3,ct_dist_mrx3);
subclasses_crosscorr_mrx3 = subclasses(ct_reorder_mrx3);
ct_corrs_mrx3 = corrcoef(Xct_corr_mrx3);

% Hierarchical clustering of cell-type distributions for spatial cross-correlation
Xct_corr_s = Tasic_ng606;
ct_dist_s = pdist(Xct_corr_s.','cityblock');
ct_linkage_s = linkage(ct_dist_s,'average');
ct_reorder_s = optimalleaforder(ct_linkage_s,ct_dist_s);
subclasses_crosscorr_s = subclasses(ct_reorder_s);
ct_corrs_s = corrcoef(Xct_corr_s); 

if plotting
    % Dendrogram, on spatial data
    figure('Units','inches','Position',[0 0 10 2])
    D = dendrogram(ct_linkage_s,'Reorder',ct_reorder_s);
    for i = 1:length(D)
        D(i).Color = [0 0 0];
    end
    set(D,'LineWidth',1.5);
    set(gca,'visible','off');
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('Dendrogram_Tasic_Spatial','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    % Dendrogram, all-gene expression data
    figure('Units','inches','Position',[0 0 10 2])
    D = dendrogram(ct_linkage_g,'Reorder',ct_reorder_g);
    for i = 1:length(D)
        D(i).Color = [0 0 0];
    end
    set(D,'LineWidth',1.5);
    set(gca,'visible','off');
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('Dendrogram_Tasic_ByGene','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end   

    % Dendrogram, MRx3 gene expression data
    figure('Units','inches','Position',[0 0 10 2])
    D = dendrogram(ct_linkage_mrx3,'Reorder',ct_reorder_mrx3);
    for i = 1:length(D)
        D(i).Color = [0 0 0];
    end
    set(D,'LineWidth',1.5);
    set(gca,'visible','off');
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('Dendrogram_Tasic_ByMRx3','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    % Heatmap, spatial correlation, spatial clustering 
    figure('Units','inches','Position',[0 0 14 10])
    I = imagesc(ct_corrs_s(ct_reorder_s,ct_reorder_s),[-0.6,0.6]); colormap redblue;
    set(gca,'TickLength',[0 0],'DataAspectRatio',[1 1 1],'XTick',1:25,...
        'XTickLabelRotation',90,'XTickLabel',subclasses_crosscorr_s,...
        'YTickLabel',{''},'FontName','Times','FontSize',20)
    for i = 1:length(subclasses_crosscorr_s)
        text(25.75,i,subclasses_crosscorr_s{i},'FontSize',19,'FontName','Times')
    end
    cb = colorbar; a = get(cb); a = a.Position;
    set(cb,'Position',[a(1)+0.115,a(2),2*a(3),a(4)])
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('CrossCorr_Tasic_SpatialCorrs_SpatialClust','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    % Heatmap, spatial correlation, all-gene clustering 
    figure('Units','inches','Position',[0 0 14 10])
    I = imagesc(ct_corrs_s(ct_reorder_g,ct_reorder_g),[-0.6,0.6]); colormap redblue;
    set(gca,'TickLength',[0 0],'DataAspectRatio',[1 1 1],'XTick',1:25,...
        'XTickLabelRotation',90,'XTickLabel',subclasses_crosscorr_g,...
        'YTickLabel',{''},'FontName','Times','FontSize',20)
    for i = 1:length(subclasses_crosscorr_g)
        text(25.75,i,subclasses_crosscorr_g{i},'FontSize',19,'FontName','Times')
    end
    cb = colorbar; a = get(cb); a = a.Position;
    set(cb,'Position',[a(1)+0.115,a(2),2*a(3),a(4)])
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('CrossCorr_Tasic_SpatialCorrs_AllGeneClust','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    % Heatmap, spatial correlation, all-gene clustering 
    figure('Units','inches','Position',[0 0 14 10])
    I = imagesc(ct_corrs_s(ct_reorder_mrx3,ct_reorder_mrx3),[-0.6,0.6]); colormap redblue;
    set(gca,'TickLength',[0 0],'DataAspectRatio',[1 1 1],'XTick',1:25,...
        'XTickLabelRotation',90,'XTickLabel',subclasses_crosscorr_mrx3,...
        'YTickLabel',{''},'FontName','Times','FontSize',20)
    for i = 1:length(subclasses_crosscorr_mrx3)
        text(25.75,i,subclasses_crosscorr_mrx3{i},'FontSize',19,'FontName','Times')
    end
    cb = colorbar; a = get(cb); a = a.Position;
    set(cb,'Position',[a(1)+0.115,a(2),2*a(3),a(4)])
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('CrossCorr_Tasic_SpatialCorrs_MRx3Clust','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end
    
    % Heatmap, all-gene correlation, all-gene clustering
    cmap_corrs_g = [[ones(50,1), 0.9*linspace(1,0.5,50).', 0.9*linspace(1,0,50).'];...
        [ones(50,1), linspace(0.5,0,50).', 0*ones(50,1)]];
    corrlim = 0.75;
    figure('Units','inches','Position',[0 0 14 10])
    I = imagesc(ct_corrs_g(ct_reorder_g,ct_reorder_g),[corrlim,1]); colormap(cmap_corrs_g);
    set(gca,'TickLength',[0 0],'DataAspectRatio',[1 1 1],'XTick',1:25,...
        'XTickLabelRotation',90,'XTickLabel',subclasses_crosscorr_g,...
        'YTickLabel',{''},'FontName','Times','FontSize',20)
    for i = 1:length(subclasses_crosscorr_g)
        text(25.75,i,subclasses_crosscorr_g{i},'FontSize',19,'FontName','Times')
    end
    cb = colorbar; a = get(cb); a = a.Position;
    cbnumticks = int8(1 + (1-corrlim)/0.05);
    cbticks = linspace(corrlim,1,cbnumticks);
    cbticklabels = cell(1,cbnumticks);
    for j = 1:cbnumticks
        if j == 1 && corrlim ~= 0 
            if mod(100*cbticks(j),10) == 0
                cbticklabels{j} = sprintf('< %.1f',cbticks(j));
            else
                cbticklabels{j} = sprintf('< %.2f',cbticks(j));
            end
        else
            cbticklabels{j} = num2str(cbticks(j),2);
        end
    end
    set(cb,'Position',[a(1)+0.115,a(2),2*a(3),a(4)],'Ticks',cbticks,'TickLabels',...
        cbticklabels);
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('CrossCorr_Tasic_AllGeneCorrs_AllGeneClust','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    % Heatmap, MRx3 correlation, MRx3 clustering
    figure('Units','inches','Position',[0 0 14 10])
    I = imagesc(ct_corrs_mrx3(ct_reorder_mrx3,ct_reorder_mrx3),[-1,1]); colormap redblue;
    set(gca,'TickLength',[0 0],'DataAspectRatio',[1 1 1],'XTick',1:25,...
        'XTickLabelRotation',90,'XTickLabel',subclasses_crosscorr_mrx3,...
        'YTickLabel',{''},'FontName','Times','FontSize',20)
    for i = 1:length(subclasses_crosscorr_mrx3)
        text(25.75,i,subclasses_crosscorr_mrx3{i},'FontSize',19,'FontName','Times')
    end
    cb = colorbar; a = get(cb); a = a.Position;
    set(cb,'Position',[a(1)+0.115,a(2),2*a(3),a(4)]);
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('CrossCorr_Tasic_MRx3Corrs_MRx3Clust','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    % Heatmap, MRx3 correlation, spatial clustering
    figure('Units','inches','Position',[0 0 14 10])
    I = imagesc(ct_corrs_mrx3(ct_reorder_s,ct_reorder_s),[-1,1]); colormap redblue;
    set(gca,'TickLength',[0 0],'DataAspectRatio',[1 1 1],'XTick',1:25,...
        'XTickLabelRotation',90,'XTickLabel',subclasses_crosscorr_s,...
        'YTickLabel',{''},'FontName','Times','FontSize',20)
    for i = 1:length(subclasses_crosscorr_s)
        text(25.75,i,subclasses_crosscorr_s{i},'FontSize',19,'FontName','Times')
    end
    cb = colorbar; a = get(cb); a = a.Position;
    set(cb,'Position',[a(1)+0.115,a(2),2*a(3),a(4)]);
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('CrossCorr_Tasic_MRx3Corrs_SpatialClust','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end
end

%% 4. Linear models, top correlations
numsigtypes = 5;
Y_ends = cell(1,length(datsetnames));
sigtypes_corr_topn = cell(1,length(datsetnames));
mdls_corr_topn = cell(1,length(datsetnames));
Y_preds_corr_topn = mdls_corr_topn;
subclasses_corr_topn = mdls_corr_topn;

sigtypes_corr_bic = cell(1,length(datsetnames));
mdls_corr_bic = cell(1,length(datsetnames));
Y_preds_corr_bic = mdls_corr_bic;
subclasses_corr_bic = mdls_corr_bic;

sigtypes_parcorr_topn = cell(1,length(datsetnames));
mdls_parcorr_topn = cell(1,length(datsetnames));
Y_preds_parcorr_topn = mdls_parcorr_topn;
subclasses_parcorr_topn = mdls_parcorr_topn;

sigtypes_parcorr_bic = cell(1,length(datsetnames));
mdls_parcorr_bic = cell(1,length(datsetnames));
Y_preds_parcorr_bic = mdls_parcorr_bic;
subclasses_parcorr_bic = mdls_parcorr_bic;

plotting = 1;
savenclose = 0;
for i = 1:length(datsetnames)
    % Establish top N types by corr. and partial corr., also variable order
    [~,corrsortinds] = sort(abs(corrmat_all(i,:)),'descend');
    corrinds_topn = sort(corrsortinds(1:numsigtypes));
    sigtypes_corr_topn{i} = corrinds_topn;
    [~,parcorrsortinds] = sort(abs(parcorrmat_all(i,:)),'descend');
    parcorrinds_topn = sort(parcorrsortinds(1:numsigtypes));
    sigtypes_parcorr_topn{i} = parcorrinds_topn;
    
    % Fit models for top N types
    datsetname = datsetnames{i};
    Y = outputs_R2studies.(datsetname).data;
    naninds = isnan(Y(:,1));
    Y_end = Y(~naninds,end);
    Y_ends{i} = Y_end;
    Xct = Tasic_ng606(~naninds,:);
    Xct_sums = repmat(max(Xct),size(Xct,1),1);
    Xct_norm = Xct./Xct_sums;
    Xct_norm_corr_topn = Xct_norm(:,corrinds_topn);
    mdls_corr_topn{i} = fitlm(Xct_norm_corr_topn,Y_end);
    Y_preds_corr_topn{i} = [ones(length(Y_end),1), Xct_norm_corr_topn]...
        * mdls_corr_topn{i}.Coefficients.Estimate;
    subclasses_corr_topn{i} = subclasses(corrinds_topn);
    Xct_norm_parcorr_topn = Xct_norm(:,parcorrinds_topn);
    mdls_parcorr_topn{i} = fitlm(Xct_norm_parcorr_topn,Y_end);
    Y_preds_parcorr_topn{i} = [ones(length(Y_end),1), Xct_norm_parcorr_topn]...
        * mdls_parcorr_topn{i}.Coefficients.Estimate;
    subclasses_parcorr_topn{i} = subclasses(parcorrinds_topn);
    
    % Fit models based on BIC
    bic_calc = @(x,y,z) -2*x + (y * log(z));
    bics_corr = zeros(1,size(corrmat_all,2));
    bics_parcorr = zeros(1,size(corrmat_all,2));
    for j = 1:size(corrmat_all,2)
        corrinds_j = corrsortinds(1:j);
        Xct_norm_corr_j = Xct_norm(:,sort(corrinds_j));
        mdls_corr_j = fitlm(Xct_norm_corr_j,Y_end);
        bics_corr(j) = bic_calc(mdls_corr_j.LogLikelihood,(j+1),length(Y_end));
        parcorrinds_j = sort(parcorrsortinds(1:j));
        Xct_norm_parcorr_j = Xct_norm(:,parcorrinds_j);
        mdls_parcorr_j = fitlm(Xct_norm_parcorr_j,Y_end);
        bics_parcorr(j) = bic_calc(mdls_parcorr_j.LogLikelihood,(j+1),length(Y_end));
    end
    figure; scatter(1:25,bics_parcorr);
    [~, modelind_corr] = min(bics_corr);
    corrinds_bic = sort(corrsortinds(1:modelind_corr));
    sigtypes_corr_bic{i} = corrinds_bic;
    Xct_norm_corr_bic = Xct_norm(:,corrinds_bic);
    mdls_corr_bic{i} = fitlm(Xct_norm_corr_bic,Y_end);
    Y_preds_corr_bic{i} = [ones(length(Y_end),1), Xct_norm_corr_bic]...
        * mdls_corr_bic{i}.Coefficients.Estimate;
    subclasses_corr_bic{i} = subclasses(corrinds_bic);

    [~, modelind_parcorr] = min(bics_parcorr);
    parcorrinds_bic = sort(parcorrsortinds(1:modelind_parcorr));
    sigtypes_parcorr_bic{i} = parcorrinds_bic;
    Xct_norm_parcorr_bic = Xct_norm(:,parcorrinds_bic);
    mdls_parcorr_bic{i} = fitlm(Xct_norm_parcorr_bic,Y_end);
    Y_preds_parcorr_bic{i} = [ones(length(Y_end),1), Xct_norm_parcorr_bic]...
        * mdls_parcorr_bic{i}.Coefficients.Estimate;
    subclasses_parcorr_bic{i} = subclasses(parcorrinds_bic);
end

if plotting
    % Corr, top N
    f = figure('Name','CorrTopN','units','inch','position',[0 0 15 7]);
    cmap = hsv(length(datsetnames));
    for i = 1:length(datsetnames)
        subplot(2,4,i); % fix later for not 8 datasets
        scatter(Y_preds_corr_topn{i},Y_ends{i},50,cmap(i,:),shapes{i},'filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds_corr_topn{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds_corr_topn{i}); min_y = min(Y_ends{i});
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
        if round(min_x,2) ~= 0 
            xticklabels({num2str(min_x,'%.2f'),num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        else
            xticklabels({'0',num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        end
        if strcmp(datsetnames{i},'Clavaguera')
            yticklabels({num2str(min_y,'%.3f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
        else
            if round(min_y,2) ~= 0 
                yticklabels({num2str(min_y,'%.2f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            else
                yticklabels({'0',num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            end
        end
        % ax = gca; ax.YAxis.Exponent = -1; ax.XAxis.Exponent = -1;
        text(0.05*(max_x-min_x)+min_x, 0.925*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls_corr_topn{i}.Rsquared.Adjusted),...
            'FontName','Times','FontSize',16);
        title(datsetlabels{i},'FontSize',20);
        set(gca,'FontSize',16,'FontName','Times');
    end
    han=axes(f,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylh = ylabel(han,'Observed Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.75);
    xlh = xlabel(han,'Predicted Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.5);
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print(sprintf('LinearModelPlots_Corrs_Top%dTypes',numsigtypes),'-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end
    % Partial corr, top N
    f = figure('Name','ParCorrTopN','units','inch','position',[0 0 15 7]);
    cmap = hsv(length(datsetnames));
    for i = 1:length(datsetnames)
        subplot(2,4,i); % fix later for not 8 datasets
        scatter(Y_preds_parcorr_topn{i},Y_ends{i},50,cmap(i,:),shapes{i},'filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds_parcorr_topn{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds_parcorr_topn{i}); min_y = min(Y_ends{i});
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
        if round(min_x,2) ~= 0 
            xticklabels({num2str(min_x,'%.2f'),num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        else
            xticklabels({'0',num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        end
        if strcmp(datsetnames{i},'Clavaguera')
            yticklabels({num2str(min_y,'%.3f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
        else
            if round(min_y,2) ~= 0 
                yticklabels({num2str(min_y,'%.2f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            else
                yticklabels({'0',num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            end
        end
        % ax = gca; ax.YAxis.Exponent = -1; ax.XAxis.Exponent = -1;
        text(0.05*(max_x-min_x)+min_x, 0.925*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls_parcorr_topn{i}.Rsquared.Adjusted),...
            'FontName','Times','FontSize',16);
        title(datsetlabels{i},'FontSize',20);
        set(gca,'FontSize',16,'FontName','Times');
    end
    han=axes(f,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylh = ylabel(han,'Observed Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.75);
    xlh = xlabel(han,'Predicted Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.5);
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print(sprintf('LinearModelPlots_ParCorrs_Top%dTypes',numsigtypes),'-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    % Corr, BIC
    f = figure('Name','CorrBIC','units','inch','position',[0 0 15 7]);
    cmap = hsv(length(datsetnames));
    for i = 1:length(datsetnames)
        subplot(2,4,i); % fix later for not 8 datasets
        scatter(Y_preds_corr_bic{i},Y_ends{i},50,cmap(i,:),shapes{i},'filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds_corr_bic{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds_corr_bic{i}); min_y = min(Y_ends{i});
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
        if round(min_x,2) ~= 0 
            xticklabels({num2str(min_x,'%.2f'),num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        else
            xticklabels({'0',num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        end
        if strcmp(datsetnames{i},'Clavaguera')
            yticklabels({num2str(min_y,'%.3f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
        else
            if round(min_y,2) ~= 0 
                yticklabels({num2str(min_y,'%.2f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            else
                yticklabels({'0',num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            end
        end
        % ax = gca; ax.YAxis.Exponent = -1; ax.XAxis.Exponent = -1;
        text(0.05*(max_x-min_x)+min_x, 0.925*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls_corr_bic{i}.Rsquared.Adjusted),...
            'FontName','Times','FontSize',16);
        title(datsetlabels{i},'FontSize',20);
        set(gca,'FontSize',16,'FontName','Times');
    end
    han=axes(f,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylh = ylabel(han,'Observed Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.75);
    xlh = xlabel(han,'Predicted Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.5);
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('LinearModelPlots_Corrs_BIC','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end

    % Partial Corr, BIC
    f = figure('Name','ParCorrBIC','units','inch','position',[0 0 15 7]);
    cmap = hsv(length(datsetnames));
    for i = 1:length(datsetnames)
        subplot(2,4,i); % fix later for not 8 datasets
        scatter(Y_preds_parcorr_bic{i},Y_ends{i},50,cmap(i,:),shapes{i},'filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds_parcorr_bic{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds_parcorr_bic{i}); min_y = min(Y_ends{i});
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
        if round(min_x,2) ~= 0 
            xticklabels({num2str(min_x,'%.2f'),num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        else
            xticklabels({'0',num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        end
        if strcmp(datsetnames{i},'Clavaguera')
            yticklabels({num2str(min_y,'%.3f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
        else
            if round(min_y,2) ~= 0 
                yticklabels({num2str(min_y,'%.2f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            else
                yticklabels({'0',num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            end
        end
        % ax = gca; ax.YAxis.Exponent = -1; ax.XAxis.Exponent = -1;
        text(0.05*(max_x-min_x)+min_x, 0.925*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls_parcorr_bic{i}.Rsquared.Adjusted),...
            'FontName','Times','FontSize',16);
        title(datsetlabels{i},'FontSize',20);
        set(gca,'FontSize',16,'FontName','Times');
    end
    han=axes(f,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylh = ylabel(han,'Observed Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.75);
    xlh = xlabel(han,'Predicted Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.5);
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('LinearModelPlots_ParCorrs_BIC','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end    
end



%% 4. Ridge-based cell-type subset selection w.r.t. t(end) tau pathology 
rng(0);
lambdas = 2.^linspace(0,12,200);
lambdas_min = zeros(1,length(datsetnames));
lambdas_1se = lambdas_min;
sigtypeinds_min = cell(1,length(datsetnames));
sigtypeinds_1se = sigtypeinds_min; sigtypeinds_0 = sigtypeinds_min;
B_mins = sigtypeinds_min; B_1ses = sigtypeinds_min; B_0s = sigtypeinds_min;
kfold = 10;

plotting = 0;
savenclose = 0;
if plotting
    f1 = figure('Units','inches','Position',[0,0,20,10]);
end
for i = 1:length(datsetnames)
    % K-fold cross-validation for ridge regression
    datsetname = datsetnames{i};
    fprintf('Finding Optimal Lambda, %s \n',datsetname)
    Y = outputs_R2studies.(datsetname).data;
    naninds = isnan(Y(:,1));
    Y_end = Y(~naninds,end);
    Xct = Tasic_ng606(~naninds,:);
    Xct_sums = repmat(mean(Xct),size(Xct,1),1);
    Xct_norm = Xct./Xct_sums;
    c = cvpartition(length(Y_end),'KFold',kfold);
    residuals_sq = zeros(kfold,length(lambdas));
    for k = 1:kfold
        testinds = test(c,k);
        Y_k_test = Y_end(testinds); Xct_k_test = Xct_norm(testinds,:);
        traininds = ~testinds;
        Y_k_train = Y_end(traininds); Xct_k_train = Xct_norm(traininds,:);
        B = ridge(Y_k_train,Xct_k_train,lambdas,1);
        Y_k_test_mat = repmat(Y_k_test,1,length(lambdas));
        Y_k_test_mat_pred = Xct_k_test * B;
        Y_k_train_mat_pred = Xct_k_train * B;
        Y_diff_sq = (Y_k_test_mat - Y_k_test_mat_pred).^2;
        residuals_sq(k,:) = sum(Y_diff_sq);
    end
    % Determination of minimum and minimum + 1se lambda values
    mean_residuals_sq = mean(residuals_sq);
    std_residuals_sq = std(residuals_sq);
    [minmse,minind] = min(mean_residuals_sq);
    lambdas_min(i) = lambdas(minind);
    minmse_plus1se = minmse + std_residuals_sq(minind);
    mean_subtract_abs = abs(mean_residuals_sq - minmse_plus1se);
    mean_subtract_abs(lambdas < lambdas_min(i)) = Inf;
    [~,min1seind] = min(mean_subtract_abs);
    lambdas_1se(i) = lambdas(min1seind);
    
    % Bootstrapping for lambda = 0 (OLS), lambda_min, lambda_1se
    fprintf('Bootstrap-based Coefficient Estimation, %s \n',datsetname)
    niters = 1000;
    B_matrix_1se = zeros(size(Xct_norm,2),niters);
    B_matrix_min = B_matrix_1se;
    B_matrix_0 = B_matrix_1se;
    for j = 1:niters
        traininds = datasample(1:length(Y_end),length(Y_end),'Replace',true);
        Y_k_train = Y_end(traininds); Xct_k_train = Xct_norm(traininds,:);
        B_matrix_min(:,j) = ridge(Y_k_train,Xct_k_train,lambdas_min(i),1);
        B_matrix_1se(:,j) = ridge(Y_k_train,Xct_k_train,lambdas_1se(i),1);
        B_matrix_0(:,j) = ridge(Y_k_train,Xct_k_train,0,1);
    end

    sigthresh = 3; % Defines CI around mean to be +/- sigthresh * std
    B_mean_min = mean(B_matrix_min,2); B_std_min = std(B_matrix_min,[],2);
    B_upper_min = B_mean_min + sigthresh*B_std_min; B_lower_min = B_mean_min - sigthresh*B_std_min;
    sigtypeinds_min{i} = find(sign(B_upper_min) == sign(B_lower_min));
    B_mins{i} = B_matrix_min;

    B_mean_1se = mean(B_matrix_1se,2); B_std_1se = std(B_matrix_1se,[],2);
    B_upper_1se = B_mean_1se + sigthresh*B_std_1se; B_lower_1se = B_mean_1se - sigthresh*B_std_1se;
    sigtypeinds_1se{i} = find(sign(B_upper_1se) == sign(B_lower_1se));
    B_1ses{i} = B_matrix_1se;
    
    B_mean_0 = mean(B_matrix_0,2); B_std_0 = std(B_matrix_0,[],2);
    B_upper_0 = B_mean_0 + sigthresh*B_std_0; B_lower_0 = B_mean_0 - sigthresh*B_std_0;
    sigtypeinds_0{i} = find(sign(B_upper_0) == sign(B_lower_0));
    B_0s{i} = B_matrix_0;
    B_means_mat = [B_mean_0,B_mean_min,B_mean_1se];

    if plotting
        figure(f1)
        subplot(2,4,i); hold on;
        inds = 1:2:(length(lambdas)-1);
        errorbar(-log2(lambdas(inds)),mean_residuals_sq(inds),std_residuals_sq(inds),'o'); 
        plot([-log2(lambdas_min(i)), -log2(lambdas_min(i))],...
            [min(mean_residuals_sq)-1.1*max(std_residuals_sq),max(mean_residuals_sq)+1.1*max(std_residuals_sq)],...
            'LineStyle',':','LineWidth',2,'Color','b');
        plot([-log2(lambdas_1se(i)), -log2(lambdas_1se(i))],...
            [min(mean_residuals_sq)-1.1*max(std_residuals_sq),max(mean_residuals_sq)+1.1*max(std_residuals_sq)],...
            'LineStyle',':','LineWidth',2,'Color','g');
        if i == 1
            legend({sprintf('%d-fold CV',kfold),'Min MSE', 'Min MSE + 1SE'},'Location','southeast',...
                'FontSize',16);
        end
        ylim([min(mean_residuals_sq)-1.1*max(std_residuals_sq),max(mean_residuals_sq)+1.1*max(std_residuals_sq)])
        xlim([-log2(lambdas(end)),-log2(lambdas(1))])
        xlabel('-log2(\lambda)'); 
        if i == 1 || i == 5
            ylabel(sprintf('MSE (%d-fold CV)',kfold)); 
        end
        title(datsetlabels{i});
        set(gca,'FontSize',20,'FontName','Times')
    end
    cmap = hsv(length(subclasses));
    if plotting 
        figure('Position',[0 0 1000 1500]); 
        subplot(3,1,1); hold on;
        b = boxplot(B_matrix_0.',g,'Colors',cmap,'Symbol','');
        set(b,{'linew'},{1})
        plot([0,length(classkey_tasic)+1],[0,0],'k--','LineWidth',2);
        set(gca,'XTick',1:length(subclasses));
        set(gca,'XTickLabel',{})
        ylim_ = max(abs([max(B_matrix_0(:)),min(B_matrix_0(:))]));
        set(gca,'YTick',[-ylim_,ylim_])
        ytickformat('%0.2f')
        text(-0.45,0,'0','FontSize',16,'FontName','Times')
        ylim([-1.1*ylim_,1.1*ylim_]);
        ylabel('Regression Coefficients');
        xlabel('');
        set(gca,'TickLength',[0 0])
        title(['OLS Coefficients, ' datsetlabels{i}],'FontSize',28);
        set(gca, 'FontSize', 16, 'FontName', 'Times');
        subplot(3,1,2); hold on;
        b = boxplot(B_matrix_min.',g,'Colors',cmap,'Symbol','');
        set(b,{'linew'},{1})
        plot([0,length(classkey_tasic)+1],[0,0],'k--','LineWidth',2);
        set(gca,'XTick',1:length(subclasses));
        set(gca,'XTickLabel',{});
        ylim_ = max(abs([max(B_matrix_min(:)),min(B_matrix_min(:))]));
        set(gca,'YTick',[-ylim_,ylim_])
        ytickformat('%0.2f')
        text(-0.45,0,'0','FontSize',16,'FontName','Times')
        ylim([-1.1*ylim_,1.1*ylim_]);
        ylabel('Regression Coefficients');
        xlabel('');
        set(gca,'TickLength',[0 0])
        title(['\lambda_m_i_n Ridge Coefficients, ' datsetlabels{i}],'FontSize',28);
        set(gca, 'FontSize', 16, 'FontName', 'Times');
        subplot(3,1,3); hold on;
        b = boxplot(B_matrix_1se.',g,'Colors',cmap,'Symbol','');
        set(b,{'linew'},{1})
        plot([0,length(classkey_tasic)+1],[0,0],'k--','LineWidth',2);
        set(gca,'XTick',1:length(subclasses));
        set(gca,'XTickLabel',subclasses,'XTickLabelRotation',45);
        ylim_ = max(abs([max(B_matrix_1se(:)),min(B_matrix_1se(:))]));
        set(gca,'YTick',[-ylim_,ylim_])
        text(-0.45,0,'0','FontSize',16,'FontName','Times')
        ytickformat('%0.2f')
        ylim([-1.1*ylim_,1.1*ylim_]);
        ylabel('Regression Coefficients');
        xlabel('');
        set(gca,'TickLength',[0 0])
        title(['\lambda_1_s_e Ridge Coefficients, ' datsetlabels{i}],'FontSize',28);
        set(gca, 'FontSize', 16, 'FontName', 'Times');
        if savenclose
            cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
            print(sprintf('RidgeCoefBar_%s',datsetname),'-dtiffn','-r300'); close;
            cd /Users/justintorok/Documents/MATLAB/Nexis/;
        end
    end
end
if savenclose
    cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
    print(sprintf('Bootstrap_Lcurves_%d',kfold),'-dtiffn','-r300'); close;
    cd /Users/justintorok/Documents/MATLAB/Nexis/;
end

%% 5. Linear models for significant cell types
sigtypeinds_parcorr = sigtypes_partialp001;
mdls_parcorr_topn = cell(1,length(datsetnames));
Y_preds_parcorr_topn = mdls_parcorr_topn;
subclasses_parcorr_topn = mdls_parcorr_topn;
sigtypeinds_bootstrap = sigtypeinds_min; % OLS - can change, but OLS is sparsest and justifiable
mdls_bootstrap = cell(1,length(datsetnames));
Y_preds_bootstrap = mdls_bootstrap;
subclasses_bootstrap = mdls_bootstrap;
Y_ends = cell(1,length(datsetnames));

plotting = 1;
savenclose = 0;
for i = 1:length(datsetnames)
    datsetname = datsetnames{i};
    Y = outputs_R2studies.(datsetname).data;
    naninds = isnan(Y(:,1));
    Y_end = Y(~naninds,end);
    Xct = Tasic_ng606(~naninds,:);
    Xct_sums = repmat(max(Xct),size(Xct,1),1);
    Xct_norm = Xct./Xct_sums;
    Xct_norm_parcorr_topn = Xct_norm(:,sigtypeinds_parcorr{i});
    mdls_parcorr_topn{i} = fitlm(Xct_norm_parcorr_topn,Y_end);
    Y_preds_parcorr_topn{i} = [ones(length(Y_end),1), Xct_norm_parcorr_topn]...
        * mdls_parcorr_topn{i}.Coefficients.Estimate;
    subclasses_parcorr_topn{i} = subclasses(sigtypeinds_parcorr{i});
    Xct_norm_bootstrap = Xct_norm(:,sigtypeinds_bootstrap{i});
    mdls_bootstrap{i} = fitlm(Xct_norm_bootstrap,Y_end);
    Y_preds_bootstrap{i} = [ones(length(Y_end),1), Xct_norm_bootstrap]...
        * mdls_bootstrap{i}.Coefficients.Estimate;
    subclasses_bootstrap{i} = subclasses(sigtypeinds_bootstrap{i});
    Y_ends{i} = Y_end;
end

if plotting
    f = figure('units','inch','position',[0 0 15 7]);
    cmap = hsv(length(datsetnames));
    for i = 1:length(datsetnames)
        subplot(2,4,i); % fix later for not 8 datasets
        scatter(Y_preds_parcorr_topn{i},Y_ends{i},50,cmap(i,:),shapes{i},'filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds_parcorr_topn{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds_parcorr_topn{i}); min_y = min(Y_ends{i});
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
        if round(min_x,2) ~= 0 
            xticklabels({num2str(min_x,'%.2f'),num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        else
            xticklabels({'0',num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        end
        if strcmp(datsetnames{i},'Clavaguera')
            yticklabels({num2str(min_y,'%.3f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
        else
            if round(min_y,2) ~= 0 
                yticklabels({num2str(min_y,'%.2f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            else
                yticklabels({'0',num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            end
        end
        % ax = gca; ax.YAxis.Exponent = -1; ax.XAxis.Exponent = -1;
        text(0.05*(max_x-min_x)+min_x, 0.925*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls_parcorr_topn{i}.Rsquared.Adjusted),...
            'FontName','Times','FontSize',16);
        title(datsetlabels{i},'FontSize',20);
        set(gca,'FontSize',16,'FontName','Times');
    end
    han=axes(f,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylh = ylabel(han,'Observed Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.75);
    xlh = xlabel(han,'Predicted Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.5);
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('LinearModelPlots_PartialCorrs','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end
    
    f = figure('units','inch','position',[0 0 15 7]);
    cmap = hsv(length(datsetnames));
    for i = 1:length(datsetnames)
        subplot(2,4,i); % fix later for not 8 datasets
        scatter(Y_preds_bootstrap{i},Y_ends{i},50,cmap(i,:),shapes{i},'filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds_bootstrap{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds_bootstrap{i}); min_y = min(Y_ends{i});
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
        if round(min_x,2) ~= 0 
            xticklabels({num2str(min_x,'%.2f'),num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        else
            xticklabels({'0',num2str((max_x+min_x)/2,'%.2f'),...
                num2str(max_x,'%.2f')});
        end
        if strcmp(datsetnames{i},'Clavaguera')
            yticklabels({num2str(min_y,'%.3f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
        else
            if round(min_y,2) ~= 0 
                yticklabels({num2str(min_y,'%.2f'),num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            else
                yticklabels({'0',num2str((max_y+min_y)/2,'%.2f'),...
                    num2str(max_y,'%.2f')});
            end
        end
        % ax = gca; ax.YAxis.Exponent = -1; ax.XAxis.Exponent = -1;
        text(0.05*(max_x-min_x)+min_x, 0.925*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls_bootstrap{i}.Rsquared.Adjusted),...
            'FontName','Times','FontSize',16);
        title(datsetlabels{i},'FontSize',20);
        set(gca,'FontSize',16,'FontName','Times');
    end
    han=axes(f,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylh = ylabel(han,'Observed Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.75);
    xlh = xlabel(han,'Predicted Pathology at End Time Point','FontSize',20,...
        'FontName','Times','FontWeight','bold');
    xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.5);
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('LinearModelPlots_BootstrapOLS','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end
end

% %% 
% Y = outputs_R2studies.Hurtado.data;
% naninds = isnan(Y(:,1));
% Y_end = Y(~naninds,end);
% Xct = Tasic_ng606(~naninds,:);
% mdl = mvregress(Xct,Y_end);
% 
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
% plotting = 'false';
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
% plotting = 'false';
% if logical(plotting)
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
% 
% %% Stepwise regression w.r.t. t(end) tau pathology 
% 
% % mdls_stepwise = cell(1,length(datsetnames));
% % top_celltypes_stepwise = mdls_stepwise;
% % sortindsmat = zeros(length(datsetnames),length(classkey_tasic));
% % for i = 1:length(datsetnames)
% %     datsetname = datsetnames{i};
% %     Y = outputs_R2studies.(datsetname).data;
% %     naninds = isnan(Y(:,1));
% %     Y_end = Y(~naninds,end);
% %     Xct = Tasic_ng606(~naninds,:);
% %     Xct_sums = repmat(max(Xct),size(Xct,1),1);
% %     Xct_norm = Xct./Xct_sums;
% %     mdls_stepwise{i} = stepwiselm(Xct_norm,Y_end,'constant','Upper','linear','VarNames',...
% %         [subclasses 'Y_end'],'PEnter',0.01/length(subclasses),...
% %         'PRemove',0.05/length(subclasses));
% %     nonzeroinds = zeros(1,size(Xct,2));
% % end
% 
% %% Linear models for top three cell types for each study (Correlations)
% mdls = cell(1,length(datsetnames));
% Y_ends = mdls;
% Y_preds = mdls;
% subclasses_top3_corrs = mdls;
% for i = 1:length(datsetnames)
%     datsetname = datsetnames{i};
%     corrs_i = corrmat_all(i,:);
%     [~, inds_i] = sort(corrs_i,'descend');
%     subclasses_top3_corrs{i} = subclasses(inds_i(1:3));
%     Y = outputs_R2studies.(datsetname).data;
%     naninds = isnan(Y(:,1));
%     Y_end = Y(~naninds,end);
%     Xct = Tasic_ng606(~naninds,:);
%     Xct_sums = repmat(max(Xct),size(Xct,1),1);
%     Xct_norm = Xct./Xct_sums;
%     Xct_norm_top3 = Xct_norm(:,inds_i(1:3));
%     mdls{i} = fitlm(Xct_norm_top3,Y_end);
%     Y_preds{i} = [ones(length(Y_end),1), Xct_norm_top3] * mdls{i}.Coefficients.Estimate;
%     Y_ends{i} = Y_end;
% end
% 
% if logical(plotting)
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
%     print('Top3Scatters_ByCorr','-dpng','-r300'); close;
% end
% 
% %% Cell Type Cross-Correlation Matrix
