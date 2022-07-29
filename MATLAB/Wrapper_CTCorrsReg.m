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
        'w_dir',0,'volcorrect',1);
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

%% Create correlation plots with the Tasic, et al. dataset
addpath('/Users/justintorok/Documents/MATLAB/Nexis/MATLAB/lib_eNDM_general')
datapath = '/Users/justintorok/Documents/MATLAB/Nexis/raw_data_mouse';

% Set plot argument (still need to run this section for dependencies)
plotting = 0;
savenclose = 0;
rng(0)

% Load datasets of interest
load([datapath filesep 'Tasic_CTMaps.mat']);
load([datapath filesep 'classkey_tasic.mat']);

% Create x-axis labels
subclasses = classkey_tasic;
strepper1 = @(x) strrep(x,'_','-');
subclasses = cellfun(strepper1,subclasses,'UniformOutput',0);
strepper2 = @(x) strrep(x,'Macro','Micro');
subclasses = cellfun(strepper2,subclasses,'UniformOutput',0);
strepper3 = @(x) strrep(x,'L2-3','L2/3');
subclasses = cellfun(strepper3,subclasses,'UniformOutput',0);

% Calculate and store correlations
resultsstruct = struct;
corrmat_all = zeros(length(datsetnames),length(classkey_tasic));

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
end

% Hierarchical clustering of cell-type distributions for cross-correlation
Xct_corr = Tasic_ng606;
ct_dist = pdist(Xct_corr.','cosine');
ct_linkage = linkage(ct_dist,'average');
ct_reorder = optimalleaforder(ct_linkage,ct_dist);
subclasses_crosscorr = subclasses(ct_reorder);
ct_corrs = corrcoef(Xct(ct_reorder,ct_reorder));

plotting = 1;
savenclose = 0;
if plotting
    cmap = hsv(length(classkey_tasic)); %#ok<UNRCH> 
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
    plothands_leg = plothands_leg(1,:);
    datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames,'UniformOutput',0);
    hleg = legend(plothands_leg,datsetlabels,'Location','northeast','NumColumns',4,'FontSize',20);
    set(gca,'XTick',1:length(subclasses),'XTickLabel',subclasses,'XTickLabelRotation',45);
    set(gca,'YTick',[-1 -0.5 0 0.5 1])
    ylabel("Pearson's R");
    xlabel('');
    set(gca,'TickLength',[0 0])
%     title('Correlations with End Tau Pathology','FontSize',32);
    set(gca, 'FontSize', 28, 'FontName', 'Times');
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('CorrelationsBarChart_Tasic','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end
    
    figure('Units','inches','Position',[0 0 10 2])
    D = dendrogram(ct_linkage,'Reorder',ct_reorder);
    for i = 1:length(D)
        D(i).Color = [0 0 0];
    end
    set(D,'LineWidth',1.5);
    set(gca,'visible','off');
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('Dendrogram_Tasic','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end
   
    figure('Units','inches','Position',[0 0 14 10])
    I = imagesc(ct_corrs - eye(length(subclasses))); colormap redblue;
    set(gca,'TickLength',[0 0],'DataAspectRatio',[1 1 1],'XTick',1:25,...
        'XTickLabelRotation',90,'XTickLabel',subclasses_crosscorr,...
        'YTickLabel',{''},'FontName','Times','FontSize',20)
    for i = 1:length(subclasses_crosscorr)
        text(25.75,i,subclasses_crosscorr{i},'FontSize',19,'FontName','Times')
    end
    cb = colorbar; a = get(cb); a = a.Position;
    set(cb,'Position',[a(1)+0.115,a(2),2*a(3),a(4)])
    if savenclose
        cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
        print('CrossCorr_Tasic','-dtiffn','-r300'); close;
        cd /Users/justintorok/Documents/MATLAB/Nexis/;
    end
end

%% Ridge-based cell-type rankings w.r.t. t(end) tau pathology 
rng(0);
lambdas = 2.^linspace(-2,12,200);
lambdas_min = zeros(1,length(datsetnames));
lambdas_1se = lambdas_min;
datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames,'UniformOutput',0);
sigtypeinds_min = cell(1,length(datsetnames));
sigtypeinds_1se = sigtypeinds_min; sigtypeinds_0 = sigtypeinds_min;
B_mins = sigtypeinds_min; B_1ses = sigtypeinds_min; B_0s = sigtypeinds_min;
plotting = 1;
savenclose = 0;

for i = 1:length(datsetnames)
    datsetname = datsetnames{i};
    fprintf('Finding Optimal Lambda, %s \n',datsetname)
    Y = outputs_R2studies.(datsetname).data;
    naninds = isnan(Y(:,1));
    Y_end = Y(~naninds,end);
    Xct = Tasic_ng606(~naninds,:);
    Xct_sums = repmat(mean(Xct),size(Xct,1),1);
    Xct_norm = Xct./Xct_sums;
    kfold = 5;
    c = cvpartition(length(Y_end),'KFold',kfold);
    residuals_sq = zeros(kfold,length(lambdas));
    for k = 1:kfold
        testinds = test(c,k);
%         Y_k_test = Y_end(testinds); Xct_k_test = [ones(sum(testinds),1),Xct_norm(testinds,:)];
        Y_k_test = Y_end(testinds); Xct_k_test = Xct_norm(testinds,:);
        traininds = ~testinds;
        Y_k_train = Y_end(traininds); Xct_k_train = Xct_norm(traininds,:);
        B = ridge(Y_k_train,Xct_k_train,lambdas,1);
        Y_k_test_mat = repmat(Y_k_test,1,length(lambdas));
        Y_k_test_mat_pred = Xct_k_test * B;
        Y_k_train_mat_pred = Xct_k_train * B;
%         Y_k_train_mat_pred = [ones(size(Xct_k_train,1),1) Xct_k_train] * B;
        Y_diff_sq = (Y_k_test_mat - Y_k_test_mat_pred).^2;
        residuals_sq(k,:) = sum(Y_diff_sq);
%         corrs_lambda1(k) = corr(Y_k_test,Y_k_test_mat_pred(:,1));
    end
%     fprintf('Test R = %.2f \n',mean(corrs_lambda1))
    mean_residuals_sq = mean(residuals_sq);
    std_residuals_sq = std(residuals_sq);
    [minmse,minind] = min(mean_residuals_sq);
    lambdas_min(i) = lambdas(minind);
    minmse_plus1se = minmse + std_residuals_sq(minind);
    mean_subtract_abs = abs(mean_residuals_sq - minmse_plus1se);
    mean_subtract_abs(lambdas < lambdas_min(i)) = Inf;
    [~,min1seind] = min(mean_subtract_abs);
    lambdas_1se(i) = lambdas(min1seind);

    fprintf('Bootstrap-based Coefficient Estimation, %s \n',datsetname)
%     fprintf('Jackknife-based Coefficient Estimation, %s \n',datsetname)
%     c = cvpartition(length(Y_end),'LeaveOut');
    niters = 200;
    B_matrix_1se = zeros(size(Xct_norm,2),niters);
    B_matrix_min = B_matrix_1se;
    B_matrix_0 = B_matrix_1se;
%     B_matrix = zeros(size(Xct_norm,2)+1,length(Y_end));
    for j = 1:niters
%         testind = test(c,j);
%         Y_k_test = Y_end(testinds); Xct_k_test = [ones(sum(testinds),1),Xct_norm(testinds,:)];
%         traininds = ~testind;
        traininds = datasample(1:length(Y_end),length(Y_end),'Replace',true);
        Y_k_train = Y_end(traininds); Xct_k_train = Xct_norm(traininds,:);
        B_matrix_min(:,j) = ridge(Y_k_train,Xct_k_train,lambdas_min(i),1);
        B_matrix_1se(:,j) = ridge(Y_k_train,Xct_k_train,lambdas_1se(i),1);
        B_matrix_0(:,j) = ridge(Y_k_train,Xct_k_train,0,1);
%         corrs_lambda1(k) = corr(Y_k_test,Y_k_test_mat_pred(:,1));
    end
%     B_matrix = B_matrix(2:end,:);
    sigthresh = 3.35;
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
    
%     figure; bar(B_means_mat); title(datsetlabels{i});
%     cmap = hsv(length(classkey_tasic));
%     g = 1:length(classkey_tasic);
    if plotting
        figure('Position',[0 0 1000 1500]); 
        subplot(3,1,1); hold on;
        b = boxplot(B_matrix_0.',g,'Colors',cmap,'Symbol','');
        set(b,{'linew'},{1})
        plot([0,length(classkey_tasic)+1],[0,0],'k--','LineWidth',2);
        set(gca,'XTick',1:length(subclasses));
        set(gca,'XTickLabel',{})
        ylim_ = max(abs([max(B_matrix_0(:)),min(B_matrix_0(:))]));
        set(gca,'YTick',[-ylim_,0,ylim_])
        ytickformat('%.2d')
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
        set(gca,'YTick',[-ylim_,0,ylim_])
        ytickformat('%.2d')
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
        set(gca,'YTick',[-ylim_,0,ylim_])
        ytickformat('%.2d')
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

if plotting
    figure('Units','inches','Position',[0,0,20,10])
    for i = 1:length(datsetnames)
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
        xlabel('-log2(\lambda)'); 
        if i == 1 || i == 5
            ylabel(sprintf('MSE (%d-fold CV)',kfold)); 
        end
        title(datsetlabels{i});
        set(gca,'FontSize',20,'FontName','Times')
    end
end
% if plotting && savenclose
%     cd /Users/justintorok/Documents/MATLAB/Nexis/Figures/;
%     print(sprintf('Bootstrap_Lcurves_%d',kfold),'-dtiffn','-r300'); close;
%     cd /Users/justintorok/Documents/MATLAB/Nexis/;
% end
%% 
ranksmat = zeros(length(datsetnames),length(subclasses));
for i = 1:length(datsetnames)
    sortinds = sortindsmat(i,:);
    for j = 1:length(subclasses)
        ranksmat(i,j) = find(sortinds == j);
    end
end
ranksmat_mean = mean(ranksmat);
[~,ctranks] = sort(ranksmat_mean);
ranksmat_rankorder = ranksmat(:,ctranks);
subclasses_rankorder = subclasses(ctranks);

plotting = 'false';
if logical(plotting)
    figure('units','inch','position',[0 0 5 9]); 
    cmap = flipud(pink(25));
    xs = [1 length(datsetnames)]; ys = [1 length(subclasses)];
    imagesc(xs,ys,ranksmat_rankorder.'); 
    h = colorbar; set(h,'YDir','reverse');
    colormap(cmap);
    set(gca,'ytick',1:length(subclasses));
    figlabs = subclasses_rankorder;
    set(gca,'YTickLabel',figlabs);
    set(gca,'xtick',1:length(datsetnames));
    set(gca,'TickLength',[0 0]);
    set(gca,'XAxisLocation','top');
    set(gca,'XTickLabelRotation',45);
    datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames,'UniformOutput',0);
    set(gca,'XTickLabel',datsetlabels);
%     ylabel('Cell Types','FontSize',20);
    set(gca,'FontSize',16,'FontName','Times');
end
print('LASSORankings','-dpng','-r300'); close;

%% Linear models for top three cell types for each study (LASSO)
mdls = cell(1,length(datsetnames));
Y_ends = mdls;
Y_preds = mdls;
subclasses_top3 = mdls;
for i = 1:length(datsetnames)
    datsetname = datsetnames{i};
    B = Bs{i}; FitInfo = FitInfos{i};
    ind_top3 = find(FitInfo.DF==3,1,'last');
    lambda = FitInfo.Lambda(ind_top3);
    ct_top3 = top3_types{i};
    ctinds_top3 = zeros(1,length(ct_top3));
    for j = 1:length(ctinds_top3)
        ctinds_top3(j) = find(ismember(classkey_tasic,ct_top3{j}));
    end
    subclasses_top3{i} = subclasses(ctinds_top3);
    Y = outputs_R2studies.(datsetname).data;
    naninds = isnan(Y(:,1));
    Y_end = Y(~naninds,end);
    Xct = Tasic_ng606(~naninds,:);
    Xct_sums = repmat(max(Xct),size(Xct,1),1);
    Xct_norm = Xct./Xct_sums;
    Xct_norm_top3 = Xct_norm(:,ctinds_top3);
    mdls{i} = fitlm(Xct_norm_top3,Y_end);
    Y_preds{i} = [ones(length(Y_end),1), Xct_norm_top3] * mdls{i}.Coefficients.Estimate;
    Y_ends{i} = Y_end;
end
plotting = 'false';
if logical(plotting)
    f = figure('units','inch','position',[0 0 15 7]);
    cmap = hsv(length(datsetnames));
    datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames,'UniformOutput',0);
    for i = 1:length(datsetnames)
        subplot(2,4,i); % fix later for not 8 datasets
        scatter(Y_preds{i},Y_ends{i},50,cmap(i,:),shapes{i},'filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds{i}); min_y = min(Y_ends{i});
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
        % xtickformat('%.1d'); ytickformat('%.1d'); 
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
        subs_top3 = subclasses_top3{i};
        text(0.05*(max_x-min_x)+min_x, 0.925*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls{i}.Rsquared.Adjusted),...
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
    print('Top3Scatters_ByLasso','-dpng','-r300'); close;
end

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
mdls = cell(1,length(datsetnames));
Y_ends = mdls;
Y_preds = mdls;
subclasses_top3_corrs = mdls;
for i = 1:length(datsetnames)
    datsetname = datsetnames{i};
    corrs_i = corrmat_all(i,:);
    [~, inds_i] = sort(corrs_i,'descend');
    subclasses_top3_corrs{i} = subclasses(inds_i(1:3));
    Y = outputs_R2studies.(datsetname).data;
    naninds = isnan(Y(:,1));
    Y_end = Y(~naninds,end);
    Xct = Tasic_ng606(~naninds,:);
    Xct_sums = repmat(max(Xct),size(Xct,1),1);
    Xct_norm = Xct./Xct_sums;
    Xct_norm_top3 = Xct_norm(:,inds_i(1:3));
    mdls{i} = fitlm(Xct_norm_top3,Y_end);
    Y_preds{i} = [ones(length(Y_end),1), Xct_norm_top3] * mdls{i}.Coefficients.Estimate;
    Y_ends{i} = Y_end;
end

if logical(plotting)
    f = figure('units','inch','position',[0 0 15 7]);
    cmap = hsv(length(datsetnames));
    datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames,'UniformOutput',0);
    for i = 1:length(datsetnames)
        subplot(2,4,i); % fix later for not 8 datasets
        scatter(Y_preds{i},Y_ends{i},50,cmap(i,:),shapes{i},'filled');
        h = lsline; set(h,'LineWidth',2,'Color','k','LineStyle','--');
        max_x = max(Y_preds{i}); max_y = max(Y_ends{i});
        min_x = min(Y_preds{i}); min_y = min(Y_ends{i});
        xlim([min_x max_x]); ylim([min_y max_y*1.05]);
        xticks([min_x (max_x+min_x)/2 max_x]); yticks([min_y (max_y+min_y)/2 max_y]);
        % xtickformat('%.1d'); ytickformat('%.1d'); 
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
        subs_top3 = subclasses_top3{i};
        text(0.05*(max_x-min_x)+min_x, 0.925*(max_y-min_y)+min_y, sprintf('R^2 = %.2f',mdls{i}.Rsquared.Adjusted),...
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
    print('Top3Scatters_ByCorr','-dpng','-r300'); close;
end

%% Cell Type Cross-Correlation Matrix