function [corrmat_sv,corrmat_glob] = CorrelationBarPlot_NexIS(outstruct,...
    subclasses_,type,cmap_,figtype_,savenclose_,figdir_)
if nargin < 7
    figdir_ = cd;
    if nargin < 6
        savenclose_ = 0;
        if nargin < 5
            figtype_ = 'Yao';
            if nargin < 4
                cmap_ = 'hsv';
                if nargin < 3
                    type = 'DeltaR2';
                end
            end
        end
    end
end

datsetnames_ = fieldnames(outstruct);
if (length(subclasses_) == 42) && ~strcmp(figtype_,'Random')
    nonneuronal_inds = find(ismember(subclasses_,{'Endo','Astro','Oligo',...
        'Micro-PVM','SMC-Peri','VLMC'}));
    gaba_inds = find(ismember(subclasses_,{'Lamp5','Sncg','Meis2','CR',...
        'Pvalb','Sst','Sst Chodl','Vip'}));
    glutctx_other_inds = find(ismember(subclasses_,{'Car3','L4 RSP-ACA'}));
    ctxtest = @(x) strcmp(x(end),'X');
    glutctx_inds = find(logical(cell2mat(cellfun(ctxtest,subclasses_,'UniformOutput',false))));
    glutctx_inds = sort([glutctx_other_inds; glutctx_inds]);
    gluthipp_inds = setdiff((1:length(subclasses_)).',[nonneuronal_inds;gaba_inds;glutctx_inds]);
    indcell = {glutctx_inds,gluthipp_inds,gaba_inds,nonneuronal_inds};
    indtest = [glutctx_inds;gluthipp_inds;gaba_inds;nonneuronal_inds];
    if strcmp(cmap_,'cool')
        cmap_labs = cool(length(indcell));
    elseif strcmp(cmap_,'hsv')
        cmap_labs = hsv(length(indcell));
    elseif strcmp(cmap_,'lines')
        cmap_labs = lines(length(indcell));
    else
        cmap_labs = twocolor(cmap_(1,:),cmap_(2,:),length(indcell));
    end
elseif strcmp(figtype_,'Random')
    if strcmp(cmap_,'cool')
        cmap_labs = cool(length(datsetnames_));
    elseif strcmp(cmap_,'hsv')
        cmap_labs = hsv(length(datsetnames_));
    elseif strcmp(cmap_,'lines')
        cmap_labs = lines(length(datsetnames_));
    else
        cmap_labs = twocolor(cmap_(1,:),cmap_(2,:),length(datsetnames_));
    end
    indtest = 1:length(datsetnames_);
    indcell = cell(1,length(datsetnames_));
    for i = 1:length(indcell)   
        indcell{i} = i;
    end
end

subclasses_outstruct = fieldnames(outstruct.(datsetnames_{1}));
corrvec_glob = NaN(length(datsetnames_),1);
corrmat_sv = NaN(length(datsetnames_),(length(subclasses_outstruct)-1));
for i = 1:length(datsetnames_)
    datsetname_ = datsetnames_{i};
    outstruct_i = outstruct.(datsetname_);
    for j = 1:length(subclasses_outstruct)
        if strcmp(type,'DeltaR2')
            if j == 1
                corrvec_glob(i) = outstruct_i.(subclasses_outstruct{j}).nexis_global.Full.results.lm_Rsquared_adj;
            else
                corrmat_sv(i,(j-1)) = outstruct_i.(subclasses_outstruct{j}).nexis_sv.Full.results.lm_Rsquared_adj;
            end
        elseif strcmp(type,'DeltaR')
            if j == 1
                R2val_ij = outstruct_i.(subclasses_outstruct{j}).nexis_global.Full.results.lm_Rsquared_ord;
                corrvec_glob(i) = R2val_ij^(0.5);
            else
                R2val_ij = outstruct_i.(subclasses_outstruct{j}).nexis_sv.Full.results.lm_Rsquared_ord;
                corrmat_sv(i,(j-1)) = R2val_ij^(0.5);
            end
        end
    end
end

corrmat_glob = repmat(corrvec_glob,1,size(corrmat_sv,2));
corrmat_ = corrmat_sv - corrmat_glob;
if strcmp(figtype_,'Random')
    corrmat_ = corrmat_.';
    R95 = prctile(corrmat_,99,1);
end

% sigpvals = zeros(length(subclasses_),1);
% sigbonf = 0.05/length(subclasses_); % Bonferroni MHC
% for i = 1:length(subclasses_)
%     [~,p1] = ttest(corrmat_(:,i));
%     if p1 < sigbonf
%         sigpvals(i) = 1;
%     end
% end
% subclasses_sig = subclasses_(logical(sigpvals));
cmap_ = zeros(length(indtest),3);
for i = 1:length(indtest)
    colorinds = zeros(1,length(indcell));
    for j = 1:length(indcell)
        colorinds(j) = ismember(i,indcell{j});
    end
    cmap_(i,:) = cmap_labs(logical(colorinds),:);
end
cmap_labs = cmap_(indtest,:);
corrmat_ = corrmat_(:,indtest);
subclasses_ = subclasses_(indtest);
 
shapes = {'o','s','d','^','v','<','>','p','h','+','.','x'};
g = 1:size(corrmat_,2);
xpos = zeros(size(corrmat_));
xposscatter = @(y) 0.2 * (2*rand - 1) + y;
for i = 1:size(corrmat_,1)
    xpos(i,:) = xposscatter(g);
end
datsetlabels = cellfun(@(x)strrep(x,'_',' '),datsetnames_,'UniformOutput',0);

if length(subclasses_) < 10
    figure('Units','inches','Position',[0 0 15 5]); hold on;
else
    figure('Units','inches','Position',[0 0 25 8]); hold on;
end
b = boxplot(corrmat_,g,'Colors',cmap_labs,'Symbol','');
set(b,{'linew'},{2});
plothands_leg = [];
if ~strcmp(figtype_,'Random')
    for i = 1:size(corrmat_,1)
        shape_i = shapes{i};
        size_i = 7;   
        v = gscatter(xpos(i,:),corrmat_(i,:),g,cmap_labs,shape_i,size_i,'off');
        if strcmp('Hurtado',datsetnames_{i})
            scatter(xpos(i,:),corrmat_(i,:),75,['k' shape_i],'filled')
        end
        w = gscatter(NaN(1,length(g)),NaN(1,length(g)),g,'k',shape_i,size_i,'off');
        for n = 1:length(v)
            if ~strcmp('Hurtado',datsetnames_{i})
                set(v(n), 'MarkerFaceColor', cmap_labs(n,:));
            end
            set(w(n), 'MarkerFaceColor', 'k');
        end
        plothands_leg = [plothands_leg, w]; %#ok<AGROW> 
    end
else
    randinds = randperm(size(corrmat_,1));
    randinds = randinds(1:50);
    for i = randinds
        shape_i = 'o';
        size_i = 4;
        v = gscatter(xpos(i,:),corrmat_(i,:),g,cmap_labs,shape_i,size_i,'off');
        if ~strcmp(figtype_,'Random') && strcmp('Hurtado',datsetnames_{i})
            scatter(xpos(i,:),corrmat_(i,:),75,['k' shape_i],'filled')
        end
        w = gscatter(NaN(1,length(g)),NaN(1,length(g)),g,'k',shape_i,size_i,'off');
        for n = 1:length(v)
            if ~strcmp(figtype_,'Random') && ~strcmp('Hurtado',datsetnames_{i})
                set(v(n), 'MarkerFaceColor', cmap_labs(n,:));
            elseif strcmp(figtype_,'Random')
                set(v(n), 'MarkerFaceColor', cmap_labs(n,:));
            end
            set(w(n), 'MarkerFaceColor', 'k');
        end
        plothands_leg = [plothands_leg, w]; %#ok<AGROW> 
    end
end
plot([0,length(subclasses_)+1],[0,0],'LineStyle','--','Color',[0.25 0.25 0.25],'LineWidth',1.5);
xlim([0,length(subclasses_)+1])
plothands_leg = plothands_leg(1,:);
if ~strcmp(figtype_,'Random')
    plot([0,length(subclasses_)+1],[0,0],'LineStyle','--','Color',[0.25 0.25 0.25],'LineWidth',1.5);
    xlim([0,length(subclasses_)+1])
    plothands_leg = plothands_leg(1,:);
    legend(plothands_leg,datsetlabels,'Location','northeast','NumColumns',5,'FontSize',20,'box','off');
    xlabs = cell(1,length(subclasses_));
    for i = 1:length(subclasses_)
        col = cmap_labs(i,:);
    %     col = [0 0 0];
        lenfigtype = length(figtype_);
        if (lenfigtype < 5) || ~strcmp(figtype_((lenfigtype-4):lenfigtype),'Genes')
            xlabs{i} = sprintf('\\color[rgb]{%f,%f,%f}%s',col(1),col(2),col(3),subclasses_{i});
        else 
            subclasses_it = cellfun(@(x)['{\it ' x '}'],subclasses_,'UniformOutput',false);
            xlabs{i} = sprintf('\\color[rgb]{%f,%f,%f}%s',col(1),col(2),col(3),subclasses_it{i});
        end
    end
    if length(subclasses_) > 10
        set(gca,'XTick',1:length(subclasses_),'XTickLabel',xlabs,'TickLabelInterpreter','tex',...
             'XTickLabelRotation',90);
    else
        set(gca,'XTick',1:length(subclasses_),'XTickLabel',xlabs,'TickLabelInterpreter','tex');
    end
else
    plot([0,length(datsetnames_)+1],[0,0],'LineStyle','--','Color',[0.25 0.25 0.25],'LineWidth',1.5);
    xlim([0,length(datsetnames_)+1])
    xlabs = cell(1,length(datsetnames_));
    for i = 1:length(datsetnames_)
        h = plot([i-0.25,i+0.25],[R95(i),R95(i)],'Color','r','LineStyle',':','LineWidth',4);
        col = cmap_labs(i,:);
    %     col = [0 0 0];
        xlabs{i} = sprintf('\\color[rgb]{%f,%f,%f}%s',col(1),col(2),col(3),datsetlabels{i});
    end
    legend(h,'99^t^h Percentile','Location','northeast','FontSize',20,'box','off');
    if length(datsetlabels) > 10
        set(gca,'XTick',1:length(subclasses_),'XTickLabel',xlabs,'TickLabelInterpreter','tex',...
             'XTickLabelRotation',90);
    else
        set(gca,'XTick',1:length(subclasses_),'XTickLabel',xlabs,'TickLabelInterpreter','tex');
    end
        
end

set(gca,'YTick',[-1 -0.5 0 0.25 0.5 0.75])
if ismember('BoludaDSAD',datsetnames_) && strcmp(type,'DeltaR2')
    ylim([-0.1,0.75])
elseif ~ismember('BoludaDSAD',datsetnames_)
    ylim([-0.1,0.25])
else
    ylim([-0.1,0.5])
end

if strcmp(type,'DeltaR')
    ylabel("\DeltaR");
    figstr = [figdir_ filesep 'DeltaRBarChart_NexIS_'];
elseif strcmp(type,'DeltaR2')
    ylabel("\DeltaR^2");
    figstr = [figdir_ filesep 'DeltaR2BarChart_NexIS_'];
end

if length(subclasses_) == 42
    figstr = [figstr 'Yao'];
else
    figstr = [figstr figtype_];
end
xlabel('');
set(gca,'TickLength',[0 0])
set(gca, 'FontSize', 24, 'FontName', 'Times');
if savenclose_
    print(figstr,'-dtiffn','-r600'); close;
end
end