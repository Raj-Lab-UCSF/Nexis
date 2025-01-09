function [corrmat_sv,corrmat_glob] = CorrelationBarPlot_NexIS_Random(outstruct,type,cmap_,figtype_,savenclose_,figdir_)
if nargin < 6
    figdir_ = cd;
    if nargin < 5
        savenclose_ = 0;
        if nargin < 4
            figtype_ = 'Genes';
            if nargin < 3
                cmap_ = 'hsv';
                if nargin < 2
                    type = 'DeltaR2';
                end
            end
        end
    end
end

datsetnames_ = fieldnames(outstruct);
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

corrmat_glob = repmat(corrvec_glob,1,length(subclasses_outstruct));
corrmat_ = corrmat_sv - corrmat_glob;

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
xpos = zeros(length(datsetnames_),length(g));
xposscatter = @(y) 0.2 * (2*rand - 1) + y;
for i = 1:length(datsetnames_)
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
for i = 1:length(datsetnames_)
    v = gscatter(xpos(i,:),corrmat_(i,:),g,cmap_labs,shapes{i},7,'off');
    if strcmp('Hurtado',datsetnames_{i})
        scatter(xpos(i,:),corrmat_(i,:),75,['k' shapes{i}],'filled')
    end
    w = gscatter(NaN(1,length(g)),NaN(1,length(g)),g,'k',shapes{i},7,'off');
    for n = 1:length(v)
        if ~strcmp('Hurtado',datsetnames_{i})
          set(v(n), 'MarkerFaceColor', cmap_labs(n,:));
        end
        set(w(n), 'MarkerFaceColor', 'k');
    end
    plothands_leg = [plothands_leg, w]; %#ok<AGROW> 
end
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

set(gca,'YTick',[-1 -0.5 0 0.25 0.5 0.75])
% ylim([-0.1,0.7])
ylim([-0.1,0.4])
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