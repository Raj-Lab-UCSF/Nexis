function [classes,classmeans,pvals_2sample,pvals_1sample] = ClassViolinPlot_NexIS(outstruct,...
    subclasses_,type,cmap_,savenclose,figdir_)
if nargin < 6
    figdir_ = cd;
    if nargin < 5
        savenclose = 0;
        if nargin < 4
            cmap_ = 'hsv';
            if nargin < 3
                type = 'DeltaR2';
            end
        end
    end
end

fisher_rtoz = @(r) 0.5*(log(1+r) - log(1-r));
datsetnames_ = fieldnames(outstruct);
subclasses_outstruct = fieldnames(outstruct.(datsetnames_{1}));
corrvec_glob = NaN(length(datsetnames_),1);
corrmat_sv = NaN(length(datsetnames_),length(subclasses_));
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
corrmat_glob = repmat(corrvec_glob,1,length(subclasses_));
corrmat_ = corrmat_sv - corrmat_glob;

if length(subclasses_) == 25
    % Will finish for Tasic at a later time

elseif length(subclasses_) == 42
    classes = {'Cortical','Hippocampal','GABAergic','Non-Neuronal'};
    classlabels = [classes; {'Glutamatergic','Glutamatergic','',''}];
%     classlabels = strjust(pad(classlabels),'center');
%     classlabels = strtrim(sprintf('%s\\newline%s\n', classlabels{:}));
    classcorrcell = cell(1,length(classes)); 
    classzscorecell = cell(1,length(classes));
    classmeans = zeros(1,length(classes));
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
    for i = 1:length(classes)
        classcorr = corrmat_(:,indcell{i});
        classcorr = classcorr(:);
        classcorrcell{i} = classcorr;        
        classzscorecell{i} = fisher_rtoz(classcorr);
        classmeans(i) = mean(classcorr);
    end
end

pvals = zeros(length(classes));
pvals_1sample = zeros(1,length(classes));
for i = 1:length(classes)
    [~,p1] = ttest(classzscorecell{i});
    pvals_1sample(i) = p1 * length(classes); % Bonferroni MHC
    for j = i:length(classes)
        if i == j
            pvals(i,j) = 0.5;
        else
            [~,p2] = ttest2(classzscorecell{i},classzscorecell{j});
            pvals(i,j) = p2;
        end
    end
end
pvals = pvals + pvals.';
pvals_2sample = pvals * (length(classes)*(length(classes)-1)/2); % Bonferroni MHC
if strcmp(cmap_,'cool')
    cmap_ = cool(length(indcell));
elseif strcmp(cmap_,'hsv')
    cmap_ = hsv(length(indcell));
elseif strcmp(cmap_,'lines')
    cmap_ = lines(length(indcell));
else
    cmap_ = twocolor(cmap_(1,:),cmap_(2,:),length(indcell));
end

figure('Units','inches','Position',[0 0 16 9]); hold on;
violin(classcorrcell,'facecolor',cmap_,'medc',[]);
plot([0.5,length(classes)+0.5],[0,0],'LineStyle','--','Color',[0.25 0.25 0.25],'LineWidth',1.5);
xlim([0.5,length(classes)+0.5])

if ismember('BoludaDSAD',datsetnames_)
    ylim([-0.1,0.7])
else
    ylim([-0.1,0.4])
end
ax = gca;
for i = 1:length(classlabels)
    text(i, ax.YLim(1), sprintf('%s\n%s\n%s', classlabels{:,i}), ...
        'horizontalalignment', 'center', 'verticalalignment', 'top', 'FontName',...
        'Times','FontSize',26,'Color',cmap_(i,:));    
end
set(gca,'XTick',1:length(classes),'XTickLabel',[]);  
set(gca,'YTick',[-1 -0.5 0 0.25 0.5])
hLegend = findobj(gcf, 'Type', 'Legend');
hLegend.String = {'Class Mean'}; hLegend.FontSize = 22; hLegend.Box = 'on';
if strcmp(type,'DeltaR')
    ylabel("\DeltaR");
    figstr = [figdir_ filesep 'DeltaRClassViolin_NexIS_'];
elseif strcmp(type,'DeltaR2')
    ylabel("\DeltaR^2");
    figstr = [figdir_ filesep 'DeltaR2ClassViolin_NexIS_'];
end
if length(subclasses_) == 25
    figstr = [figstr 'Tasic'];
elseif length(subclasses_) == 42
    figstr = [figstr 'Yao'];
end
xlabel('');
set(gca,'TickLength',[0 0])
set(gca, 'FontSize', 26, 'FontName', 'Times');
if savenclose
    print(figstr,'-dtiffn','-r600'); close;
end
end