function CellTypeFrequencyPlot_NexIS(outstruct_ct,outstruct_null,subclasses_,...
    type_,thresh_pval_,cmap_,savenclose_,figdir_)

datsetnames_ = fieldnames(outstruct_ct);
if (length(subclasses_) == 42)
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
else
    % Figure out later for Tasic?
    % if strcmp(cmap_,'cool')
    %     cmap_labs = cool(length(datsetnames_));
    % elseif strcmp(cmap_,'hsv')
    %     cmap_labs = hsv(length(datsetnames_));
    % elseif strcmp(cmap_,'lines')
    %     cmap_labs = lines(length(datsetnames_));
    % else
    %     cmap_labs = twocolor(cmap_(1,:),cmap_(2,:),length(datsetnames_));
    % end
    % indtest = 1:length(datsetnames_);
    % indcell = cell(1,length(datsetnames_));
    % for i = 1:length(indcell)   
    %     indcell{i} = i;
    % end
end

subclasses_outstruct_ct = fieldnames(outstruct_ct.(datsetnames_{1}));
subclasses_outstruct_null = fieldnames(outstruct_null.(datsetnames_{1}));
corrvec_glob_ct = NaN(length(datsetnames_),1);
corrmat_sv_ct = NaN(length(datsetnames_),(length(subclasses_outstruct_ct)-1));
corrvec_glob_null = NaN(length(datsetnames_),1); % should be the same as ct
corrmat_sv_null = NaN(length(datsetnames_),(length(subclasses_outstruct_null)-1));

for i = 1:length(datsetnames_)
    datsetname_ = datsetnames_{i};
    outstruct_i_ct = outstruct_ct.(datsetname_);
    outstruct_i_null = outstruct_null.(datsetname_);
    for j = 1:length(subclasses_outstruct_ct)
        if strcmp(type_,'DeltaR2')
            if j == 1
                corrvec_glob_ct(i) = outstruct_i_ct.(subclasses_outstruct_ct{j}).nexis_global.Full.results.lm_Rsquared_adj;
            else
                corrmat_sv_ct(i,(j-1)) = outstruct_i_ct.(subclasses_outstruct_ct{j}).nexis_sv.Full.results.lm_Rsquared_adj;
            end
        elseif strcmp(type_,'DeltaR')
            if j == 1
                R2val_ij = outstruct_i_ct.(subclasses_outstruct_ct{j}).nexis_global.Full.results.lm_Rsquared_ord;
                corrvec_glob_ct(i) = R2val_ij^(0.5);
            else
                R2val_ij = outstruct_i_ct.(subclasses_outstruct_ct{j}).nexis_sv.Full.results.lm_Rsquared_ord;
                corrmat_sv_ct(i,(j-1)) = R2val_ij^(0.5);
            end
        end
    end
    for j = 1:length(subclasses_outstruct_null)
        if strcmp(type_,'DeltaR2')
            if j == 1
                corrvec_glob_null(i) = outstruct_i_null.(subclasses_outstruct_null{j}).nexis_global.Full.results.lm_Rsquared_adj;
            else
                corrmat_sv_null(i,(j-1)) = outstruct_i_null.(subclasses_outstruct_null{j}).nexis_sv.Full.results.lm_Rsquared_adj;
            end
        elseif strcmp(type_,'DeltaR')
            if j == 1
                R2val_ij = outstruct_i_null.(subclasses_outstruct_null{j}).nexis_global.Full.results.lm_Rsquared_ord;
                corrvec_glob_null(i) = R2val_ij^(0.5);
            else
                R2val_ij = outstruct_i_null.(subclasses_outstruct_null{j}).nexis_sv.Full.results.lm_Rsquared_ord;
                corrmat_sv_null(i,(j-1)) = R2val_ij^(0.5);
            end
        end
    end
end

corrmat_glob_ct = repmat(corrvec_glob_ct,1,size(corrmat_sv_ct,2));
corrmat_ct = corrmat_sv_ct - corrmat_glob_ct;
corrmat_glob_null = repmat(corrvec_glob_null,1,size(corrmat_sv_null,2));
corrmat_null = corrmat_sv_null - corrmat_glob_null;

corrmat_null = corrmat_null.';
pval = thresh_pval_;
thresh = 100*(1 - pval);
Rthresh = prctile(corrmat_null,thresh,1);
for i = 1:length(datsetnames_)
    corrvec_ct = corrmat_ct(i,:);
    corrvec_ct(corrvec_ct < Rthresh(i)) = 0;
    corrmat_ct(i,:) = corrvec_ct;
end
corrmat_ct = logical(corrmat_ct);
corrmat_ct = corrmat_ct(:,indtest);

subclasses_ = subclasses_(indtest);
subclasses_1d = {};
for i = 1:length(datsetnames_)
    Rinds = corrmat_ct(i,:);
    subclasses_sig_i = subclasses_(Rinds);
    subclasses_1d = [subclasses_1d; subclasses_sig_i]; %#ok<AGROW> 
end
subclasses_sig = unique(subclasses_1d);
ctfreq = zeros(1,length(subclasses_sig));
for i = 1:length(ctfreq)
    ctfreq(i) = sum(ismember(subclasses_1d,subclasses_sig{i}));
end
[ctfreq_sort,sortinds] = sort(ctfreq,'descend');
subsigs_sort = subclasses_sig(sortinds);

cmap_ct = zeros(length(indtest),3);
for i = 1:length(indtest)
    colorinds = zeros(1,length(indcell));
    for j = 1:length(indcell)
        colorinds(j) = ismember(i,indcell{j});
    end
    cmap_ct(i,:) = cmap_labs(logical(colorinds),:);
end
cmap_labs = cmap_ct(indtest,:);

figure('Units','inches','Position',[0 0 8 15]); 
b = barh(ctfreq_sort/length(datsetnames_),'EdgeColor',[0.25 0.25 0.25],'LineWidth',1);
b.FaceColor = 'flat';
ylabs = cell(1,length(subsigs_sort));
for i = 1:length(subsigs_sort)
    col = cmap_labs(ismember(subclasses_,subsigs_sort{i}),:);
    ylabs{i} = sprintf('\\color[rgb]{%f,%f,%f}%s',col(1),col(2),col(3),subsigs_sort{i});
    b.CData(i,:) = col;
end
maxx = max(ctfreq/length(datsetnames_));
xlabel('')
set(gca,'YTick',1:length(subsigs_sort),'YTickLabel',ylabs,'TickLength',[0 1],...
    'XTick',[0 maxx/2 maxx],'XTickLabel',{'0',num2str(maxx/2,'%.2f'),...
    num2str(maxx,'%.2f')},'TickLabelInterpreter','tex','FontName',...
    'Times','FontSize',20,'YDir','reverse');
if strcmp(type_,'DeltaR2') 
    xlabel('NexIS Selection Frequency, \DeltaR^2','FontWeight','bold','FontName','Times','FontSize',24); 
elseif strcmp(type_,'DeltaR') 
    xlabel('NexIS Selection Frequency, \DeltaR','FontWeight','bold','FontName','Times','FontSize',24); 
end
threshtext = sprintf('p < %.1d',thresh_pval_);
text(0.7,0.05,threshtext,'Units','normalized','FontSize',20,'FontName','Times');
if savenclose_
    if strcmp(type_,'DeltaR2') 
        print([figdir_ filesep 'CellTypeFrequencyPlot_R2_Yao'],'-dtiffn','-r600'); close;
    elseif strcmp(type_,'DeltaR') 
        print([figdir_ filesep 'CellTypeFrequencyPlot_R_Yao'],'-dtiffn','-r600'); close;
    end
end
end