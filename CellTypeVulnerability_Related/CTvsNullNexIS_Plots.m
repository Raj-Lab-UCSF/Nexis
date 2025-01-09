function xcell = CTvsNullNexIS_Plots(outstruct_ct,outstruct_null,subclasses_,type,cmap_,savenclose_,figdir_)

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
        if strcmp(type,'DeltaR2')
            if j == 1
                corrvec_glob_ct(i) = outstruct_i_ct.(subclasses_outstruct_ct{j}).nexis_global.Full.results.lm_Rsquared_adj;
            else
                corrmat_sv_ct(i,(j-1)) = outstruct_i_ct.(subclasses_outstruct_ct{j}).nexis_sv.Full.results.lm_Rsquared_adj;
            end
        elseif strcmp(type,'DeltaR')
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
        if strcmp(type,'DeltaR2')
            if j == 1
                corrvec_glob_null(i) = outstruct_i_null.(subclasses_outstruct_null{j}).nexis_global.Full.results.lm_Rsquared_adj;
            else
                corrmat_sv_null(i,(j-1)) = outstruct_i_null.(subclasses_outstruct_null{j}).nexis_sv.Full.results.lm_Rsquared_adj;
            end
        elseif strcmp(type,'DeltaR')
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
pval = 0.01;
thresh = 100*(1 - pval);
Rthresh = prctile(corrmat_null,thresh,1);

cmap_ct = zeros(length(indtest),3);
for i = 1:length(indtest)
    colorinds = zeros(1,length(indcell));
    for j = 1:length(indcell)
        colorinds(j) = ismember(i,indcell{j});
    end
    cmap_ct(i,:) = cmap_labs(logical(colorinds),:);
end
cmap_labs = cmap_ct(indtest,:);
corrmat_ct = corrmat_ct(:,indtest);
subclasses_ = subclasses_(indtest);
cmap_hist = hsv(size(corrmat_null,2));
% cmap_hist = hsv(size(corrmat_null,2)+1);
% cmap_hist(1,:) = [];

if length(datsetnames_) == 9
    f = figure('Units','inches','Position',[0,0,11,10]);
else
    f = figure('Units','inches','Position',[0,0,17,15]);
end
datsetlabels = cellfun(@(x) strrep(x,'_',' '), datsetnames_,'UniformOutput',0);
xcell = cell(1,length(datsetnames_));
for i = 1:length(datsetnames_)
    if length(datsetnames_) == 9
        subplot(3,3,i); hold on;
    else
        subplot(3,4,i); hold on;
    end
    % plothands_leg = [];
    % subclasses_leg = {};
    % xvec_sig = Rthresh(i);
    xvec_all = corrmat_null(:,i).';
    histogram(corrmat_null(:,i),10,'FaceAlpha',0.4,'FaceColor',cmap_hist(i,:));
    if length(corrmat_null(:,i)) == 100
        maxylim = 30;
    elseif length(corrmat_null(:,i)) == 500
        maxylim = 200;
    end
    w = plot([Rthresh(i),Rthresh(i)],[0,maxylim],...
                'Color','r','LineStyle',':','LineWidth',3);
    if strcmp(type,'DeltaR2')
        legstr = '\DeltaR^2';
    else
        legstr = '\DeltaR';
    end
    legend(w,[legstr ' = ' num2str(Rthresh(i),'%.2f')],'Location','northwest','FontSize',16,'box','off');
    % for j = 1:length(subclasses_)
        % if Rthresh(i) < corrmat_ct(i,j)
        %     w = plot([corrmat_ct(i,j),corrmat_ct(i,j)],[0,17],...
        %         'Color',cmap_labs(j,:),'LineStyle',':','LineWidth',3);
            % plothands_leg = [plothands_leg, w];
            % subclasses_leg = [subclasses_leg, subclasses_{j}];
            % xvec_sig = [xvec_sig, corrmat_ct(i,j)];
            % xvec_all = [xvec_all, corrmat_ct(i,j)];
        % end
        % xcell{i} = xvec_sig;
    % end
    % legend(plothands_leg,subclasses_leg,'Location','northwest','FontSize',16,'box','off');
    max_x = 1.1*max(xvec_all); min_x = 0.9*min(xvec_all); xlim([min_x max_x]); 
    xticks([min_x (max_x+min_x)/2 max_x]); 
    max_y = maxylim; min_y = 0; ylim([min_y max_y]); yticks([min_y max_y]);
    if round(min_x,1) ~= 0 
        xticklabels({num2str(min_x,'%.2f'),num2str((max_x+min_x)/2,'%.2f'),...
            num2str(max_x,'%.2f')});
    else
        xticklabels({'0',num2str((max_x+min_x)/2,'%.2f'),...
            num2str(max_x,'%.2f')});
    end
    % if round(min_y,1) ~= 0 
    %     yticklabels({num2str(min_y,'%.1f'),num2str((max_y+min_y)/2,'%.1f'),...
    %         num2str(max_y,'%.1f')});
    % else
    %     yticklabels({'0',num2str((max_y+min_y)/2,'%.1f'),...
    %         num2str(max_y,'%.1f')});
    % end

    % text(0.05*(max_x-min_x)+min_x, 0.85*(max_y-min_y)+min_y,...
    %     {sprintf('R^2 = %.2f%s',mdls_bic{i}.Rsquared.Adjusted,pstr),...
    %     sprintf('n = %d',length(subclasses_bic{i}))},'FontName','Times','FontSize',18);

    title(datsetlabels{i},'FontSize',24);
    set(gca,'FontSize',20,'FontName','Times','box','on');
end
han=axes(f,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
if strcmp(type,'DeltaR2')
    xlabstr = '\DeltaR^2';
else
    xlabstr = '\DeltaR';
end
ylh = ylabel(han,'Frequency','FontSize',24,...
    'FontName','Times','FontWeight','bold');
ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 1);
xlh = xlabel(han,xlabstr,'FontSize',24,...
    'FontName','Times','FontWeight','bold');
xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 1);
if savenclose_
    if strcmp(tpt_,'end')
        print([figdir_ filesep 'LinearModelPlots_Coeffs_' critstr '_' datasource '_t_' tpt_],'-dtiffn','-r600'); close;
    else
        print([figdir_ filesep 'LinearModelPlots_Coeffs_' critstr '_' datasource '_t' num2str(tpt_)],'-dtiffn','-r600'); close;
    end
end
end