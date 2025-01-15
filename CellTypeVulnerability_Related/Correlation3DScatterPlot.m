function Correlation3DScatterPlot(corrmat_ct_,corrmat_res_,corrmat_nx_,subclasses_,...
    corrtype_nx_,cmap_,savenclose_,figdir_)

if size(corrmat_ct_,2) == 42
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
    % indtest = [glutctx_inds;gluthipp_inds;gaba_inds;nonneuronal_inds];
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
    if strcmp(cmap_,'cool')
        cmap_labs = cool(size(corrmat_ct_,1));
    elseif strcmp(cmap_,'hsv')
        cmap_labs = hsv(size(corrmat_ct_,1));
    elseif strcmp(cmap_,'lines')
        cmap_labs = lines(size(corrmat_ct_,1));
    else
        cmap_labs = twocolor(cmap_(1,:),cmap_(2,:),size(corrmat_ct_,1));
    end
    % indtest = 1:length(subclasses_);
    indcell = cell(1,size(corrmat_ct_,1));
    for i = 1:length(indcell)   
        indcell{i} = i;
    end
end

corrvec_ct_mean = mean(corrmat_ct_);
corrvec_res_mean = mean(corrmat_res_);
corrvec_nx_mean = mean(corrmat_nx_);

% cmap_ = zeros(length(indtest),3);
% for i = 1:length(indtest)
%     colorinds = zeros(1,length(indcell));
%     for j = 1:length(indcell)
%         colorinds(j) = ismember(i,indcell{j});
%     end
%     cmap_(i,:) = cmap_labs(logical(colorinds),:);
% end
% cmap_labs = cmap_(indtest,:);
% corrmat_ = corrmat_(:,indtest);
% subclasses_ = subclasses_(indtest);

cmap_ = cmap_labs;
figure('Units','inches','Position',[0 0 8 8]); hold on;
for i = 1:size(cmap_,1)
    ct_inds = indcell{i};
    scatter3(corrvec_ct_mean(ct_inds),corrvec_res_mean(ct_inds),corrvec_nx_mean(ct_inds),...
        75,cmap_(i,:),'filled');
end
if strcmp(corrtype_nx_,'DeltaR2')
    corrtype_nx_str = '\DeltaR^2';
elseif strcmp(corrtype_nx_,'DeltaR')
    corrtype_nx_str = '\DeltaR';
end
xlabel('R_o_b_s'); ylabel('R_r_e_s'); zlabel(corrtype_nx_str);
view([135,35.2644]);
set(gca,'FontSize',24,'FontName','Times','box','on');
if savenclose_
    figstr = [figdir_ filesep '3d_scatter_' corrtype_nx_];
    print(figstr,'-dtiffn','-r600'); close;
end
end