function [modelpreds,modelfits] = RvstPlots(outstruct,tpt_fit_,use_fits,matdir_)

studynames_ = fieldnames(outstruct);
studynames_(ismember(studynames_,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
studylabels = cellfun(@(x)strrep(x,'_',' '),studynames_,'UniformOutput',0);
% studynames_ = {'DS9_110'};
modelnames = fieldnames(outstruct.(studynames_{1}));
if ~use_fits
    modelnames = setdiff(modelnames,'fit_s');
end
modelpreds = cell(length(studynames_),length(modelnames));
modelfits = modelpreds;
tranges = modelpreds;
ntsim = 100;
for i = 1:length(studynames_)
    ts_i = outstruct.(studynames_{i}).(modelnames{1}).nexis_global.Full.time_stamps;
    C_ = outstruct.(studynames_{i}).(modelnames{1}).nexis_global.Full.init.C;
    U_ = zeros(size(C_,1),1);
    seed_i = outstruct.(studynames_{i}).(modelnames{1}).nexis_global.Full.init.seed;
    if isnan(seed_i)
        ts_i = ts_i - 2; % Hurtado correction
    end
    trange_sim = linspace(0, 1.5*ts_i(end), ntsim);
    datafit_i = outstruct.(studynames_{i}).(modelnames{1}).nexis_global.Full.data(:,tpt_fit_);
    for j = 1:length(modelnames)
        resstruct_ij = outstruct.(studynames_{i}).(modelnames{j}).nexis_global.Full;
        params_ij = resstruct_ij.param_fit;
        if isnan(seed_i)
            params_ij(1) = 1; % should be redundant
            seed_i = resstruct_ij.baseline;
            t0 = 2;
        else
            t0 = 0;
        end
        seed_i_ccf = DataToCCF(seed_i,studynames_{i},matdir_);
        seed_i_ccf(isnan(seed_i_ccf)) = 0;
        preds_ij_ccf = NexIS_fun(C_,U_,trange_sim,seed_i_ccf,params_ij,'analytic',1,matdir_);
        preds_ij = CCFToData(preds_ij_ccf,studynames_{i},matdir_);
        modelpreds{i,j} = preds_ij;
        corrXmat_ij = [datafit_i, preds_ij];
        Rvals_ij = corr(corrXmat_ij);
        modelfits{i,j} = Rvals_ij(1,2:end);
        tranges{i,j} = trange_sim + t0;
    end
end

cmap = hsv(length(modelnames));
figure('Units','inches','Position',[0 0 20 20]); 
tiledlayout(3,4,'TileSpacing','compact','Padding','tight');
for i = 1:length(studynames_)
    nexttile; hold on;
    mins_i = NaN(1,length(modelnames)); maxs_i = mins_i;
    for j = 1:length(modelnames)
        plot(tranges{i,j},modelfits{i,j},'Color',cmap(j,:),'LineWidth',3);
        mins_i(j) = min(modelfits{i,j}); maxs_i(j) = max(modelfits{i,j});
    end
    plotmax_i = 1.1*max(maxs_i);
    if sign(min(mins_i)) == 1
        plotmin_i = 0.9*min(mins_i);
    else
        plotmin_i = 1.1*min(mins_i);
    end
    ts_i = outstruct.(studynames_{i}).(modelnames{1}).nexis_global.Full.time_stamps;
    tfit_i = ts_i(tpt_fit_);
    plot([tfit_i,tfit_i], [plotmin_i,plotmax_i],'LineWidth',2,'LineStyle',':','Color','k');
    legnams = {'Ret','Ant','N.D.'};
    if use_fits
        s_i = outstruct.(studynames_{i}).(modelnames{1}).nexis_global.Full.param_fit(4);
        fits_str = sprintf('s = %.2f',s_i);
        legnams = [fits_str, legnams];
    end
    legend(legnams,'Location','southwest')
    title(studylabels{i});
    set(gca,'FontName','Times','FontSize',16);
end