function [svals,svals_adl,svals_nadl] = CompareDirPlots_s(outstruct,pertimepoint)
    
    rng(0);
    studynames = fieldnames(outstruct);
    studynames(ismember(studynames,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
    modelnames = fieldnames(outstruct.(studynames{1}));
    isadl = ismember(studynames,{'BoludaDSAD','DS4','Hurtado'});    

    if ~pertimepoint
        svals = NaN(length(studynames),1); % fit_s only
        for i = 1:length(svals)
            resstruct = outstruct.(studynames{i}).('fit_s').nexis_global.Full;
            svals(i) = resstruct.param_fit(4);
        end
        svals_adl = svals(isadl); svals_nadl = svals(~isadl);
       
        figure; hold on;
        cmap_boxplot = [[1 0 0]; [0 0 1]];
        g = [1,2];
        xposscatter = @(y) 0.2 * (2*rand - 1) + y;
        xpos_adl = NaN(length(svals_adl),1); gvec_adl = xpos_adl;
        xpos_nadl = NaN(length(svals_nadl),1); gvec_nadl = xpos_nadl;
        for i = 1:length(xpos_adl)
            xpos_adl(i) = xposscatter(g(1));
            gvec_adl(i) = g(1);
        end
        for i = 1:length(xpos_nadl)
            xpos_nadl(i) = xposscatter(g(2));
            gvec_nadl(i) = g(2);
        end
        allsvals = [svals_adl; svals_nadl];
        gvec = [gvec_adl; gvec_nadl];
        b = boxplot(allsvals,gvec,'Colors',cmap_boxplot,'Symbol','');
        set(b,{'linew'},{2});
        scatter(xpos_adl,svals_adl,[],cmap_boxplot(1,:),'filled');
        scatter(xpos_nadl,svals_nadl,[],cmap_boxplot(2,:),'filled');
        xticks([1,2]); xlim([0.5,2.5]); xticklabels({'AD like','Not AD like'});
        yplotmax = 1.1; yplotmin = -0.1; ylim([yplotmin,yplotmax]); 
        yticks([0,0.5,1]); yticklabels({'0','0.5','1'});
        ylabel('s'); title('Longitudinal Fit');
        set(gca,'FontSize',20,'FontName','Times');

    else
        tptnames = fieldnames(outstruct.(studynames{1}).(modelnames{1}));
        svals = NaN(length(studynames),length(tptnames)); % 4 models
        for i = 1:size(svals,1)
            resstruct = outstruct.(studynames{i}).('fit_s');
            for k = 1:length(tptnames)
                resstruct_tpt = resstruct.(tptnames{k}).nexis_global.Full;
                svals(i,k) = resstruct_tpt.param_fit(4);
            end
        end
        svals_adl = svals(isadl,:); svals_nadl = svals(~isadl,:);
        svals_adl = svals_adl(:); svals_nadl = svals_nadl(:);
       
        figure; hold on;
        cmap_boxplot = [[1 0 0]; [0 0 1]];
        g = [1,2];
        xposscatter = @(y) 0.2 * (2*rand - 1) + y;
        xpos_adl = NaN(length(svals_adl),1); gvec_adl = xpos_adl;
        xpos_nadl = NaN(length(svals_nadl),1); gvec_nadl = xpos_nadl;
        for i = 1:length(xpos_adl)
            xpos_adl(i) = xposscatter(g(1));
            gvec_adl(i) = g(1);
        end
        for i = 1:length(xpos_nadl)
            xpos_nadl(i) = xposscatter(g(2));
            gvec_nadl(i) = g(2);
        end
        allsvals = [svals_adl; svals_nadl];
        gvec = [gvec_adl; gvec_nadl];
        b = boxplot(allsvals,gvec,'Colors',cmap_boxplot,'Symbol','');
        set(b,{'linew'},{2});
        scatter(xpos_adl,svals_adl,[],cmap_boxplot(1,:),'filled');
        scatter(xpos_nadl,svals_nadl,[],cmap_boxplot(2,:),'filled');
        xticks([1,2]); xlim([0.5,2.5]); xticklabels({'AD like','Not AD like'});
        yplotmax = 1.1; yplotmin = -0.1; ylim([yplotmin,yplotmax]); 
        yticks([0,0.5,1]); yticklabels({'0','0.5','1'});
        ylabel('s'); title('All Timepoints');
        set(gca,'FontSize',20,'FontName','Times');
    end
end