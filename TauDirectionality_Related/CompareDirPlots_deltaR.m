function Rmat = CompareDirPlots_deltaR(outstruct,pertimepoint)
    
    rng(0);
    studynames = fieldnames(outstruct);
    studynames(ismember(studynames,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
    modelnames = fieldnames(outstruct.(studynames{1}));
    isadl = ismember(studynames,{'BoludaDSAD','DS4','Hurtado'});    

    if ~pertimepoint

        Rmat = NaN(length(studynames),length(modelnames)); % 4 models
        for i = 1:size(Rmat,1)
            for j = 1:size(Rmat,2)
                resstruct = outstruct.(studynames{i}).(modelnames{j}).nexis_global.Full;
                % if strcmp(RvR2,'R2')
                %     Rmat(i,j) = resstruct.results.lm_Rsquared_adj;
                % else
                    datavec = resstruct.data(:);
                    predvec = resstruct.predicted(:);
                    Rmat(i,j) = corr(datavec,predvec,'rows','complete');
                % end
            end
        end
        % fitsind = find(ismember(modelnames,'fit_s'));
        retind = find(ismember(modelnames,'ret'));
        antind = find(ismember(modelnames,'ant'));
        % ndind = find(ismember(modelnames,'nd'));
        Rdiffs = Rmat(:,retind) - Rmat(:,antind); %#ok<FNDSB>
        adlRdiffs = Rdiffs(isadl); nadlRdiffs = Rdiffs(~isadl);
       
        figure; hold on;
        cmap_boxplot = [[1 0 0]; [0 0 1]];
        g = [1,2];
        xposscatter = @(y) 0.2 * (2*rand - 1) + y;
        xpos_adl = NaN(length(adlRdiffs),1); gvec_adl = xpos_adl;
        xpos_nadl = NaN(length(nadlRdiffs),1); gvec_nadl = xpos_nadl;
        for i = 1:length(xpos_adl)
            xpos_adl(i) = xposscatter(g(1));
            gvec_adl(i) = g(1);
        end
        for i = 1:length(xpos_nadl)
            xpos_nadl(i) = xposscatter(g(2));
            gvec_nadl(i) = g(2);
        end
        alldiffs = [adlRdiffs; nadlRdiffs];
        gvec = [gvec_adl; gvec_nadl];
        b = boxplot(alldiffs,gvec,'Colors',cmap_boxplot,'Symbol','');
        set(b,{'linew'},{2});
        scatter(xpos_adl,adlRdiffs,[],cmap_boxplot(1,:),'filled');
        scatter(xpos_nadl,nadlRdiffs,[],cmap_boxplot(2,:),'filled');
        xticks([1,2]); xlim([0.5,2.5]); xticklabels({'AD like','Not AD like'});
        yplotmax1 = max(adlRdiffs(:)); yplotmax2 = max(nadlRdiffs(:));
        yplotmax = max([yplotmax1, yplotmax2]) + 0.05;
        yplotmin1 = min(adlRdiffs(:)); yplotmin2 = min(nadlRdiffs(:));
        yplotmin = min([yplotmin1, yplotmin2]) - 0.05;
        ylim([yplotmin,yplotmax]); yticks([yplotmin,(yplotmin+yplotmax)/2,yplotmax]);
        ytickformat('%.2f');
        ylabel('R_r_e_t - R_a_n_t'); title('Longitudinal Fit');
        set(gca,'FontSize',20,'FontName','Times');

    else
        tptnames = fieldnames(outstruct.(studynames{1}).(modelnames{1}));
        Rmat = NaN(length(studynames),length(modelnames)*length(tptnames)); % 4 models
        for i = 1:size(Rmat,1)
            for j = 1:length(modelnames)
                resstruct = outstruct.(studynames{i}).(modelnames{j});
                for k = 1:length(tptnames)
                    resstruct_tpt = resstruct.(tptnames{k}).nexis_global.Full;
                    Rind = k + length(tptnames)*(j-1);
                    % if strcmp(RvR2,'R2')
                    %     Rmat(i,Rind) = resstruct_tpt.results.lm_Rsquared_adj;
                    % else
                        datavec = resstruct_tpt.data(:);
                        predvec = resstruct_tpt.predicted(:);
                        Rmat(i,Rind) = corr(datavec,predvec,'rows','complete');
                    % end
                end
            end
        end
        modelnames_ind = repmat(modelnames.',length(tptnames),1);
        modelnames_ind = modelnames_ind(:);
        % fitsinds = find(ismember(modelnames_ind,'fit_s'));
        retinds = find(ismember(modelnames_ind,'ret'));
        antinds = find(ismember(modelnames_ind,'ant'));
        % ndind = find(ismember(modelnames_ind,'nd'));
        Rdiffs = Rmat(:,retinds) - Rmat(:,antinds); %#ok<FNDSB>
        adlRdiffs = Rdiffs(isadl,:); nadlRdiffs = Rdiffs(~isadl,:);
        adlRdiffs = adlRdiffs(:); nadlRdiffs = nadlRdiffs(:);
       
        figure; hold on;
        cmap_boxplot = [[1 0 0]; [0 0 1]];
        g = [1,2];
        xposscatter = @(y) 0.2 * (2*rand - 1) + y;
        xpos_adl = NaN(length(adlRdiffs),1); gvec_adl = xpos_adl;
        xpos_nadl = NaN(length(nadlRdiffs),1); gvec_nadl = xpos_nadl;
        for i = 1:length(xpos_adl)
            xpos_adl(i) = xposscatter(g(1));
            gvec_adl(i) = g(1);
        end
        for i = 1:length(xpos_nadl)
            xpos_nadl(i) = xposscatter(g(2));
            gvec_nadl(i) = g(2);
        end
        alldiffs = [adlRdiffs; nadlRdiffs];
        gvec = [gvec_adl; gvec_nadl];
        b = boxplot(alldiffs,gvec,'Colors',cmap_boxplot,'Symbol','');
        set(b,{'linew'},{2});
        scatter(xpos_adl,adlRdiffs,[],cmap_boxplot(1,:),'filled');
        scatter(xpos_nadl,nadlRdiffs,[],cmap_boxplot(2,:),'filled');
        xticks([1,2]); xlim([0.5,2.5]); xticklabels({'AD like','Not AD like'});
        yplotmax1 = max(adlRdiffs(:)); yplotmax2 = max(nadlRdiffs(:));
        yplotmax = max([yplotmax1, yplotmax2]) + 0.05;
        yplotmin1 = min(adlRdiffs(:)); yplotmin2 = min(nadlRdiffs(:));
        yplotmin = min([yplotmin1, yplotmin2]) - 0.05;
        ylim([yplotmin,yplotmax]); yticks([yplotmin,(yplotmin+yplotmax)/2,yplotmax]);
        ytickformat('%.2f');
        ylabel('R_r_e_t - R_a_n_t'); title('All Timepoints');
        set(gca,'FontSize',20,'FontName','Times');



    end
end