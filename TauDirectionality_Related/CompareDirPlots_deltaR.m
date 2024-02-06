function Rmat = CompareDirPlots_deltaR(outstruct,pertimepoint,RvR2)

    studynames = fieldnames(outstruct);
    studynames(ismember(studynames,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
    modelnames = fieldnames(outstruct.(studynames{1}));
    isadl = ismember(studynames,{'BoludaDSAD','DS4','Hurtado'});    

    if ~pertimepoint

        Rmat = NaN(length(studynames),length(modelnames)); % 4 models
        for i = 1:size(Rmat,1)
            for j = 1:size(Rmat,2)
                resstruct = outstruct.(studynames{i}).(modelnames{j}).nexis_global.Full;
                if strcmp(RvR2,'R2')
                    Rmat(i,j) = resstruct.results.lm_Rsquared_adj;
                else
                    datavec = resstruct.data(:);
                    predvec = resstruct.predicted(:);
                    Rmat(i,j) = corr(datavec,predvec,'rows','complete');
                end
            end
        end
        % fitsind = find(ismember(modelnames,'fit_s'));
        retind = find(ismember(modelnames,'ret'));
        antind = find(ismember(modelnames,'ant'));
        % ndind = find(ismember(modelnames,'nd'));
        Rdiffs = Rmat(:,retind) - Rmat(:,antind); %#ok<FNDSB>
        adlRdiffs = Rdiffs(isadl); nadlRdiffs = Rdiffs(~isadl);
        figure; hold on;
        scatter(ones(1,length(adlRdiffs)),adlRdiffs,'bo','filled');
        scatter(2*ones(1,length(nadlRdiffs)),nadlRdiffs,'ro','filled');
        xticks([1,2]); xticklabels({'AD like','Not AD like'});
        xlim([0.5,2.5]); set(gca,'FontSize',16,'FontName','Times')

    else
        tptnames = fieldnames(outstruct.(studynames{1}).(modelnames{1}));
        Rmat = NaN(length(studynames),length(modelnames)*length(tptnames)); % 4 models
        for i = 1:size(Rmat,1)
            for j = 1:length(modelnames)
                resstruct = outstruct.(studynames{i}).(modelnames{j});
                for k = 1:length(tptnames)
                    resstruct_tpt = resstruct.(tptnames{k}).nexis_global.Full;
                    Rind = k + length(tptnames)*(j-1);
                    if strcmp(RvR2,'R2')
                        Rmat(i,Rind) = resstruct_tpt.results.lm_Rsquared_adj;
                    else
                        datavec = resstruct_tpt.data(:);
                        predvec = resstruct_tpt.predicted(:);
                        Rmat(i,Rind) = corr(datavec,predvec,'rows','complete');
                    end
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
        scatter(ones(1,length(adlRdiffs)),adlRdiffs,'bo','filled');
        scatter(2*ones(1,length(nadlRdiffs)),nadlRdiffs,'ro','filled');
        xticks([1,2]); xticklabels({'AD like','Not AD like'});
        xlim([0.5,2.5]); set(gca,'FontSize',16,'FontName','Times')


    end
end