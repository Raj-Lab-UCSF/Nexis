function CorrComparePlot(outstruct,pertimepoint)
     
    studynames = fieldnames(outstruct);
    studynames(ismember(studynames,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
    studylabels = cellfun(@(x)strrep(x,'_',' '),studynames,'UniformOutput',0);
    modelnames = fieldnames(outstruct.(studynames{1}));
    isadl = ismember(studynames,{'BoludaDSAD','DS4','Hurtado'});    

    if ~pertimepoint
        Rmat = NaN(length(studynames),length(modelnames)); % 4 models
        for i = 1:size(Rmat,1)
            for j = 1:size(Rmat,2)
                resstruct = outstruct.(studynames{i}).(modelnames{j}).nexis_global.Full;
                datavec = resstruct.data(:);
                predvec = resstruct.predicted(:);
                Rmat(i,j) = corr(datavec,predvec,'rows','complete');
            end
        end
        shapes = {'o','s','d','^'};
        cmap = hsv(length(studynames));
        xpos = 1:length(studynames); g = xpos;
        offsets = [-0.3,-0.1,0.1,0.3];
        figure('Units','inches','Position',[0 0 15 8]); hold on;
        legstr = []; xlabs = cell(1,length(xpos));
        for i = 1:length(shapes)
            gplot1 = gscatter((xpos.' + offsets(i)),Rmat(:,i),g,zeros(size(cmap)),shapes{i},7,'doleg','off');
            gplot2 = gscatter((xpos.' + offsets(i)),Rmat(:,i),g,cmap,shapes{i},7,'doleg','off');
            for j = xpos
                col = cmap(j,:);
                gplot1(j).MarkerFaceColor = [0 0 0];
                gplot2(j).MarkerFaceColor = col;
                xlabs{j} = sprintf('\\color[rgb]{%f,%f,%f}%s',col(1),col(2),col(3),studylabels{j});
            end
            legstr = [legstr gplot1(1)];
        end
        xticks(xpos); xlim([min(xpos)-0.5,max(xpos)+0.5]); xticklabels(xlabs);
        xlabel([]);
        yplotmax = max(Rmat(:)); yplotmin = 0;
        ylim([yplotmin,1.1*yplotmax]); yticks([yplotmin,(yplotmin+yplotmax)/2,yplotmax]);
        yticklabels({'0',num2str((yplotmin+yplotmax)/2,'%.2f'),num2str(yplotmax,'%.2f')})
        ylabel('R'); title('Longitudinal Fit');
        legend(legstr,{'fit s','ret','ant','nd'},'Location','southeast','NumColumns',4,'FontSize',22,'box','off')
        set(gca,'FontSize',26,'FontName','Times');

    else
        tptnames = fieldnames(outstruct.(studynames{1}).(modelnames{1}));
        Rmat = NaN(length(studynames),length(modelnames)*length(tptnames)); % 4 models
        for i = 1:size(Rmat,1)
            for j = 1:length(modelnames)
                resstruct = outstruct.(studynames{i}).(modelnames{j});
                for k = 1:length(tptnames)
                    resstruct_tpt = resstruct.(tptnames{k}).nexis_global.Full;
                    Rind = k + length(tptnames)*(j-1);
                    datavec = resstruct_tpt.data(:);
                    predvec = resstruct_tpt.predicted(:);
                    Rmat(i,Rind) = corr(datavec,predvec,'rows','complete');
                end
            end
        end

        shapes = {'o','s','d','^'};
        cmap = hsv(length(studynames));
        cmapfull = NaN(size(cmap,1)*length(tptnames),3);
        Rvecs = NaN(size(cmap,1)*length(tptnames),length(modelnames));
        for i = 1:size(cmap,1)
            inds = (1:length(tptnames)) + (i-1)*length(tptnames);
            cmapfull(inds,:) = repmat(cmap(i,:),length(tptnames),1);
        end
        for i = 1:length(modelnames)
            inds = (1:length(tptnames)) + (i-1)*length(tptnames);
            Rvec = Rmat(:,inds).';
            Rvecs(:,i) = Rvec(:);
        end
        xpos = 1:length(studynames); 
        xposvec = repmat(xpos,length(tptnames),1); 
        xposvec = xposvec(:); g = xposvec;
        offsets = [-0.3,-0.1,0.1,0.3];
        figure('Units','inches','Position',[0 0 15 8]); hold on;
        legstr = []; xlabs = cell(1,length(xpos));
        for i = 1:length(shapes)
            gplot1 = gscatter((xposvec.' + offsets(i)),Rvecs(:,i),g,zeros(size(cmapfull)),shapes{i},7,'doleg','off');
            gplot2 = gscatter((xposvec.' + offsets(i)),Rvecs(:,i),g,cmapfull,shapes{i},7,'doleg','off');
            for j = 1:length(xpos)
                col = cmap(j,:);
                gplot1(j).MarkerFaceColor = [0 0 0];
                gplot2(j).MarkerFaceColor = col; gplot2(j).MarkerEdgeColor = col;
                xlabs{j} = sprintf('\\color[rgb]{%f,%f,%f}%s',col(1),col(2),col(3),studylabels{j});
            end
            legstr = [legstr gplot1(1)];
        end
        xticks(xpos); xlim([min(xposvec)-0.5,max(xposvec)+0.5]); xticklabels(xlabs);
        xlabel([]);
        yplotmax = max(Rmat(:)); yplotmin = 0;
        ylim([yplotmin,1.1*yplotmax]); yticks([yplotmin,(yplotmin+yplotmax)/2,yplotmax]);
        yticklabels({'0',num2str((yplotmin+yplotmax)/2,'%.2f'),num2str(yplotmax,'%.2f')})
        ylabel('R'); title('Per Timepoint');
        legend(legstr,{'fit s','ret','ant','nd'},'Location','southeast','NumColumns',4,'FontSize',22,'box','off')
        set(gca,'FontSize',26,'FontName','Times');
        % figure; hold on;
        % cmap_boxplot = [[1 0 0]; [0 0 1]];
        % g = [1,2];
        % xposscatter = @(y) 0.2 * (2*rand - 1) + y;
        % xpos_adl = NaN(length(adlRdiffs),1); gvec_adl = xpos_adl;
        % xpos_nadl = NaN(length(nadlRdiffs),1); gvec_nadl = xpos_nadl;
        % for i = 1:length(xpos_adl)
        %     xpos_adl(i) = xposscatter(g(1));
        %     gvec_adl(i) = g(1);
        % end
        % for i = 1:length(xpos_nadl)
        %     xpos_nadl(i) = xposscatter(g(2));
        %     gvec_nadl(i) = g(2);
        % end
        % alldiffs = [adlRdiffs; nadlRdiffs];
        % gvec = [gvec_adl; gvec_nadl];
        % b = boxplot(alldiffs,gvec,'Colors',cmap_boxplot,'Symbol','');
        % set(b,{'linew'},{2});
        % scatter(xpos_adl,adlRdiffs,[],cmap_boxplot(1,:),'filled');
        % scatter(xpos_nadl,nadlRdiffs,[],cmap_boxplot(2,:),'filled');
        % xticks([1,2]); xlim([0.5,2.5]); xticklabels({'AD like','Not AD like'});
        % yplotmax1 = max(adlRdiffs(:)); yplotmax2 = max(nadlRdiffs(:));
        % yplotmax = max([yplotmax1, yplotmax2]) + 0.05;
        % yplotmin1 = min(adlRdiffs(:)); yplotmin2 = min(nadlRdiffs(:));
        % yplotmin = min([yplotmin1, yplotmin2]) - 0.05;
        % ylim([yplotmin,yplotmax]); yticks([yplotmin,(yplotmin+yplotmax)/2,yplotmax]);
        % ytickformat('%.2f');
        % ylabel('R_r_e_t - R_a_n_t'); title('All Timepoints');
        % set(gca,'FontSize',20,'FontName','Times');
    end

end