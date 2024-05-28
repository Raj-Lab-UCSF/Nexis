function CorrComparePlot_Combined(outstruct_all,outstruct_tpt,usepertpt_,savenclose_,figdir_)
 
studynames = fieldnames(outstruct_all);
studynames(ismember(studynames,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
studylabels = cellfun(@(x)strrep(x,'_',' '),studynames,'UniformOutput',0);
modelnames = fieldnames(outstruct_all.(studynames{1}));

Rmat_all = NaN(length(studynames),length(modelnames)); % 4 models
if ~usepertpt_
    for i = 1:size(Rmat_all,1)
        for j = 1:size(Rmat_all,2)
            resstruct = outstruct_all.(studynames{i}).(modelnames{j}).nexis_global.Full;
            Rmat_all(i,j) = (resstruct.results.lm_Rsquared_ord)^(0.5);
        end
    end
end

if usepertpt_
    tptnames = fieldnames(outstruct_tpt.(studynames{1}).(modelnames{1}));
else
    tptnames = {'t_1','t_2','t_3'};
end
Rmat_tpt = NaN(length(studynames),length(modelnames)*length(tptnames)); % 4 models
for i = 1:size(Rmat_tpt,1)
    for j = 1:length(modelnames)
        if usepertpt_
            resstruct = outstruct_tpt.(studynames{i}).(modelnames{j});
            corrs_j = NaN(1,length(tptnames));
            for k = 1:length(tptnames)
                Rind = k + length(tptnames)*(j-1);
                resstruct_tpt = resstruct.(tptnames{k}).nexis_global.Full;
                corrs_j(k) = resstruct_tpt.results.Corrs;
                Rmat_tpt(i,Rind) = corrs_j(k);
            end
            Rmat_all(i,j) = mean(corrs_j);
        else
            resstruct = outstruct_all.(studynames{i}).(modelnames{j});
            for k = 1:length(tptnames)
                Rind = k + length(tptnames)*(j-1);
                resstruct_tpt = resstruct.nexis_global.Full;
                Rmat_tpt(i,Rind) = resstruct_tpt.results.Corrs(k);
            end
        end
    end
end
shapes = {'^','v','<','s'};
cmap = cool(length(modelnames));
for i = 1:length(modelnames)
    inds = (1:length(tptnames)) + (i-1)*length(tptnames);
    Rvec = Rmat_tpt(:,inds).';
    Rvecs(:,i) = Rvec(:);
end
R_all = NaN((size(Rvecs,1)+size(Rmat_all,1)),length(modelnames));
for i = 1:length(studynames)
    inds_tpt = (1:length(tptnames)) + ((i-1)*length(tptnames));
    ind_long = i*length(modelnames);
    inds_tpt_all = (1:length(tptnames)) + ((i-1)*length(modelnames));
    R_all(inds_tpt_all,:) = Rvecs(inds_tpt,:);
    R_all(ind_long,:) = Rmat_all(i,:);
end

xpos = 1:length(studynames); 
xposvec = repmat(xpos,(length(tptnames)+1),1);
xposvec = xposvec(:);
inds_cols = repmat(1:4,44,1);
inds_cols = inds_cols(:);
offsets = [-0.3,-0.1,0.1,0.3];
legstr = [];
figure('Units','inches','Position',[0 0 18 8]); hold on;
for i = 1:length(modelnames)
    for j = 1:length(xposvec)
        ind_shape = mod(j,(length(tptnames)+1));
        if ind_shape == 0
            ind_shape = length(tptnames)+1;
        end
        ind_cols1 = j + (i-1)*length(xposvec);
        ind_cols2 = inds_cols(ind_cols1);
        col = cmap(ind_cols2,:);
        if ind_shape == 4
            mfc = col;
            sz = 120;
        else
            mfc = 'none';
            sz = 60;
        end
        scatter((xposvec(j) + offsets(i)),R_all(j,i),sz,col,shapes{ind_shape},...
            'MarkerEdgeColor',col,'MarkerFaceColor',mfc);
    end
end
legmkrs = [{'o','o','o','o'},shapes];
legclrs = [cmap; zeros(length(shapes),3)];
legfcclrs = [cmap; ones(length(tptnames),3); zeros(1,3)];
if usepertpt_
    title('Per Timepoint')
    leglabs = {'Fit s','Ret.','Ant.','N.D.','t = 1','t = 2', 't = 3', 'Mean'};
    figstr = 'PerTpt';
else
    title('Longitudinal')
    leglabs = {'Fit s','Ret.','Ant.','N.D.','t = 1','t = 2', 't = 3', 'Long.'};
    figstr = 'Longitudinal';
end
leginds = 1:(length(tptnames)+1+length(modelnames));
for i = leginds
    s = scatter(0,0,80,legclrs(i,:),legmkrs{i},'MarkerFaceColor',legfcclrs(i,:));
    legstr = [legstr s];
end
for i = 1:length(studynames)
    plot([i,i]+0.5,[0,1],'LineStyle','-','LineWidth',1,'Color',[0.75 0.75 0.75])
end
xticks(xpos); xlim([min(xposvec)-0.5,max(xposvec)+0.5]); xticklabels(studylabels);
xlabel([]);
yplotmax = 1; yplotmin = 0;
ylim([yplotmin,yplotmax]); yticks(0:0.25:1);
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel('R');
legend(legstr,leglabs,'Location','northeast','NumColumns',2,'FontSize',20,'box','on')
set(gca,'FontSize',24,'FontName','Times','box','on','TickLength',[0 0.25]);

if savenclose_
    print([figdir_ filesep 'CorrComparePlot_' figstr],'-dtiffn','-r300'); close;
end
end