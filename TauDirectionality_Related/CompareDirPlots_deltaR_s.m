function [Rmat, svals, tstats_struct] = CompareDirPlots_deltaR_s(outstruct,...
                    pertimepoint,savenclose_,figdir_)

rng(0);
studynames = fieldnames(outstruct);
studynames(ismember(studynames,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
modelnames = fieldnames(outstruct.(studynames{1}));
fisher_rtoz = @(r) 0.5*(log(1+r) - log(1-r));
tstats_struct = struct;

if ~pertimepoint
    Rmat = NaN(length(studynames),length(modelnames)); % 4 models
    svals = NaN(length(studynames),1);
    for i = 1:size(Rmat,1)
        for j = 1:size(Rmat,2)
            resstruct = outstruct.(studynames{i}).(modelnames{j}).nexis_global.Full;
            Rmat(i,j) = (resstruct.results.lm_Rsquared_ord)^(0.5);
            if j == 1
                svals(i) = resstruct.param_fit(4);
            end
        end
    end
    retind = find(ismember(modelnames,'ret'));
    antind = find(ismember(modelnames,'ant'));
    Rdiffs = Rmat(:,retind) - Rmat(:,antind); %#ok<FNDSB>
    zdiffs = fisher_rtoz(Rdiffs);
    svals_trans = 2*(svals - 0.5); % Put on [-1,1]
    [~,pvalR,~,statsR] = ttest(zdiffs);
    [~,pvals,~,statss] = ttest(svals_trans);
    tstats_struct.Rdiff.pval = pvalR; tstats_struct.Rdiff.tstat = statsR.tstat;
    tstats_struct.s.pval = pvals; tstats_struct.s.tstat = statss.tstat;

    cmap_violinplot = [[1 0 0]; [0 0 1]];
    xposscatter = @(y) 0.2 * (2*rand - 1) + y;
    xpos_R = NaN(length(Rdiffs),1); 
    % gvec_R = xpos_R;
    xpos_s = xpos_R; 
    % gvec_s = gvec_R;
    for i = 1:length(xpos_R)
        xpos_R(i) = xposscatter(1);
        % gvec_R(i) = 1;
    end
    for i = 1:length(xpos_s)
        xpos_s(i) = xposscatter(1);
        % gvec_s(i) = 2;
    end
    allvals = [Rdiffs, svals];
    % gvec = [gvec_R, gvec_s];

    figure('Units','inches','Position',[0 0 3.5 5]); hold on;
    b = boxchart(ones(size(allvals,1),1),allvals(:,1),'BoxFaceColor',...
        cmap_violinplot(1,:),'MarkerStyle','none');
    set(b,{'linew'},{2});
    % violin(Rdiffs,'facecolor',cmap_violinplot(1,:),'medc',[]);
    scatter(xpos_R,Rdiffs,[],cmap_violinplot(1,:),'filled');
    plot([0.5,1.5],[0 0],'k:','LineWidth',1)
    % hLegend = findobj(gcf, 'Type', 'Legend'); hLegend.Visible = 'off';
    xticks(1); xlim([0.5,1.5]); xticklabels({'$\Delta R_{dir}$'});    
    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex';
    yplotmax = 0.35; yplotmin = -0.05;
    ylim([yplotmin,yplotmax]); yticks([0,0.1,0.2,0.3]);
    yticklabels({'0','0.1','0.2','0.3'});
    title('Longitudinal');
    set(gca,'FontSize',20,'FontName','Times');
    if savenclose_
        print([figdir_ filesep 'DeltaRViolin_Longitudinal'],'-dtiffn','-r300'); close;
    end

    figure('Units','inches','Position',[0 0 3.5 5]); hold on;
    b = boxchart(ones(size(allvals,1),1),allvals(:,2),'BoxFaceColor',...
        cmap_violinplot(2,:),'MarkerStyle','none');
    set(b,{'linew'},{2});
    % violin(svals,'facecolor',cmap_violinplot(2,:),'medc',[]);
    scatter(xpos_s,svals,[],cmap_violinplot(2,:),'filled');
    plot([0.5,1.5],[0.5 0.5],'k:','LineWidth',1)
    % hLegend = findobj(gcf, 'Type', 'Legend'); hLegend.Visible = 'off';
    xticks(1); xlim([0.5,1.5]); xticklabels({'$s$'});    
    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex';
    yplotmax = 0.9; yplotmin = 0.3;
    ylim([yplotmin,yplotmax]); yticks([0.4,0.6,0.8]);
    yticklabels({'0.4','0.6','0.8'});
    title('Longitudinal');
    set(gca,'FontSize',20,'FontName','Times');
    if savenclose_
        print([figdir_ filesep 'sViolin_Longitudinal'],'-dtiffn','-r300'); close;
    end

else
    tptnames = fieldnames(outstruct.(studynames{1}).(modelnames{1}));
    Rmat = NaN(length(studynames),length(modelnames)*length(tptnames)); % 4 models
    svals = NaN(length(studynames),length(tptnames));
    for i = 1:size(Rmat,1)
        for j = 1:length(modelnames)
            resstruct = outstruct.(studynames{i}).(modelnames{j});
            for k = 1:length(tptnames)
                resstruct_tpt = resstruct.(tptnames{k}).nexis_global.Full;
                Rind = k + length(tptnames)*(j-1);
                Rmat(i,Rind) = (resstruct_tpt.results.lm_Rsquared_ord)^(0.5);
                if j == 1
                    svals(i,k) = resstruct_tpt.param_fit(4);
                end
            end
        end
    end

    modelnames_ind = repmat(modelnames.',length(tptnames),1);
    modelnames_ind = modelnames_ind(:);
    retinds = find(ismember(modelnames_ind,'ret'));
    antinds = find(ismember(modelnames_ind,'ant'));
    Rdiffs = Rmat(:,retinds) - Rmat(:,antinds); %#ok<FNDSB>
    Rdiffs = Rdiffs(:);
    zdiffs = fisher_rtoz(Rdiffs);
    svals = svals(:);
    svals_trans = 2*(svals - 0.5); % Put on [-1,1]
    [~,pvalR,~,statsR] = ttest(zdiffs(:));
    [~,pvals,~,statss] = ttest(svals_trans(:));
    tstats_struct.Rdiff.pval = pvalR; tstats_struct.Rdiff.tstat = statsR.tstat;
    tstats_struct.s.pval = pvals; tstats_struct.s.tstat = statss.tstat;

    cmap_violinplot = [[1 0 0]; [0 0 1]];
    xposscatter = @(y) 0.2 * (2*rand - 1) + y;
    xpos_R = NaN(length(Rdiffs),1); 
    % gvec_R = xpos_R;
    xpos_s = xpos_R; 
    % gvec_s = gvec_R;
    for i = 1:length(xpos_R)
        xpos_R(i) = xposscatter(1);
        % gvec_R(i) = 1;
    end
    for i = 1:length(xpos_s)
        xpos_s(i) = xposscatter(1);
        % gvec_s(i) = 2;
    end
    allvals = [Rdiffs, svals];
    % gvec = [gvec_R, gvec_s];

    figure('Units','inches','Position',[0 0 3.5 5]); hold on;
    b = boxchart(ones(size(allvals,1),1),allvals(:,1),'BoxFaceColor',cmap_violinplot(1,:));
    set(b,{'linew'},{2});
    % violin(Rdiffs,'facecolor',cmap_violinplot(1,:),'medc',[]);
    scatter(xpos_R,Rdiffs,[],cmap_violinplot(1,:),'filled');
    plot([0.5,1.5],[0 0],'k:','LineWidth',1)
    % hLegend = findobj(gcf, 'Type', 'Legend'); hLegend.Visible = 'off';
    xticks(1); xlim([0.5,1.5]); xticklabels({'$\Delta R_{dir}$'});    
    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex';
    yplotmax = 0.4; yplotmin = -0.15;
    ylim([yplotmin,yplotmax]); yticks([-0.1,0.1,0.3]);
    yticklabels({'-0.1','0.1','0.3'});
    title('Per Timepoint');
    set(gca,'FontSize',20,'FontName','Times');
    if savenclose_
        print([figdir_ filesep 'DeltaRViolin_PerTpt'],'-dtiffn','-r300'); close;
    end

    figure('Units','inches','Position',[0 0 3.5 5]); hold on;
    b = boxchart(ones(size(allvals,1),1),allvals(:,2),'BoxFaceColor',...
        cmap_violinplot(2,:),'MarkerStyle','none');
    set(b,{'linew'},{2});
    % violin(svals,'facecolor',cmap_violinplot(2,:),'medc',[]);
    scatter(xpos_s,svals,[],cmap_violinplot(2,:),'filled');
    plot([0.5,1.5],[0.5 0.5],'k:','LineWidth',1)
    % hLegend = findobj(gcf, 'Type', 'Legend'); hLegend.Visible = 'off';
    xticks(1); xlim([0.5,1.5]); xticklabels({'$s$'});    
    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex';
    yplotmax = 1.1; yplotmin = -0.1;
    ylim([yplotmin,yplotmax]); yticks([0,0.5,1]);
    yticklabels({'0','0.5','1'});
    title('Per Timepoint');
    set(gca,'FontSize',20,'FontName','Times');
    if savenclose_
        print([figdir_ filesep 'sViolin_PerTpt'],'-dtiffn','-r300'); close;
    end

end
end