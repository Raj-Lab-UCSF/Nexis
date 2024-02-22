function PerTimepointPlot_sbeta(outstruct,usebeta)
    
    rng(0);
    studynames = fieldnames(outstruct);
    studynames(ismember(studynames,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
    studylabels = cellfun(@(x)strrep(x,'_',' '),studynames,'UniformOutput',0);
    modelnames = fieldnames(outstruct.(studynames{1}));   

    tptnames = fieldnames(outstruct.(studynames{1}).(modelnames{1}));
    sbetavals = NaN(length(studynames),length(tptnames)); % 4 models
    tpts = sbetavals;
    for i = 1:size(sbetavals,1)
        % studyname = studynames{i};
        resstruct = outstruct.(studynames{i}).('fit_s');
        for k = 1:length(tptnames)
            resstruct_tpt = resstruct.(tptnames{k}).nexis_global.Full;
            if usebeta
                sbetavals(i,k) = resstruct_tpt.param_fit(3);
            else
                sbetavals(i,k) = resstruct_tpt.param_fit(4);
            end
            tpt = resstruct_tpt.time_stamps;
            % if strcmp(studyname(1:2),'DS')
            %     tpt = tpt/4; % weeks to months for Kaufman
            % end
            tpts(i,k) = tpt;
        end
    end
    
    cmap = hsv(length(studynames));
    shapes = {'o','s','d','^','v','<','>','p','h','+','x'};
    figure('Units','inches','Position',[0 0 11 10]); hold on;
    plothands = {};
    for i = 1:length(studynames)
        plot(tpts(i,:),sbetavals(i,:),'LineStyle','none','Color',cmap(i,:));
        s = scatter(tpts(i,:),sbetavals(i,:),75,shapes{i},...
            'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),...
            'MarkerFaceAlpha',0.3);
        plothands = [plothands,s];
    end
    xticks([0,3,6,9]); xlim([0,10]); xlabel('Time (Months)')
    if usebeta
        ylabel('\beta');
        yplotmax = 1.1*max(sbetavals(:)); yplotmin = 0; 
        ylim([yplotmin,yplotmax]); yticks([0,yplotmax/2,yplotmax]); 
        yticklabels({'0',num2str(yplotmax/2,'%.2f'),num2str(yplotmax,'%.2f')}); 
        title('Progression of diffusivity constant')
        loc = 'northeast';
    else
        ylabel('s');
        yplotmax = 1.1; yplotmin = -0.1; ylim([yplotmin,yplotmax]); 
        yticks([0,0.5,1]); yticklabels({'0','0.5','1'}); 
        title('Progression of directionality bias')
        loc = 'southeast';
    end
    legend(plothands,studylabels,'Location',loc,'NumColumns',3,'FontSize',20,'box','off');
    set(gca,'FontSize',24,'FontName','Times');

end