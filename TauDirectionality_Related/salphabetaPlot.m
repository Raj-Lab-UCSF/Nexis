function [alpha_mat,beta_mat,s_mat] = salphabetaPlot(outstruct,comparetype_,savenclose_,figdir_)
     
studynames = fieldnames(outstruct);
studynames(ismember(studynames,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
studylabels = cellfun(@(x)strrep(x,'_',' '),studynames,'UniformOutput',0);
modelnames = fieldnames(outstruct.(studynames{1}));

tptnames = fieldnames(outstruct.(studynames{1}).(modelnames{1}));
s_mat = NaN(length(studynames),length(tptnames)); % 4 models
alpha_mat = s_mat;
beta_mat = s_mat;
for i = 1:size(s_mat,1)
    resstruct = outstruct.(studynames{i}).('fit_s');
    for k = 1:length(tptnames)
        resstruct_tpt = resstruct.(tptnames{k}).nexis_global.Full;
        alpha_mat(i,k) = resstruct_tpt.param_fit(2);
        beta_mat(i,k) = resstruct_tpt.param_fit(3);
        s_mat(i,k) = resstruct_tpt.param_fit(4);
    end
end

cmap = hsv(length(studynames));
shapes = {'o','s','d','^','v','<','>','p','h','+','x'};
figure('Units','inches','Position',[0 0 11 10]); hold on;
plothands = {};
switch comparetype_
    case 'alpha_s'
        for i = 1:length(studynames)
            s = scatter(alpha_mat(i,:),s_mat(i,:),75,shapes{i},...
                'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),...
                'MarkerFaceAlpha',0.3);
            plothands = [plothands,s];
        end
        lm = fitlm(alpha_mat(:),s_mat(:));
        x_lm = linspace(0,1.5*max(alpha_mat(:)),100).'; 
        [y_lm, y_ci] = predict(lm, x_lm);
        plot(x_lm,y_ci(:,1),'k:'); plot(x_lm,y_ci(:,2),'k:'); 
        fill([x_lm; flipud(x_lm)],[y_ci(:,1); flipud(y_ci(:,2))],[1 0 0.25],...
            'EdgeColor','none','FaceAlpha',0.15);
        plot(x_lm,y_lm,'k','LineWidth',2);
        Rptext = sprintf(['R^2 = %.2f,' newline 'p = %.1d'],...
                    lm.Rsquared.Adjusted,lm.ModelFitVsNullModel.Pvalue);
        
        xplotmax = 1.1*max(alpha_mat(:)); xplotmin = 0;
        xticks([0,max(alpha_mat(:))/2,max(alpha_mat(:))]);
        xticklabels({'0',num2str(max(alpha_mat(:))/2,'%.2f'),num2str(max(alpha_mat(:)),'%.2f')}); 
        xlim([xplotmin,xplotmax]);
        yplotmax = 1.1; yplotmin = -0.1; ylim([yplotmin,yplotmax]); 
        yticks([0,0.5,1]); yticklabels({'0','0.5','1'}); 
        loc = 'southeast';
        text(0.75,0.95,Rptext,'FontSize',20,'FontName','Times','Units','normalized');
        ylabel('s'); xlabel('\alpha'); title('Bias vs. Accumulation Parameters');
        legend(plothands,studylabels,'Location',loc,'NumColumns',3,'FontSize',20);
        set(gca,'FontSize',24,'FontName','Times');
        
        if savenclose_
            print([figdir_ filesep 'salpha_plot'],'-dtiffn','-r300'); close;
        end
        
    case 'beta_s'
        for i = 1:length(studynames)
            s = scatter(beta_mat(i,:),s_mat(i,:),75,shapes{i},...
                'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),...
                'MarkerFaceAlpha',0.3);
            plothands = [plothands,s];
        end
        lm = fitlm(beta_mat(:),s_mat(:));
        x_lm = linspace(0,1.5*max(beta_mat(:)),100).'; 
        [y_lm, y_ci] = predict(lm, x_lm);
        plot(x_lm,y_ci(:,1),'k:'); plot(x_lm,y_ci(:,2),'k:'); 
        fill([x_lm; flipud(x_lm)],[y_ci(:,1); flipud(y_ci(:,2))],[1 0 0.25],...
            'EdgeColor','none','FaceAlpha',0.15);
        plot(x_lm,y_lm,'k','LineWidth',2);
        Rptext = sprintf(['R^2 = %.2f,' newline 'p = %.1d'],...
                    lm.Rsquared.Adjusted,lm.ModelFitVsNullModel.Pvalue);
        
        xplotmax = 1.1*max(beta_mat(:)); xplotmin = 0;
        xticks([0,max(beta_mat(:))/2,max(beta_mat(:))]);
        xticklabels({'0',num2str(max(beta_mat(:))/2,'%.1f'),num2str(max(beta_mat(:)),'%.1f')}); 
        xlim([xplotmin,xplotmax]);
        yplotmax = 1.1; yplotmin = -0.1; ylim([yplotmin,yplotmax]); 
        yticks([0,0.5,1]); yticklabels({'0','0.5','1'}); 
        loc = 'southeast';
        text(0.75,0.95,Rptext,'FontSize',20,'FontName','Times','Units','normalized');
        
        ylabel('s'); xlabel('\beta'); title('Bias vs. Spread Parameters');
        legend(plothands,studylabels,'Location',loc,'NumColumns',3,'FontSize',20);
        set(gca,'FontSize',24,'FontName','Times');
        
        if savenclose_
            print([figdir_ filesep 'sbeta_plot'],'-dtiffn','-r300'); close;
        end

    case 'alpha_beta'
        for i = 1:length(studynames)
            s = scatter(alpha_mat(i,:),beta_mat(i,:),75,shapes{i},...
                'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),...
                'MarkerFaceAlpha',0.3);
            plothands = [plothands,s];
        end
        lm = fitlm(alpha_mat(:),beta_mat(:));
        x_lm = linspace(0,1.5*max(alpha_mat(:)),100).'; 
        [y_lm, y_ci] = predict(lm, x_lm);
        plot(x_lm,y_ci(:,1),'k:'); plot(x_lm,y_ci(:,2),'k:'); 
        fill([x_lm; flipud(x_lm)],[y_ci(:,1); flipud(y_ci(:,2))],[1 0 0.25],...
            'EdgeColor','none','FaceAlpha',0.15);
        plot(x_lm,y_lm,'k','LineWidth',2);
        Rptext = sprintf(['R^2 = %.2f,' newline 'p = %.1d'],...
                    lm.Rsquared.Adjusted,lm.ModelFitVsNullModel.Pvalue);
        
        xplotmax = 1.1*max(alpha_mat(:)); xplotmin = 0;
        xticks([0,max(alpha_mat(:))/2,max(alpha_mat(:))]);
        xticklabels({'0',num2str(max(alpha_mat(:))/2,'%.1f'),num2str(max(alpha_mat(:)),'%.1f')}); 
        xlim([xplotmin,xplotmax]);
        yplotmax = 1.1*max(beta_mat(:)); yplotmin = 0;
        yticks([0,max(beta_mat(:))/2,max(beta_mat(:))]);
        yticklabels({'0',num2str(max(beta_mat(:))/2,'%.1f'),num2str(max(beta_mat(:)),'%.1f')}); 
        ylim([yplotmin,yplotmax]);
        loc = 'northwest';
        text(0.75,0.1,Rptext,'FontSize',20,'FontName','Times','Units','normalized');
        
        ylabel('\beta'); xlabel('\alpha'); title('Spread vs. Accumulation Parameters');
        legend(plothands,studylabels,'Location',loc,'NumColumns',3,'FontSize',19);
        set(gca,'FontSize',24,'FontName','Times');
        
        if savenclose_
            print([figdir_ filesep 'alpha_beta_plot'],'-dtiffn','-r300'); close;
        end
end
end