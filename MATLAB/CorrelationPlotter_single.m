function CorrelationPlotter_single(predicted,data,tpts,study,sv_types,...
                                    geneorct,ispca,wdir,color,shape,exclinit,...
                                    savenclose)

if nargin < 12
    savenclose = 0;
    if nargin < 11
        exclinit = 0;
        if nargin < 10
            shape = 'o';
            if nargin < 9
                color = [1 0 0];
                if nargin < 8
                    wdir = 0;
                    if nargin < 7
                        ispca = 0;
                        if nargin < 6
                            geneorct = 'gene';
                            if nargin < 5
                                sv_types = {};
                            end
                        end
                    end
                end
            end
        end
    end
end

notnaninds = ~isnan(data(:,1));
data = data(notnaninds,:);
predicted = predicted(notnaninds,:);
if exclinit
    data = data(:,2:end);
    predicted = predicted(:,2:end);
    tpts = tpts(2:end);
end
plot_max_data = max(max(data));
plot_max_pred = max(max(predicted));

num_plots = size(data,2);
if isempty(sv_types)
    if logical(wdir)       
        flnamestr = 'Nexis_global_wdir';
        fig = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8],...
        'Units','inch','Position',[0 0 5*num_plots 4.85]);
        sgtitle(['Nexis:global with Directionality - Study: ' study],'FontName','Times',...
            'FontSize',32,'FontWeight','bold');
    else
        flnamestr = 'Nexis_global_nodir';
        fig = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8],...
        'Units','inch','Position',[0 0 5*num_plots 4.85]);
        sgtitle(['Nexis:global without Directionality - Study: ' study],'FontName','Times',...
            'FontSize',32,'FontWeight','bold');
    end
else
    if logical(wdir)
        flnamestr = 'Nexis_sv_wdir';
        fig = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8],...
        'Units','inch','Position',[0 0 5*num_plots 5.45]);
        if strcmp(geneorct,'gene')
            if length(sv_types) > 1
                dattypestr = 'Genes: ';
            else
                dattypestr = 'Gene: ';
            end
        else
            if length(sv_types) > 1
                dattypestr = 'Cell Types: ';
            else
                dattypestr = 'Cell Type: ';
            end            
        end
        namestr = [dattypestr sv_types{1}];
        if length(sv_types) > 1
            for i = 2:length(sv_types)
                namestr = [namestr ', ' sv_types{i}];
            end
        end
        if logical(ispca)
            namestr = [namestr ' (1st PC)'];
        end
        sgtitle({['Nexis:sv with Directionalty - Study: ' study ', '], namestr},'FontName','Times',...
            'FontSize',32,'FontWeight','bold');
    else
        flnamestr = 'Nexis_sv_nodir';
        fig = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8],...
        'Units','inch','Position',[0 0 5*num_plots 5.45]);
        if strcmp(geneorct,'gene')
            if length(sv_types) > 1
                dattypestr = 'Genes: ';
            else
                dattypestr = 'Gene: ';
            end
        else
            if length(sv_types) > 1
                dattypestr = 'Cell Types: ';
            else
                dattypestr = 'Cell Type: ';
            end            
        end
        namestr = [dattypestr sv_types{1}];
        if length(sv_types) > 1
            for i = 2:length(sv_types)
                namestr = [namestr ', ' sv_types{i}];
            end
        end
        if logical(ispca)
            namestr = [namestr ' (1st PC)'];
        end
        sgtitle({['Nexis:sv without Directionality - Study: ' study ', '], namestr},'FontName','Times',...
            'FontSize',32,'FontWeight','bold');
    end
end

for j = 1:num_plots
    subplot(1,num_plots,j);
    if logical(wdir)
        scatter(predicted(:,j),data(:,j),50,shape,'MarkerFaceColor',color,...
            'MarkerEdgeColor',color,'LineWidth',1)
    else
        scatter(predicted(:,j),data(:,j),50,shape,'MarkerFaceColor','w',...
            'MarkerEdgeColor',color,'LineWidth',1)
    end
    box on; hold on;
    p = polyfit(predicted(:,j),data(:,j),1);
    lsx = linspace(0,1.1*plot_max_pred);
    lsy = polyval(p,lsx);
    plot(lsx,lsy,'--','LineWidth',3);
    h1 = lsline;
    h1.Color = 'k'; h1.LineWidth = 3; h1.LineStyle = '--';
    legend(h1,sprintf('R = %.2f',corr(data(:,j),predicted(:,j),...
        'rows','complete')),'Location','northwest')
    xticks([plot_max_pred/2,plot_max_pred]);
    xtickformat('%.2f');
    xlim([0 plot_max_pred*1.1]);
    yticks([plot_max_data/2,plot_max_data]);
    ytickformat('%.2f');
    ylim([0 plot_max_data*1.1]);
    text(-0.02*plot_max_pred*1.1,-0.08*plot_max_data*1.1,'0','FontSize',20,'FontName','Times')
    if j == 1
        ylabel('Observed Pathology');
        text(-0.08*plot_max_pred*1.1,0,'0','FontSize',20,'FontName','Times')
    else
        yticklabels({});
    end
    set(gca,'FontSize',20,'FontName','Times')
    title(['Time = ' num2str(tpts(j))]);
    axis('square')
end  
h = axes(fig,'visible','off');
h.Title.Visible='on';
h.XLabel.Visible='on';
xlab = xlabel(h,'Modeled Pathology','FontSize',20,'FontName','Times');
xlab.Position(2) = 0;
if savenclose
    if isempty(sv_types)
        print([flnamestr '_' study '_scatterplot'],'-dtiffn','-r300');
    else
        if ~logical(ispca)
            print([flnamestr '_' study ...
                '_' sv_types{1} '_scatterplot'],'-dtiffn','-r300');
        else
            print([flnamestr '_' tempstruct.Full.init.study ...
                '_' sv_types{1} '_PC1_scatterplot'],'-dtiffn','-r300');
        end
    end
    close;
end
end