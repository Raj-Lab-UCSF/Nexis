function CorrelationPlotter_single(predicted,data,tpts,study,sv_types,...
                                    geneorct,ispca,savenclose)

if nargin < 8
    savenclose = 0;
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

notnaninds = ~isnan(data(:,1));
data = data(notnaninds,:);
predicted = predicted(notnaninds,:);
plot_max_data = 1.1*max(max(data));
plot_max_pred = 1.1*max(max(predicted));

num_plots = size(data,2);
if isempty(sv_types)
    flnamestr = 'Nexis_global';
    fig = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8],...
    'Units','inch','Position',[0 0 5*num_plots 4.85]);
    sgtitle(['Nexis:global - Study: ' study],'FontName','Times',...
        'FontSize',32,'FontWeight','bold');
else
    flnamestr = 'Nexis_sv';
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
    sgtitle({['Nexis:sv - Study: ' study ', '], namestr},'FontName','Times',...
        'FontSize',32,'FontWeight','bold');
end

for j = 1:num_plots
    subplot(1,num_plots,j);
    scatter(data(:,j),predicted(:,j),30,'o','MarkerFaceColor','k',...
        'MarkerEdgeColor','k')
    box on;
    axis([0 plot_max_data 0 plot_max_pred])
    h1 = lsline;
    h1.Color = 'r'; h1.LineWidth = 3;
    legend(h1,sprintf('R = %.2f',corr(data(:,j),predicted(:,j),...
        'rows','complete')),'Location','northwest')
    xticks([0,plot_max_data/2,plot_max_data]);
    xtickformat('%,.1f');
    yticks([0,plot_max_pred/2,plot_max_pred]);
    ytickformat('%,.1f');
    if j == 1
        ylabel('Modeled Pathology');
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
xlab = xlabel(h,'Observed Pathology','FontSize',20,'FontName','Times');
xlab.Position(2) = 0;
if savenclose
    if isempty(sv_types)
        print([flnamestr '_' study '_scatterplot'],'-dtiffn');
    else
        if ~logical(ispca)
            print([flnamestr '_' study ...
                '_' sv_types{1} '_scatterplot'],'-dtiffn');
        else
            print([flnamestr '_' tempstruct.Full.init.study ...
                '_' sv_types{1} '_PC1_scatterplot'],'-dtiffn');
        end
    end
    close;
end
end