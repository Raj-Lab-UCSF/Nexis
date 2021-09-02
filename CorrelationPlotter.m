function CorrelationPlotter(outputs)

fldnames = fieldnames(outputs);

% Rvals across time points for optimal model
for k = 1:length(fldnames)
    tempstruct = outputs.(fldnames{k});
    data = tempstruct.Full.data;
    predicted = tempstruct.Full.predicted;
    notnaninds = ~isnan(data(:,1));
    data = data(notnaninds,:);
    predicted = predicted(notnaninds,:);
    plot_max_data = 1.1*max(max(data));
    plot_max_pred = 1.1*max(max(predicted));
    num_plots = size(data,2);
    figure('Name',[tempstruct.Full.init.study '_' fldnames{k}],'Units','inch','Position',[0 0 20 7]);
    if strcmp('ndm',fldnames{k})
        sgtitle(['NDM - Study: ' tempstruct.Full.init.study],'FontName','Times','FontSize',32,'FontWeight','bold');
    else
        datsettype = tempstruct.Full.init.datatype_endm;
        namelist = tempstruct.Full.init.datalist_endm;
        if isnumeric(namelist)
            namelist = IndexName(namelist,datsettype);
        end
        if strcmp(datsettype,'gene')
            if length(namelist) > 1
                dattypestr = 'Genes: ';
            else
                dattypestr = 'Gene: ';
            end
        else
            if length(namelist) > 1
                dattypestr = 'Cell Types: ';
            else
                dattypestr = 'Cell Type: ';
            end            
        end
        namestr = [dattypestr namelist{1}];
        if length(namelist) > 1
            for i = 2:length(namelist)
                namestr = [namestr ', ' namelist{i}];
            end
        end
        if logical(tempstruct.Full.init.datapca_endm)
            namestr = [namestr ' (1st PC)'];
        end
        sgtitle(['eNDM - Study: ' tempstruct.Full.init.study ', ' namestr],'FontName','Times','FontSize',32,'FontWeight','bold');
    end
    for j = 1:num_plots
        subplot(1,num_plots,j)
        scatter(data(:,j),predicted(:,j),30,'o','MarkerFaceColor','k','MarkerEdgeColor','k')
        axis([0 plot_max_data 0 plot_max_pred])
        h1 = lsline;
        h1.Color = 'r'; h1.LineWidth = 3;
        legend(h1,sprintf('R = %.2f',corr(data(:,j),predicted(:,j),'rows','complete')),'Location','northwest')
        xlabel('Observed Pathology')
        ylabel('Modeled Pathology')
        set(gca,'FontSize',18,'FontName','Times')
        title(['Time = ' num2str(tempstruct.Full.time_stamps(j))]);
    end  
    clear tempstruct
end

if logical(outputs.ndm.Full.init.bootstrapping) || logical(outputs.endm.Full.init.bootstrapping_endm)
    if length(fldnames) > 1
        corrscell = cell(1,length(fldnames));
    else
        corrscell = cell(1,2);
    end
    for k = 1:length(fldnames)
        tempstruct = outputs.(fldnames{k});
        subfldnames = fieldnames(tempstruct);
        if length(subfldnames) > 1
            niters = length(subfldnames)-1;
            corrs = zeros(niters,1);
            for i = 1:niters
                corrs(i) = tempstruct.(subfldnames{i}).results.lm_Rsquared_adj;
            end
            corrscell{k} = corrs;
        end
        clear tempstruct corrs
    end
    if ismember('endm',fldnames)
        datsettype = outputs.endm.Full.init.datatype_endm;
        namelist = outputs.endm.Full.init.datalist_endm;
        if isnumeric(namelist)
            namelist = IndexName(namelist,datsettype);
        end
        if strcmp(datsettype,'gene')
            if length(namelist) > 1
                dattypestr = 'Genes: ';
            else
                dattypestr = 'Gene: ';
            end
        else
            if length(namelist) > 1
                dattypestr = 'Cell Types: ';
            else
                dattypestr = 'Cell Type: ';
            end            
        end
        namestr = [dattypestr namelist{1}];
        if length(namelist) > 1
            for i = 2:length(namelist)
                namestr = [namestr ', ' namelist{i}];
            end
        end
        if logical(outputs.endm.Full.init.datapca_endm)
            namestr = [namestr ' (1st PC)'];
        end
        xlabs = {'NDM',['\begin{tabular}{c} eNDM \\' namestr '\end{tabular}']};
    else
        xlabs = {'NDM','eNDM'};
    end
    cmap = [[1 0 0]; [0 0 1]];
    xlabinds = 1:length(xlabs);
    for k = 1:length(corrscell)
        if isempty(corrscell{k})
            xlabinds(k) = NaN;
        end
    end
    
    corrscell = corrscell(~isnan(xlabinds));
    cmap = cmap(~isnan(xlabinds),:);
    xlabs = xlabs(~isnan(xlabinds));

    rng(0);
    corrsfull = zeros(1,length(corrscell));
    if length(corrscell) > 1
        corrs = [];
        for k = 1:length(corrscell)
            corrs = [corrs; corrscell{k}];
            corrsfull(k) = outputs.(fldnames{k}).Full.results.lm_Rsquared_adj;
        end
        xpos = zeros(length(corrs),1);
    else
        corrs = corrscell{1};
        corrsfull = outputs.(fldnames{xlabinds(~isnan(xlabinds))}).Full.results.lm_Rsquared_adj;
        xpos = zeros(length(corrs),1);
    end
    g = xpos;
    xposscatter = @(x,y) 0.2 * (2*rand(length(x),1) - 1) + y;
    for j = 1:length(corrscell)
        if j == 1
            curinds = (1:length(corrscell{j}));
        else
            curinds = (1:length(corrscell{j})) + (j-1)*length(corrscell{j-1});
        end
        g(curinds) = j;
        xpos(curinds) = xposscatter(curinds,j);
    end
    figure('Position',[0 0 1000 500]); hold on; box on;
    gscatter(xpos,corrs,g,cmap,[],15,'off');
    b = boxplot(corrs,g,'Colors',cmap,'Symbol','');
    for k = 1:length(corrscell)
       h1 = plot([k-0.2, k+0.2],[corrsfull(k),corrsfull(k)],'m--','LineWidth',3); 
    end
    set(b,{'linew'},{3})
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    set(gca, 'XTick', 1:length(corrscell), 'XTickLabel', xlabs);
    ylim([0.9*min([min(corrs),min(corrsfull)]),1.1*max([max(corrs),max(corrsfull)])]);
    legend(h1,'Full Dataset with Optimal Parameters','Location','northwest');
    % set(gca, 'YTick', [0.5 2 3.5])
    ylabel('R^{2}');
    xlabel('');
    set(gca,'TickLength',[0 0])
    title(sprintf('Bootstrapped Performance for Study %s',outputs.ndm.Full.init.study),'FontName','Times');
    set(gca, 'FontSize', 24, 'LineWidth', 0.75,'FontName','Times');
end

    function names = IndexName(indices,dattypeendm)
        if strcmp(dattypeendm,'gene')
            load([cd filesep 'raw_data_mouse' filesep 'gene_names_trans.mat'],'gene_names_trans');
            namescell = gene_names_trans;
        elseif strcmp(dattypeendm,'ct_tasic')
            load([cd filesep 'raw_data_mouse' filesep 'classkey_tasic.mat'],'classkey_tasic');
            namescell = classkey_tasic;
        elseif strcmp(dattypeendm,'ct_zeisel')
            load([cd filesep 'raw_data_mouse' filesep 'classkey_zeisel.mat'],'classkey_zeisel');
            namescell = classkey_zeisel;
        end
        
        names = namescell(indices);
    end



end