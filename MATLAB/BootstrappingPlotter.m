function BootstrappingPlotter(outputs,savenclose)

if nargin < 2
    savenclose = 0;
end

fldnames = fieldnames(outputs);
for k = 1:length(fldnames)
    tempstruct = outputs.(fldnames{k});
    subfldnames = fieldnames(tempstruct);
    if length(subfldnames) > 1
        niters = length(subfldnames)-1;
        param_fit = zeros(niters,length(tempstruct.Full.param_fit));
        for i = 1:niters
            param_fit(i,:) = tempstruct.(subfldnames{i}).param_fit;
        end
        if ismember('endm',fldnames) && ~logical(outputs.endm.Full.init.datapca_endm)
            ntypes = length(outputs.endm.Full.init.datalist_endm);
        else
            ntypes = 1;
        end
        if ntypes == 1
            xlabs = {'$\gamma$','$\alpha$','$\beta$','$s$','$a$','$b$','$p$'};
        else
            for i = 1:ntypes
                acell{i} = sprintf('$a_{%d}$',i);
                bcell{i} = sprintf('$b_{%d}$',i);
                pcell{i} = sprintf('$p_{%d}$',i);
            end
            xlabs = cat(2,{'$\gamma$','$\alpha$','$\beta$','$s$'},acell,bcell,pcell);
        end
        cmap = hsv(length(xlabs));
        xlabinds = 1:size(param_fit,2);
        if ~logical(tempstruct.Full.init.w_dir)
            xlabinds(4) = NaN;
        end
        if strcmp(fldnames{k},'ndm')
            xlabinds(5:end) = NaN;
        else
            xlabinds(5:(ntypes+4)) = NaN;
        end
        
        param_fit = param_fit(:,~isnan(xlabinds));
        cmap = cmap(~isnan(xlabinds),:);
        xlabs = xlabs(~isnan(xlabinds));

        rng(0);
        xpos = zeros((niters*size(param_fit,2)),1);
        g = xpos;
        xposscatter = @(x,y) 0.2 * (2*rand(length(x),1) - 1) + y;
        for j = 1:size(param_fit,2)
            curinds = (1:size(param_fit,1)) + (j-1)*size(param_fit,1);
            g(curinds) = j;
            xpos(curinds) = xposscatter(curinds,j);
        end

        figure('Position',[0 0 1000 500]); hold on; box on;
        gscatter(xpos,param_fit(:),g,cmap(1:size(param_fit,2),:),[],15,'off');
        b = boxplot(param_fit,'Colors',cmap(1:size(param_fit,2),:),'Symbol','');
        set(b,{'linew'},{3})
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        set(gca, 'XTick', 1:size(param_fit,2), 'XTickLabel', xlabs(1:size(param_fit,2)));
        % set(gca, 'YTick', [0.5 2 3.5])
        ylabel('');
        xlabel('');
        set(gca,'TickLength',[0 0])
        if strcmp(fldnames{k},'endm')
            title({sprintf('Nexis:sv Bootstrapping for Study %s:\nn = %d, resample rate = %.1f',...
                tempstruct.Full.init.study,niters,...
                tempstruct.Full.init.resample_rate_endm)},'FontName','Times');
        elseif strcmp(fldnames{k},'ndm')
            title({sprintf('Nexis:global Bootstrapping for Study %s:\nn = %d, resample rate = %.1f',...
                tempstruct.Full.init.study,niters,...
                tempstruct.Full.init.resample_rate)},'FontName','Times');
        end
        set(gca, 'FontSize', 30, 'LineWidth', 0.75);
        if savenclose
            if strcmp('ndm',fldnames{k})
                print(['Nexis_global_' tempstruct.Full.init.study '_scatterplot'],'-dtiffn');
            else
                namelist = outputs.endm.Full.init.datalist_endm;
                if isnumeric(namelist)
                    namelist = IndexName(namelist,...
                        outputs.endm.Full.init.datatype_endm);
                end
                if ~tempstruct.Full.init.datapca_endm
                    print(['Nexis_sv_' tempstruct.Full.init.study ...
                        '_' namelist{1} '_scatterplot'],'-dtiffn');
                else
                    print(['Nexis_sv_' tempstruct.Full.init.study ...
                        '_' namelist{1} '_PC1_scatterplot'],'-dtiffn');
                end
            end
            close;
        end
    end
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