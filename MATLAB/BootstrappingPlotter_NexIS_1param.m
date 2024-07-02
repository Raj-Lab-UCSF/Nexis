function BootstrappingPlotter_NexIS_1param(outputs_bs,outputs_full,paramname,...
                                    fitsname,savenclose,figdir_)
% Need to fix for NexIS sv

rng(0);
studynames_ = fieldnames(outputs_bs);
studynames_plot = cellfun(@(x)strrep(x,'_',' '),studynames_,'UniformOutput', false);
paramnames = {'gamma','alpha','beta','s','b','p','R'};
parambool = ismember(paramnames,paramname);
param_fit = [];
param_fit_final = NaN(1,length(studynames_));
for m = 1:length(studynames_)
    outputs_bs_m = outputs_bs.(studynames_{m}).(fitsname);
    outputs_full_m = outputs_full.(studynames_{m}).(fitsname).nexis_global;
    fldnames = fieldnames(outputs_bs_m);
    for k = 1:length(fldnames)
        tempstruct = outputs_bs_m.(fldnames{k});
        subfldnames = fieldnames(tempstruct);
        % Following "full" stuff only works for NexIS global! needs fix for SV
        param_fit_full = outputs_full_m.Full.param_fit;
        R2_full = outputs_full_m.Full.results.lm_Rsquared_ord;
        param_fit_full = [param_fit_full sqrt(R2_full)];
        param_fit_final(m) = param_fit_full(parambool);
        if length(subfldnames) > 1
            niters = length(subfldnames)-1;
            niters_true = tempstruct.(subfldnames{k}).init.niters;
            param_fit_k = NaN(niters_true,length(tempstruct.Full.param_fit)+1); % takes into account failed sims
            for i = 1:niters
                param_fit_k(i,1:(end-1)) = tempstruct.(subfldnames{i}).param_fit;
                R2val_n = tempstruct.(subfldnames{i}).results.lm_Rsquared_ord;
                param_fit_k(i,end) = sqrt(R2val_n);
            end
            param_fit = [param_fit, param_fit_k(:,parambool)];
        end
    end
end

xlabs = studynames_plot;
cmap = hsv(length(xlabs));
xpos = zeros(niters_true,size(param_fit,2));
xposscatter = @(x,y) 0.2 * (2*rand(x,1) - 1) + y;
for j = 1:size(param_fit,2)
    xpos(:,j) = xposscatter(niters_true,j);
end
xlabinds = 1:length(studynames_plot);

figure('Position',[0 0 1500 500]); hold on; box on;
for j = 1:size(param_fit,2)
    scatter(xpos(:,j),param_fit(:,j),25,'MarkerEdgeColor',cmap(j,:));
end
% gscatter(xpos,param_fit(:),g,cmap,[],15,'off');
b = boxplot(param_fit,'Colors',cmap,'Symbol','');
set(b,{'linew'},{1.5})
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
    cmapind = 12 - j;
    patch(get(h(j),'XData'),get(h(j),'YData'),cmap(cmapind,:),'FaceAlpha',0.25);
end
b = boxplot(param_fit,'Colors',cmap,'Symbol','');
for j = 1:size(param_fit,2)
    plot([j-0.25, j+0.25],[param_fit_final(j),param_fit_final(j)],...
        'k-','LineWidth',4);
end
set(b,{'linew'},{1.5})

ax = gca;
ax.TickLabelInterpreter = 'latex';
set(gca, 'XTick', 1:size(param_fit,2), 'XTickLabel', xlabs(xlabinds));
if strcmp(paramname,'s')
    yticks([0 0.5 1]); ylim([-0.1 1.1]);
else
    ymax = max(param_fit(:)); ymin = min(param_fit(:));
    if ymin < 0
        yplotmin = 1.1*ymin;
    else
        yplotmin = 0.9*ymin;
    end
    yplotmax = 1.1*ymax; 
    yticks([ymin mean([ymin,ymax]) ymax]); ylim([yplotmin,yplotmax]); 
    ytickformat('%.2f');
end
xlim([0.5,length(studynames_)+0.5])
if ismember(paramname,{'gamma','alpha','beta'})
    paramnamestr = ['\' paramname];
else
    paramnamestr = paramname;
end
ylabel(paramnamestr,'FontWeight','normal','FontName','Times');
xlabel('');
set(gca,'TickLength',[0 0])
if strcmp(fldnames{k},'nexis_sv') % fix later
    title({sprintf('Nexis:sv Bootstrapping for Study %s:\nn = %d, resample rate = %.1f',...
        tempstruct.Full.init.study,niters,...
        tempstruct.Full.init.resample_rate_endm)},'FontName','Times');
elseif strcmp(fldnames{k},'nexis_global')
    title({sprintf('Nexis:global Bootstrapping for %s, resample rate = %.1f',...
        paramnamestr,...
        tempstruct.Full.init.resample_rate)},'FontName','Times');
end
set(gca, 'FontSize', 30, 'LineWidth', 0.75);
if savenclose
    if strcmp('nexis_global',fldnames{k})
        print([figdir_ filesep 'Nexis_global_' paramname '_boxplot'],'-dtiffn');
    else % fix for sv
        namelist = outputs_bs_m.endm.Full.init.datalist_endm;
        if isnumeric(namelist)
            namelist = IndexName(namelist,...
                outputs_bs_m.endm.Full.init.datatype_endm);
        end
        if ~tempstruct.Full.init.datapca_endm
            print([figdir filesep 'Nexis_sv_' tempstruct.Full.init.study ...
                '_' namelist{1} '_parameterplot'],'-dtiffn');
        else
            print([figdir filesep 'Nexis_sv_' tempstruct.Full.init.study ...
                '_' namelist{1} '_PC1_parameterplot'],'-dtiffn');
        end
    end
    close;
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