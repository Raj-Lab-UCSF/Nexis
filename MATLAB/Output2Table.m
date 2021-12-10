function summarytable = Output2Table(outputs,writeout,filename,filepath)

fldnames = fieldnames(outputs);
if nargin < 4
    filepath = cd;
    if nargin < 3
        datestr = date;
        if ~ismember('endm',fldnames)
            filename = ['summary_Nexis_mouse_' outputs.ndm.Full.init.study...
                '_global_' datestr];             
        else
            typename = outputs.endm.Full.init.datalist_endm(1);
            if isnumeric(typename)
                typename = IndexName(typename,outputs.endm.Full.init.datatype_endm);
            end
            filename = ['summary_Nexis_mouse_' outputs.endm.Full.init.study...
                '_' typename{1} '_' datestr];  
        end
        if nargin < 2
            writeout = 0;
        end
    end
end
            
rownames = fldnames.';
columnnames = cell(1,1); vartypes = columnnames;
if ismember('endm',fldnames)
    if strcmp(outputs.endm.Full.init.datatype_endm,'gene')
        colstr1 = 'Gene';
    else
        colstr1 = 'Cell Type';
    end
    typenames = outputs.endm.Full.init.datalist_endm;
    for i = 1:length(typenames)
        columnnames{i} = sprintf([colstr1 ' %d'],i);
        vartypes{i} = 'string';
    end
else
   columnnames{1} = 'Gene 1'; 
   vartypes{1} = 'string';
end
columnnames{end+1} = 'Uses PCA'; vartypes{end+1} = 'string';    
columnnames{end+1} = 'Cost Function'; vartypes{end+1} = 'string';
columnnames{end+1} = 'gamma (Mean)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'gamma (95% CI)'; vartypes{end+1} = 'cell';
columnnames{end+1} = 'alpha (Mean)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'alpha (95% CI)'; vartypes{end+1} = 'cell';
columnnames{end+1} = 'beta (Mean)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'beta (95% CI)'; vartypes{end+1} = 'cell';
if logical(outputs.(fldnames{1}).Full.init.w_dir)    
    columnnames{end+1} = 's (Mean)'; vartypes{end+1} = 'double';
    columnnames{end+1} = 's (95% CI)'; vartypes{end+1} = 'cell';
end
if ismember('endm',fldnames) && (length(outputs.endm.Full.init.datalist_endm)>1)...
        && ~logical(outputs.endm.Full.init.datapca_endm)
    for i = 1:length(outputs.endm.Full.init.datalist_endm)
        columnnames{end+1} = sprintf('b%d (Mean)',i); vartypes{end+1} = 'double';
        columnnames{end+1} = sprintf('b%d (95% CI)',i); vartypes{end+1} = 'cell';
    end
    for i = 1:length(outputs.endm.Full.init.datalist_endm)
        columnnames{end+1} = sprintf('p%d (Mean)',i); vartypes{end+1} = 'double';
        columnnames{end+1} = sprintf('p%d (95% CI)',i); vartypes{end+1} = 'cell';
    end
else
    columnnames{end+1} = 'b (Mean)'; vartypes{end+1} = 'double';
    columnnames{end+1} = 'b (95% CI)'; vartypes{end+1} = 'cell';
    columnnames{end+1} = 'p (Mean)'; vartypes{end+1} = 'double';
    columnnames{end+1} = 'p (95% CI)'; vartypes{end+1} = 'cell';
end

ts = outputs.(fldnames{1}).Full.time_stamps;
for i = 1:length(ts)
    columnnames{end+1} = sprintf('R, t = %d',ts(i)); vartypes{end+1} = 'double';
end
columnnames{end+1} = 'Linear Model: Log-Likelihood'; vartypes{end+1} = 'double';
columnnames{end+1} = 'Linear Model: AIC'; vartypes{end+1} = 'double';
columnnames{end+1} = 'Linear Model: BIC'; vartypes{end+1} = 'double';
columnnames{end+1} = 'Linear Model: Intercept'; vartypes{end+1} = 'double';
columnnames{end+1} = 'Linear Model: Ordinary R^2'; vartypes{end+1} = 'double';
columnnames{end+1} = 'Linear Model: Adjusted R^2'; vartypes{end+1} = 'double';

summarytable = table('Size',[length(rownames),length(columnnames)],'VariableTypes',vartypes);
summarytable.Properties.RowNames = rownames;
summarytable.Properties.VariableNames = columnnames;
for k = 1:length(rownames)
    index = 1;
    if ismember('endm',fldnames)
        typenames = outputs.endm.Full.init.datalist_endm;
        if isnumeric(typenames)
            typenames = IndexName(typenames,outputs.endm.Full.init.datatype_endm);
        end

        for i = 1:length(typenames)
            if strcmp('ndm',rownames{k})
                summarytable{k,index} = "None"; index = index + 1;
            else
                summarytable{k,index} = string(typenames{i}); index = index + 1;
            end
        end
    else
        summarytable{k,index} = "None"; index = index + 1;
    end
    
    if ismember('endm',fldnames)
        if strcmp('ndm',rownames{k})
            summarytable{k,index} = "No"; index = index + 1;
        else
            if ~logical(outputs.endm.Full.init.datapca_endm)
                summarytable{k,index} = "No"; index = index + 1;
            else
                summarytable{k,index} = "Yes"; index = index + 1;
            end
        end        
    else
        summarytable{k,index} = "No"; index = index + 1;
    end    
        
    summarytable{k,index} = string(outputs.(fldnames{k}).Full.init.costfun); index = index + 1;
    
    subfldnames = fieldnames(outputs.(fldnames{k}));
    if length(subfldnames) == 1
        params = outputs.(fldnames{k}).Full.param_fit;
        inclinds = 1:length(params);
        if ~logical(outputs.(fldnames{k}).Full.init.w_dir)
            inclinds(4) = NaN;
        end
        if strcmp('endm',fldnames{k}) && (length(outputs.endm.Full.init.datalist_endm)>1)  && ...
                ~logical(outputs.endm.Full.init.datapca_endm)
            inclinds(5:(4+length(outputs.endm.Full.init.datalist_endm))) = NaN;
        else
            inclinds(5) = NaN;
        end
        params = params(~isnan(inclinds));
        if strcmp('ndm',fldnames{k}) && ismember('endm',fldnames) && ...
                (length(outputs.endm.Full.init.datalist_endm)>1) && ...
                ~logical(outputs.endm.Full.init.datapca_endm)
            params = [params, zeros(1,2*(length(outputs.endm.Full.init.datalist_endm)-1))];
        end
        for i = 1:length(params)
            summarytable{k,index} = params(i); index = index + 1;
            summarytable{k,index} = {[params(i) params(i)]}; index = index + 1;
        end
    else
        params = zeros((length(subfldnames)-1),length(outputs.(fldnames{k}).Full.param_fit));
        for i = 1:size(params,1)
            params(i,:) = outputs.(fldnames{k}).(subfldnames{i}).param_fit;
        end
        params_mean = mean(params); 
        params_ci95_lb = prctile(params,2.5,1); params_ci95_ub = prctile(params,97.5,1);
        params_ci95 = cat(1,params_ci95_lb,params_ci95_ub);
        inclinds = 1:length(params_mean);
        if ~logical(outputs.(fldnames{k}).Full.init.w_dir)
            inclinds(4) = NaN;
        end
        if strcmp('endm',fldnames{k}) && (length(outputs.endm.Full.init.datalist_endm)>1)  && ...
                ~logical(outputs.endm.Full.init.datapca_endm)
            inclinds(5:(4+length(outputs.endm.Full.init.datalist_endm))) = NaN;
        else
            inclinds(5) = NaN;
        end
        params_mean = params_mean(~isnan(inclinds)); 
        params_ci95 = params_ci95(:,~isnan(inclinds));
        if strcmp('ndm',fldnames{k}) && ismember('endm',fldnames) &&...
                (length(outputs.endm.Full.init.datalist_endm)>1) &&...
                ~logical(outputs.endm.Full.init.datapca_endm)
            params_mean = [params_mean, zeros(1,2*(length(outputs.endm.Full.init.datalist_endm)-1))];
            params_ci95 = [params_ci95, zeros(2,2*(length(outputs.endm.Full.init.datalist_endm))-1)];
        end
        for i = 1:length(params_mean)
            summarytable{k,index} = params_mean(i); index = index + 1;
            summarytable{k,index} = {params_ci95(:,i).'}; index = index + 1;
        end
    end
    ts = outputs.(fldnames{1}).Full.time_stamps;
    for i = 1:length(ts)
        summarytable{k,index} = outputs.(fldnames{k}).Full.results.Corrs(i); index = index + 1;       
    end
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_LogL; index = index + 1; 
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_AIC; index = index + 1;  
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_BIC; index = index + 1;
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_intercept; index = index + 1;    
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_Rsquared_ord; index = index + 1; 
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_Rsquared_adj;
end

if logical(writeout)
    writetable(summarytable,[filepath filesep filename '.csv'],'WriteRowNames',true)
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