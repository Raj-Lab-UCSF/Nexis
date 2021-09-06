function summarytable = Output2Table(outputs,writeout,filename,filepath)

fldnames = fieldnames(outputs);
if nargin < 4
    filepath = cd;
    if nargin < 3
        datestr = date;
        if ismember('ndm',fldnames)
            filename = ['summary_eNDM_mouse_' outputs.ndm.Full.init.study '_' datestr];             
        elseif ismember('endm',fldnames)
            filename = ['summary_eNDM_mouse_' outputs.endm.Full.init.study '_' datestr];  
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
columnnames{end+1} = 'r (mean)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'r (std)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'alpha (mean)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'alpha (std)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'beta (mean)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'beta (std)'; vartypes{end+1} = 'double';
if logical(outputs.(fldnames{1}).Full.init.w_dir)    
    columnnames{end+1} = 's (mean)'; vartypes{end+1} = 'double';
    columnnames{end+1} = 's (std)'; vartypes{end+1} = 'double';
end
if ismember('endm',fldnames) && (length(outputs.endm.Full.init.datalist_endm)>1)...
        && ~logical(outputs.endm.Full.init.datapca_endm)
    for i = 1:length(outputs.endm.Full.init.datalist_endm)
        columnnames{end+1} = sprintf('b%d (mean)',i); vartypes{end+1} = 'double';
        columnnames{end+1} = sprintf('b%d (std)',i); vartypes{end+1} = 'double';
    end
    for i = 1:length(outputs.endm.Full.init.datalist_endm)
        columnnames{end+1} = sprintf('p%d (mean)',i); vartypes{end+1} = 'double';
        columnnames{end+1} = sprintf('p%d (std)',i); vartypes{end+1} = 'double';
    end
else
    columnnames{end+1} = 'b (mean)'; vartypes{end+1} = 'double';
    columnnames{end+1} = 'b (std)'; vartypes{end+1} = 'double';
    columnnames{end+1} = 'p (mean)'; vartypes{end+1} = 'double';
    columnnames{end+1} = 'p (std)'; vartypes{end+1} = 'double';
end

ts = outputs.(fldnames{1}).Full.time_stamps;
for i = 1:length(ts)
    columnnames{end+1} = sprintf('R, t = %d',ts(i)); vartypes{end+1} = 'double';
end
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
            summarytable{k,index} = 0; index = index + 1;
        end
    else
        params = zeros((length(subfldnames)-1),length(outputs.(fldnames{k}).Full.param_fit));
        for i = 1:size(params,1)
            params(i,:) = outputs.(fldnames{k}).(subfldnames{i}).param_fit;
        end
        params_mean = mean(params); params_std = std(params);
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
        params_std = params_std(~isnan(inclinds));
        if strcmp('ndm',fldnames{k}) && ismember('endm',fldnames) &&...
                (length(outputs.endm.Full.init.datalist_endm)>1) &&...
                ~logical(outputs.endm.Full.init.datapca_endm)
            params_mean = [params_mean, zeros(1,2*(length(outputs.endm.Full.init.datalist_endm)-1))];
            params_std = [params_std, zeros(1,2*(length(outputs.endm.Full.init.datalist_endm))-1)];
        end
        for i = 1:length(params_mean)
            summarytable{k,index} = params_mean(i); index = index + 1;
            summarytable{k,index} = params_std(i); index = index + 1;
        end
    end
    ts = outputs.(fldnames{1}).Full.time_stamps;
    for i = 1:length(ts)
        summarytable{k,index} = outputs.(fldnames{k}).Full.results.Corrs(i); index = index + 1;       
    end
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_AIC; index = index + 1;  
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_BIC; index = index + 1;
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_intercept; index = index + 1;    
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_Rsquared_ord; index = index + 1; 
    summarytable{k,index} = outputs.(fldnames{k}).Full.results.lm_Rsquared_adj;
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

if logical(writeout)
    writetable(summarytable,[filepath filesep filename '.csv'],'WriteRowNames',true)
end
end