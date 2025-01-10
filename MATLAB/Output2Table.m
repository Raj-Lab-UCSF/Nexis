function summarytable = Output2Table(outputs,writeout,filename,filepath)
% Function that unpacks output structs from NexIS_global.m or NexIS_SV.m
% into a table with relevant summary statistics and optimized parameter 
% values. NOTE: Changed from mean to median (1/2/25)

fldnames = fieldnames(outputs);
if nargin < 4
    filepath = cd;
    if nargin < 3
        datestr = datetime('today');
        if ~ismember('nexis_sv',fldnames)
            filename = ['summary_Nexis_mouse_' outputs.nexis_global.Full.init.study...
                '_global_' datestr];             
        else
            typename = outputs.nexis_sv.Full.init.datalist_nexis_sv(1);
            if isnumeric(typename)
                typename = IndexName(typename,outputs.nexis_sv.Full.init.datatype_nexis_sv);
            end
            filename = ['summary_Nexis_mouse_' outputs.nexis_sv.Full.init.study...
                '_' typename{1} '_' datestr];  
        end
        if nargin < 2
            writeout = 0;
        end
    end
end
            
rownames = fldnames.';
columnnames = cell(1,1); vartypes = columnnames;
% if ismember('nexis_sv',fldnames)
%     if strcmp(outputs.nexis_sv.Full.init.datatype_nexis_sv,'gene')
%         colstr1 = 'Gene';
%     else
%         colstr1 = 'Cell Type';
%     end
%     typenames = outputs.nexis_sv.Full.init.datalist_nexis_sv;
%     for i = 1:length(typenames)
%         columnnames{i} = sprintf([colstr1 ' %d'],i);
%         vartypes{i} = 'string';
%     end
% else
   % columnnames{1} = 'Gene 1'; 
   % vartypes{1} = 'string';
% end
columnnames{1} = 'SV Factor'; 
vartypes{1} = 'string';
columnnames{end+1} = 'Uses PCA'; vartypes{end+1} = 'string';    
columnnames{end+1} = 'Cost Function'; vartypes{end+1} = 'string';
columnnames{end+1} = 'gamma (Median)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'gamma (95% CI)'; vartypes{end+1} = 'cell';
columnnames{end+1} = 'alpha (Median)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'alpha (95% CI)'; vartypes{end+1} = 'cell';
columnnames{end+1} = 'beta (Median)'; vartypes{end+1} = 'double';
columnnames{end+1} = 'beta (95% CI)'; vartypes{end+1} = 'cell';
columnnames{end+1} = 's (Median)'; vartypes{end+1} = 'double';
columnnames{end+1} = 's (95% CI)'; vartypes{end+1} = 'cell';

if ismember('nexis_sv',fldnames) && ~strcmp(outputs.nexis_sv.Full.init.datatype_nexis_sv,'User_specified')
    if ismember('nexis_sv',fldnames) && (length(outputs.nexis_sv.Full.init.datalist_nexis_sv)>1)...
            && ~logical(outputs.nexis_sv.Full.init.datapca_nexis_sv)
        for i = 1:length(outputs.nexis_sv.Full.init.datalist_nexis_sv)
            columnnames{end+1} = sprintf('b%d (Median)',i); vartypes{end+1} = 'double';
            columnnames{end+1} = sprintf('b%d (95%% CI)',i); vartypes{end+1} = 'cell';
        end
        for i = 1:length(outputs.nexis_sv.Full.init.datalist_nexis_sv)
            columnnames{end+1} = sprintf('p%d (Median)',i); vartypes{end+1} = 'double';
            columnnames{end+1} = sprintf('p%d (95%% CI)',i); vartypes{end+1} = 'cell';
        end
    else
        columnnames{end+1} = 'b (Median)'; vartypes{end+1} = 'double';
        columnnames{end+1} = 'b (95% CI)'; vartypes{end+1} = 'cell';
        columnnames{end+1} = 'p (Median)'; vartypes{end+1} = 'double';
        columnnames{end+1} = 'p (95% CI)'; vartypes{end+1} = 'cell';
    end
else
    if ismember('nexis_sv',fldnames) && (size(outputs.nexis_sv.Full.init.datalist_nexis_sv,2)>1)...
            && ~logical(outputs.nexis_sv.Full.init.datapca_nexis_sv)
        for i = 1:size(outputs.nexis_sv.Full.init.datalist_nexis_sv,2)
            columnnames{end+1} = sprintf('b%d (Median)',i); vartypes{end+1} = 'double';
            columnnames{end+1} = sprintf('b%d (95%% CI)',i); vartypes{end+1} = 'cell';
        end
        for i = 1:size(outputs.nexis_sv.Full.init.datalist_nexis_sv,2)
            columnnames{end+1} = sprintf('p%d (Median)',i); vartypes{end+1} = 'double';
            columnnames{end+1} = sprintf('p%d (95%% CI)',i); vartypes{end+1} = 'cell';
        end
    else
        columnnames{end+1} = 'b (Median)'; vartypes{end+1} = 'double';
        columnnames{end+1} = 'b (95% CI)'; vartypes{end+1} = 'cell';
        columnnames{end+1} = 'p (Median)'; vartypes{end+1} = 'double';
        columnnames{end+1} = 'p (95% CI)'; vartypes{end+1} = 'cell';
    end

end

ts = outputs.(fldnames{1}).Full.time_stamps;
if isnan(outputs.(fldnames{1}).Full.init.seed)
    ts = ts(2:end);
end
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
    if ismember('nexis_sv',fldnames)
        typenames = outputs.nexis_sv.Full.init.datalist_nexis_sv;
        if isnumeric(typenames) && ~strcmp(outputs.nexis_sv.Full.init.datatype_nexis_sv,'User_specified')
            typenames = IndexName(typenames,outputs.nexis_sv.Full.init.datatype_nexis_sv);
        elseif strcmp(outputs.nexis_sv.Full.init.datatype_nexis_sv,'User_specified')
            typenames = {'User specified'};
        end
        
        if length(typenames) > 1
            typestr = [];
            for m = 1:length(typenames)
                typestr = [typestr typenames{m} ', '];
            end
            typestr = typestr(1:(end-2));
        else
            typestr = typenames{1};
        end
        if strcmp('nexis_global',rownames{k})
            summarytable{k,index} = "None"; index = index + 1;
        else
            summarytable{k,index} = string(typestr); index = index + 1;
        end
    else
        summarytable{k,index} = "None"; index = index + 1;
    end
    
    if ismember('nexis_sv',fldnames)
        if strcmp('nexis_global',rownames{k})
            summarytable{k,index} = "No"; index = index + 1;
        else
            if ~logical(outputs.nexis_sv.Full.init.datapca_nexis_sv)
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
        % if ~logical(outputs.(fldnames{k}).Full.init.w_dir)
        %     inclinds(4) = NaN;
        % end
        % if strcmp('nexis_sv',fldnames{k}) && (length(outputs.nexis_sv.Full.init.datalist_nexis_sv)>1)  && ...
        %         ~logical(outputs.nexis_sv.Full.init.datapca_nexis_sv)
        %     inclinds(5:(4+length(outputs.nexis_sv.Full.init.datalist_nexis_sv))) = NaN;
        % else
        %     inclinds(5) = NaN;
        % end
        params = params(~isnan(inclinds));
        if ismember('nexis_sv',fldnames) && ~isequal(typenames,{'User specified'})
            if strcmp('nexis_global',fldnames{k}) && ismember('nexis_sv',fldnames) && ...
                    (length(outputs.nexis_sv.Full.init.datalist_nexis_sv)>1) && ...
                    ~logical(outputs.nexis_sv.Full.init.datapca_nexis_sv)
                params = [params, zeros(1,2*(length(outputs.nexis_sv.Full.init.datalist_nexis_sv)-1))];
            end
        else
            if strcmp('nexis_global',fldnames{k}) && ismember('nexis_sv',fldnames) && ...
                    (size(outputs.nexis_sv.Full.init.datalist_nexis_sv,2)>1) && ...
                    ~logical(outputs.nexis_sv.Full.init.datapca_nexis_sv)
                params = [params, zeros(1,2*(size(outputs.nexis_sv.Full.init.datalist_nexis_sv,2)-1))];
            end
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
        params_median = median(params);
        params_ci95_lb = prctile(params,2.5,1); params_ci95_ub = prctile(params,97.5,1);
        params_ci95 = cat(1,params_ci95_lb,params_ci95_ub);
        inclinds = 1:length(params_median);
        % if ~logical(outputs.(fldnames{k}).Full.init.w_dir)
        %     inclinds(4) = NaN;
        % end
        % if strcmp('nexis_sv',fldnames{k}) && (length(outputs.nexis_sv.Full.init.datalist_nexis_sv)>1)  && ...
        %         ~logical(outputs.nexis_sv.Full.init.datapca_nexis_sv)
        %     inclinds(5:(4+length(outputs.nexis_sv.Full.init.datalist_nexis_sv))) = NaN;
        % else
        %     inclinds(5) = NaN;
        % end
        params_median = params_median(~isnan(inclinds)); 
        params_ci95 = params_ci95(:,~isnan(inclinds));
        if ismember('nexis_sv',fldnames) && ~isequal(typenames,{'User specified'})
            if strcmp('nexis_global',fldnames{k}) && ismember('nexis_sv',fldnames) &&...
                    (length(outputs.nexis_sv.Full.init.datalist_nexis_sv)>1) &&...
                    ~logical(outputs.nexis_sv.Full.init.datapca_nexis_sv)
                params_median = [params_median, zeros(1,2*(length(outputs.nexis_sv.Full.init.datalist_nexis_sv)-1))];
                params_ci95 = [params_ci95, zeros(2,2*(length(outputs.nexis_sv.Full.init.datalist_nexis_sv))-1)];
            end
        else
            if strcmp('nexis_global',fldnames{k}) && ismember('nexis_sv',fldnames) &&...
                    (size(outputs.nexis_sv.Full.init.datalist_nexis_sv,2)>1) &&...
                    ~logical(outputs.nexis_sv.Full.init.datapca_nexis_sv)
                params_median = [params_median, zeros(1,2*(size(outputs.nexis_sv.Full.init.datalist_nexis_sv,2)-1))];
                params_ci95 = [params_ci95, zeros(2,2*(size(outputs.nexis_sv.Full.init.datalist_nexis_sv,2))-1)];
            end
        end
        for i = 1:length(params_median)
            summarytable{k,index} = params_median(i); index = index + 1;
            summarytable{k,index} = {params_ci95(:,i).'}; index = index + 1;
        end
    end

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

    function names = IndexName(indices,dattypenexis_sv)
        if strcmp(dattypenexis_sv,'gene')
            datstruct = load([cd filesep 'raw_data_mouse' filesep 'GeneExpressionMaps.mat'],'GeneExpressionMaps');
            namescell = datstruct.GeneExpressionMaps.All.gene_names;
        else
            datstruct = load([cd filesep 'raw_data_mouse' filesep 'CellTypeMaps.mat'],'CellTypeMaps');
            namescell = datstruct.CellTypeMaps.(dattypenexis_sv).classkey;
        end
        names = namescell(indices);
    end
end