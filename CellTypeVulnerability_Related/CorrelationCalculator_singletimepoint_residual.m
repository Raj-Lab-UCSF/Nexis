function [corrmat_all,logpmat_all,Y] = CorrelationCalculator_singletimepoint_residual...
    (datsetnames_,tpt_,outputs,outstruct_nx,Xct_,corrtype_,seedcorr_,matdir_)
% Calculate and store correlations

corrmat_all = zeros(length(datsetnames_),size(Xct_,2));
logpmat_all = corrmat_all;
for j = 1:length(datsetnames_)
    datsetname = datsetnames_{j};
    Xct_j = CCFToData(Xct_,datsetname,matdir_);
    Y_all = outstruct_nx.(datsetname).global.nexis_global.Full.data; % normalized post-NexIS
    Y_pred = outstruct_nx.(datsetname).global.nexis_global.Full.predicted; % normalized post-NexIS
    if logical(seedcorr_) && (sum(outputs.(datsetname).seed)>0) % remove seed
        Y_all(logical(outputs.(datsetname).seed),:) = [];
        Y_pred(logical(outputs.(datsetname).seed),:) = [];
        Xct_j(logical(outputs.(datsetname).seed),:) = [];
    end
    if strcmp(tpt_,'end')
        t_ = length(outputs.(datsetname).time_stamps);
    else
        t_ = tpt_;
    end
    if strcmp(datsetname,'Hurtado')
        Y_d = Y_all(:,t_-1); % used baseline for global prediction
        Y_p = Y_pred(:,t_-1); % used baseline for global prediction
    else
        Y_d = Y_all(:,t_);
        Y_p = Y_pred(:,t_);
    end
    % mdl = fitlm(Y_p,Y_d,'Intercept',false);
    % Y = mdl.Residuals.Raw;
    Y = Y_d - Y_p;
    testmat = [Y Xct_j];
    if strcmp(corrtype_,'Pearson')
        [corrmat,corrp] = corrcoef(testmat);
        corrmat_all(j,:) = corrmat(1,2:end);
        logpmat_all(j,:) = -log10(corrp(1,2:end));
    elseif strcmp(corrtype_,'Partial')
        [corrmat,corrp] = partialcorr(testmat);
        corrmat_all(j,:) = corrmat(1,2:end);
        logpmat_all(j,:) = -log10(corrp(1,2:end));
    end
end
end