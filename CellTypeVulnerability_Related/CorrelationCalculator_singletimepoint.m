function [corrmat_all,logpmat_all] = CorrelationCalculator_singletimepoint...
    (datsetnames_,tpt_,outputs,Xct_,corrtype_,seedcorr_,matdir_)
% Calculate and store correlations

corrmat_all = zeros(length(datsetnames_),size(Xct_,2));
logpmat_all = corrmat_all;
for j = 1:length(datsetnames_)
    datsetname = datsetnames_{j};
    Xct_j = CCFToData(Xct_,datsetname,matdir_);
    Y_all = outputs.(datsetname).data;
    if logical(seedcorr_) && (sum(outputs.(datsetname).seed)>0)
%         Y_all = SeedCorrection(datsetname,matdir_);
        Y_all(logical(outputs.(datsetname).seed),:) = [];
        Xct_j(logical(outputs.(datsetname).seed),:) = [];
    end
    if strcmp(tpt_,'end')
        t_ = length(outputs.(datsetname).time_stamps);
    else
        t_ = tpt_;
    end
    Y = Y_all(:,t_);
    testmat = [Y Xct_j];
    if strcmp(corrtype_,'Pearson')
        [corrmat,corrp] = corrcoef(testmat);
        corrmat_all(j,:) = corrmat(1,2:end);
        logpmat_all(j,:) = -log10(corrp(1,2:end));
    elseif strcmp(corrtype_,'Spearman')
        [corrmat,corrp] = corr(testmat,'type','Spearman');
        corrmat_all(j,:) = corrmat(1,2:end);
        logpmat_all(j,:) = -log10(corrp(1,2:end));
    elseif strcmp(corrtype_,'Partial')
        [corrmat,corrp] = partialcorr(testmat);
        corrmat_all(j,:) = corrmat(1,2:end);
        logpmat_all(j,:) = -log10(corrp(1,2:end));
    end
end
end