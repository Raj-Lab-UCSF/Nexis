function data_new = DataToCCF(data_old, studyname_, rois, matdir_)

load([matdir_ filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct');
taustudies = fieldnames(mousedata_struct);
if isempty(data_old) && ismember(studyname_,taustudies)
    data_old = mousedata_struct.(studyname_).data;
end
if ismember(studyname_,taustudies) && isempty(rois)
    rois = mousedata_struct.(studyname_).regions(:,2);
end
data_new = NaN(426,size(data_old,2)); % n ROI in CCF
for i = 1:length(rois) 
    roi_inds_i = rois{i};
    for j = 1:length(roi_inds_i)
        data_new(roi_inds_i(j),:) = data_old(i,:);
    end
end
end