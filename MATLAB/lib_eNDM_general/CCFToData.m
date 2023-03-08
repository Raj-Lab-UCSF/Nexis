function data_new = CCFToData(data_old, studyname, matdir_)

load([matdir_ filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct');
load([matdir_ filesep 'regionvoxels.mat'],'voxels')
voxels_2hem = [voxels; voxels].';
rois = mousedata_struct.(studyname).regions(:,2);
data_new = NaN(length(rois),size(data_old,2)); % n ROI in CCF
for i = 1:length(rois) 
    roi_inds_i = rois{i};
    vols_i = voxels_2hem(roi_inds_i);
    data_new(i,:) = (vols_i * data_old(roi_inds_i,:)) / sum(vols_i);
end
end