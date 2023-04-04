function data_new = CCFToData(data_old, studyname, matdir_)

taustudies = {'IbaP301S','IbaHippInj','IbaStrInj','BoludaDSAD','BoludaCBD'...
    'DS1','DS4','DS6','DS6_110','DS7','DS9','DS9_110','DS10'};
if ismember(studyname,taustudies)
    load([matdir_ filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct');
    load([matdir_ filesep 'DefaultAtlas.mat'],'DefaultAtlas')
    voxels_2hem = DefaultAtlas.volumes;
    rois = mousedata_struct.(studyname).regions(:,2);
    data_new = NaN(length(rois),size(data_old,2)); % n ROI in CCF
    for i = 1:length(rois) 
        roi_inds_i = rois{i};
        vols_i = voxels_2hem(roi_inds_i);
        data_new(i,:) = (vols_i * data_old(roi_inds_i,:)) / sum(vols_i);
    end
else
    error('Function only valid for selected tau studies; invalid for study name %s',...
        studyname);
end