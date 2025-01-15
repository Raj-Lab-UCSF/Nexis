function data_new = CCFToData(data_old, studyname_, rois, matdir_)

load([matdir_ filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct');
load([matdir_ filesep 'DefaultAtlas.mat'],'DefaultAtlas')
voxels_2hem = DefaultAtlas.volumes;
taustudies = fieldnames(mousedata_struct);
if ~ismember(studyname_,taustudies) && ~isempty(rois)
    data_new = NaN(length(rois),size(data_old,2)); % n ROI in data
    for i = 1:length(rois) 
        roi_inds_i = rois{i};
        vols_i = voxels_2hem(roi_inds_i);
        data_old_i = data_old(roi_inds_i,:);
        if all(isnan(data_old_i(:)))
            data_new(i,:) = NaN;
        else
            if any(isnan(data_old_i(:))) % possible missing data for CCF-parcellated data (Zhuang, e.g.)
                nanbool = isnan(data_old_i(:,1)); % should be consistent across rows
                data_old_i(nanbool,:) = [];
                vols_i(nanbool) = [];
            end
            data_new(i,:) = (vols_i.' * data_old_i) / sum(vols_i);
        end
    end
elseif ismember(studyname_,taustudies)
    rois = mousedata_struct.(studyname_).regions(:,2);
    data_new = NaN(length(rois),size(data_old,2)); % n ROI in data
    for i = 1:length(rois) 
        roi_inds_i = rois{i};
        vols_i = voxels_2hem(roi_inds_i);
        data_old_i = data_old(roi_inds_i,:);
        if all(isnan(data_old_i(:)))
            data_new(i,:) = NaN;
        else
            if any(isnan(data_old_i(:))) % possible missing data for CCF-parcellated data (Zhuang, e.g.)
                nanbool = isnan(data_old_i(:,1)); % should be consistent across rows
                data_old_i(nanbool,:) = [];
                vols_i(nanbool) = [];
            end
            data_new(i,:) = (vols_i.' * data_old_i) / sum(vols_i);
        end
    end
else
    error('Missing ROI information!');
end
end