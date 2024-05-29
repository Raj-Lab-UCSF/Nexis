function BrainframeSeedConnectivityPlot(outstruct,datset,tptsplotind_,C_,Cdir_,matdir_,savenclose_,figdir_)

% Define major region groups and colors
brainframedir = '/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe';
addpath(brainframedir)

reggroups_ = zeros(213,1); %Chunk of code to define region_groups
amy = 1:11; cer = 12:23; sub = 24:26; hip = 27:37; hyp = 38:57;
ncx = 58:95; med = 96:120; mid = 121:141; olf = 142:149; pal = 150:157;
pon = 158:170; str = 171:178; tha = 179:213;
reggroups_(amy) = 1; reggroups_(cer) = 2; reggroups_(sub) = 3; 
reggroups_(hip) = 4; reggroups_(hyp) = 5; reggroups_(ncx) = 6;
reggroups_(med) = 7; reggroups_(mid) = 8; reggroups_(olf) = 9;
reggroups_(pal) = 10; reggroups_(pon) = 11; reggroups_(str) = 12;
reggroups_(tha) = 13;
reggroups_ = [reggroups_;reggroups_];
cmap_ = hsv(length(unique(reggroups_)));

% Transform data & seed to CCF space, obtain queried time point
datinput_data = DataToCCF([],datset,matdir_);
datinput_data(isnan(datinput_data)) = 0;
seedinput_data = outstruct.(datset).seed;
tpt = outstruct.(datset).time_stamps;
if isnan(seedinput_data)
    seedinput_data = logical(datinput_data(:,1));
    toffset = 1;
else
    seedinput_data = DataToCCF(seedinput_data,datset,matdir_);
    seedinput_data(isnan(seedinput_data)) = 0;
    seedinput_data = logical(seedinput_data);
    toffset = 0;
end
reggroups_data_ = reggroups_;
reggroups_data_(seedinput_data) = 14;
cmap_data_ = [cmap_; [1 0.5 0]];

tptsplotind_ = tptsplotind_ + toffset;
datinput_data = datinput_data(:,tptsplotind_);
tpt = tpt(tptsplotind_);

% Thresholding connectivity and data
conthresh = 90;
C_thresh_seed = zeros(size(C_));
if strcmp(Cdir_,'In')
    C_thresh_seed(:,seedinput_data) = C_(:,seedinput_data);
elseif strcmp(Cdir_,'Out')
    C_thresh_seed(seedinput_data,:) = C_(seedinput_data,:);
elseif strcmp(Cdir_,'Both')
    C_thresh_seed(:,seedinput_data) = C_(:,seedinput_data);
    C_thresh_seed(seedinput_data,:) = C_(seedinput_data,:);
end
conthresh_val = prctile(nonzeros(C_thresh_seed(:)),conthresh);
thresh_inds_C = (C_thresh_seed >= conthresh_val);
C_thresh_seed(~thresh_inds_C) = 0;

datathresh = 50;
datathresh_val = prctile(nonzeros(datinput_data),datathresh);
thresh_inds_data = (datinput_data >= datathresh_val);
datinput_data(~thresh_inds_data) = 0;
datinput_data(seedinput_data) = 1.5*max(datinput_data);

% Genereate glass brain
imglabel = sprintf('%s_%sSeedConn_t%d_datathresh%d_Cthresh%d',datset,...
    Cdir_,tpt,datathresh,conthresh);
imglabel = strrep(imglabel,'.','_');
imgview = [-90,-18];
input_struct_seedconn = brainframe_inputs_mouse(brainframedir,'conmat',C_thresh_seed,...
                                             'region_groups',reggroups_data_,...
                                             'con_regiongroups',reggroups_,...
                                             'cmap',cmap_data_,...
                                             'con_cmap',cmap_,...
                                             'xfac',4.5,...
                                             'sphere',1,...
                                             'sphere_npts',20,...
                                             'pointsize',5,...
                                             'voxUreg',1,...
                                             'iscon',1,...
                                             'conarrow_WL',[1.5 1],...
                                             'data',datinput_data,...
                                             'norm_method','max',...
                                             'bgcolor','w',...
                                             'con_rescale',20,...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figdir_,...
                                             'savenclose',savenclose_,...
                                             'con_arch',0.3);
brainframe(input_struct_seedconn);

end