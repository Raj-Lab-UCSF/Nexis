addpath('~/Documents/MATLAB/Brainframe-Dev/Brainframe/');
matdir = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
matdir1 = '/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe';
matdir2 = '/Users/justintorok/Documents/MATLAB/MISS/MISS-MatFiles';
load([matdir2 filesep 'Zeisel_Inputs.mat'],'GENGDmod','structList','structIndex','nonzerovox','classkey');
load([matdir2 filesep 'Zeisel_outstruct_nG1360.mat'],'outstruct');
load([matdir filesep 'classkey_zeisel.mat']);

hbnorind = find(ismember(classkey_zeisel,'HBNOR'));
curinput = outstruct(1).corrB(:,hbnorind);
newVoxMap = zeros(size(GENGDmod));
newVoxMap(nonzerovox) = curinput;
datinput = imresize3(newVoxMap,[133 81 115]);
datinput(datinput < 0) = 0;


input_struct_hbnor = brainframe_inputs_mouse(matdir1,'voxUreg',0,...
                                                'data',datinput,...
                                                'centered',[0,1],...
                                                'xfac',0.75,...
                                                'voxthresh',0.65,...
                                                'pointsize',0.1,...
                                                'nbin',10,...
                                                'cmap',twocolor([1 0.5 0.5],[1 0 0],10),...
                                                'savenclose',1,...
                                                'img_labels','HBNOR');

mbdop1ind = find(ismember(classkey_zeisel,'MBDOP1'));
curinput = outstruct(1).corrB(:,mbdop1ind);
newVoxMap = zeros(size(GENGDmod));
newVoxMap(nonzerovox) = curinput;
datinput = imresize3(newVoxMap,[133 81 115]);
datinput(datinput < 0) = 0;


input_struct_mbdop1 = brainframe_inputs_mouse(matdir1,'voxUreg',0,...
                                                'data',datinput,...
                                                'centered',[0,1],...
                                                'xfac',0.75,...
                                                'voxthresh',0.65,...
                                                'pointsize',0.1,...
                                                'nbin',10,...
                                                'cmap',twocolor([1 0.5 1],[1 0 1],10),...
                                                'savenclose',1,...
                                                'img_labels','MBDOP1');

brainframe(input_struct_hbnor);
brainframe(input_struct_mbdop1);
% outstruct_hbnor = struct;
% outstruct_hbnor.corrB = hbnor_zscore;
% outstruct_hbnor.nGen = 1360;
% MISS_Brainframe(outstruct_hbnor,1,'Zeisel',1,xfac,0,0.75,{[1 0 0; 1 1 1]},'hbnorcorr',matdir);
% hbnorind = find(ismember(classkey,'HBNOR'));
% MISS_Brainframe(outstruct,hbnorind,'Zeisel',1,xfac,0,0.8,{[0 0 1; 1 1 1]},'hbnorcorr',matdir);
% 
% input_struct.voxUreg = 1; 
% input_struct.nbins = 10;
% input_struct.data = zeros(426,1); input_struct.data([162,162+213]) = 1;
% input_struct.cmap = [1 0 0];
% input_struct.pointsize = 0.0001;
% input_struct.voxthresh = 0.75;
% input_struct.xfac = 5;
% input_struct.regsUbins = 0;
% input_struct.centered = [0 1];
% brainframe(input_struct);




% hbnor = genevct(:,ismember(classkey,'HBNOR'));
% hbnor_mrx3 = hbnor(geneinds(1:1360));
% voxvgene_mrx3 = voxvgene(:,geneinds(1:1360)).';
% hbnor_corr_map = corr(voxvgene_mrx3,hbnor_mrx3);
% hbnor_zscore = (hbnor_corr_map - mean(hbnor_corr_map))/std(hbnor_corr_map);
% hbnor_zscore(hbnor_zscore < 2) = 0;
% hbnor_zscore_map = zeros(size(GENGDmod));
% hbnor_zscore_map(nonzerovox) = hbnor_zscore;
% hbnor_zscore_map = imresize3(hbnor_zscore_map,[133 81 115]);
% hbnor_zscore_map(hbnor_zscore_map < 0) = 0;
% slc6a2_expr = voxvgene(:,ismember(gene_names,'Slc6a2'));
% slc6a2_expr_map = zeros(size(GENGDmod));
% slc6a2_expr_map(nonzerovox) = slc6a2_expr;
% slc6a2_expr_map = imresize3(slc6a2_expr_map,[133 81 115]);
% slc6a2_expr_map(slc6a2_expr_map < 0) = 0;
