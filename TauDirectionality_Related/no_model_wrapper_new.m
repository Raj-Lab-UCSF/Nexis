%% 0. Loading data
clear; clc; close all; rng(0);
matdir = '/Users/justintorok/Documents/MATLAB/TauDirectionality_ChrisPaper_Project/TauDirectionality/raw_data_mouse';
figdir = '/Users/justintorok/Documents/MATLAB/TauDirectionality_ChrisPaper_Project/Figures';
load([matdir filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct');
load([matdir filesep 'Connectomes.mat'],'Connectomes')
C_ret = Connectomes.default; C_ant = C_ret.'; C_nd = (C_ret + C_ant)/2;
C_ret = (C_ret - diag(diag(C_ret)));
C_ant = (C_ant - diag(diag(C_ant)));
C_nd = (C_nd - diag(diag(C_nd)));

% taudata_all = struct;
% tpts_all = struct;
% seed_all = struct;
% labnms_all = labnms;
% for i = 1:length(labnms)
%     taudata_all.(labnms{i}) = data426.(labnms{i});
%     tpts_all.(labnms{i}) = tpts.(labnms{i});
%     seed_all.(labnms{i}) = seed426.(labnms{i});
% end
% clearvars -except C_ret C_ant C_nd taudata_all tpts_all labnms_all seed_all netlabs matdir figdir   
% 
% load([matdir filesep 'KaufmanDiamond_datasets_dat&seed.mat'],'tpts','data426','seed426','labnms');
% labnms_all = [labnms_all.' labnms];
% for i = 1:length(labnms)
%     taudata_all.(labnms{i}) = data426.(labnms{i});
%     tpts_all.(labnms{i}) = tpts.(labnms{i});
%     seed_all.(labnms{i}) = seed426.(labnms{i});
% end
% clearvars -except C_ret C_ant C_nd taudata_all tpts_all labnms_all seed_all netlabs matdir figdir

%% 1. Graph analysis
L_ret = genLplcns(C_ret);
[v_ret,d_ret] = eig(L_ret); d_ret = abs(diag(d_ret)); v_ret = abs(v_ret);
[dretsort,sortinds] = sort(d_ret); vretsort = v_ret(:,sortinds);

L_ant = genLplcns(C_ant);
[v_ant,d_ant] = eig(L_ant); d_ant = abs(diag(d_ant)); v_ant = abs(v_ant);
[dantsort,sortinds] = sort(d_ant); vantsort = v_ant(:,sortinds);

%% 2.1 Figure 1, scatterplots
% Hurtado & IbaHippInj
hurtado_end = mousedata_struct.Hurtado.data(:,end);
ibahippinj_end = mousedata_struct.IbaStrInj.data(:,end);

outdeg = sum(C_ret,2); 
outdeg_hurtado = CCFToData(outdeg,'Hurtado',matdir);
outdeg_ibahippinj = CCFToData(outdeg,'IbaStrInj',matdir);

indeg = sum(C_ret,1).'; 
indeg_hurtado = CCFToData(indeg,'Hurtado',matdir);
indeg_ibahippinj = CCFToData(indeg,'IbaStrInj',matdir);

u1_ret = v_ret(:,1); 
u1_ret_hurtado = CCFToData(u1_ret,'Hurtado',matdir);
u1_ret_ibahippinj = CCFToData(u1_ret,'IbaStrInj',matdir);

u1_ant = v_ant(:,1); 
u1_ant_hurtado = CCFToData(u1_ant,'Hurtado',matdir);
u1_ant_ibahippinj = CCFToData(u1_ant,'IbaStrInj',matdir);

ylim1 = [0 max(ibahippinj_end)];
ylim2 = [0 max(hurtado_end)];
xlim1 = [0 max(outdeg_ibahippinj)];
xlim2 = [0 max(indeg_ibahippinj)];
xlim3 = [0 max(u1_ret_ibahippinj)];
xlim4 = [0 max(u1_ant_ibahippinj)];

figure('Units','inches','Position',[0 0 11 6]); 
tiledlayout(2,4,'TileSpacing','compact','Padding','tight');
nexttile;
scatter(outdeg_ibahippinj,ibahippinj_end,'bo','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim1); xlim(xlim1);
xticks([xlim1(1) (xlim1(1)+xlim1(2))/2 xlim1(2)]);
yticks([ylim1(1) (ylim1(1)+ylim1(2))/2 ylim1(2)]);
xticklabels({''});
ylabel('IbaHippInj Pathology')
text(0.6*xlim1(2),0.1*ylim1(2),sprintf('R = %.2f',corr(outdeg_ibahippinj,ibahippinj_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(indeg_ibahippinj,ibahippinj_end,'ro','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim1); xlim(xlim2); xticklabels({''}); yticklabels({''}); 
xticks([xlim2(1) (xlim2(1)+xlim2(2))/2 xlim2(2)]);
yticks([ylim1(1) (ylim1(1)+ylim1(2))/2 ylim1(2)]);
text(0.6*xlim2(2),0.1*ylim1(2),sprintf('R = %.2f',corr(indeg_ibahippinj,ibahippinj_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(u1_ret_ibahippinj,ibahippinj_end,'bd','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim1); xlim(xlim3); xticklabels({''}); yticklabels({''});
xticks([xlim3(1) (xlim3(1)+xlim3(2))/2 xlim3(2)]);
yticks([ylim1(1) (ylim1(1)+ylim1(2))/2 ylim1(2)]);
text(0.6*xlim3(2),0.1*ylim1(2),sprintf('R = %.2f',corr(u1_ret_ibahippinj,ibahippinj_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(u1_ant_ibahippinj,ibahippinj_end,'rd','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim1); xlim(xlim4); xticklabels({''}); yticklabels({''});
xticks([xlim4(1) (xlim4(1)+xlim4(2))/2 xlim4(2)]);
yticks([ylim1(1) (ylim1(1)+ylim1(2))/2 ylim1(2)]);
text(0.6*xlim4(2),0.1*ylim1(2),sprintf('R = %.2f',corr(u1_ant_ibahippinj,ibahippinj_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(outdeg_hurtado,hurtado_end,'bo','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim2); xlim(xlim1);
xticks([xlim1(1) (xlim1(1)+xlim1(2))/2 xlim1(2)]);
yticks([ylim2(1) (ylim2(1)+ylim2(2))/2 ylim2(2)]);
xtickformat('%.1d');
xlabel('Out Degree'); ylabel('Hurtado Pathology');
text(0.6*xlim1(2),0.1*ylim2(2),sprintf('R = %.2f',corr(outdeg_hurtado,hurtado_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(indeg_hurtado,hurtado_end,'ro','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim2); xlim(xlim2); yticklabels({''});
xticks([xlim2(1) (xlim2(1)+xlim2(2))/2 xlim2(2)]);
yticks([ylim2(1) (ylim2(1)+ylim2(2))/2 ylim2(2)]);
xtickformat('%.1d');
xlabel('In Degree'); 
text(0.6*xlim2(2),0.1*ylim2(2),sprintf('R = %.2f',corr(indeg_hurtado,hurtado_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(u1_ret_hurtado,hurtado_end,'bd','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim2); xlim(xlim3); yticklabels({''}); 
xticklabels({'0',num2str((xlim3(1)+xlim3(2))/2,'%.2f'),num2str(xlim3(2),'%.2f')})
xticks([xlim3(1) (xlim3(1)+xlim3(2))/2 xlim3(2)]);
yticks([ylim2(1) (ylim2(1)+ylim2(2))/2 ylim2(2)]);
xlabel('u1_r_e_t');
text(0.6*xlim3(2),0.1*ylim2(2),sprintf('R = %.2f',corr(u1_ret_hurtado,hurtado_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

nexttile;
scatter(u1_ant_hurtado,hurtado_end,'rd','filled'); 
l = lsline; l.LineWidth = 2; l.Color = 'k';
ylim(ylim2); xlim(xlim4); yticklabels({''});
xticklabels({'0',num2str((xlim4(1)+xlim4(2))/2,'%.2f'),num2str(xlim4(2),'%.2f')})
xticks([xlim4(1) (xlim4(1)+xlim4(2))/2 xlim4(2)]);
yticks([ylim2(1) (ylim2(1)+ylim2(2))/2 ylim2(2)]);
xlabel('u1_a_n_t'); 
text(0.6*xlim4(2),0.1*ylim2(2),sprintf('R = %.2f',corr(u1_ant_hurtado,hurtado_end)),...
    'FontSize',16);
set(gca,'FontSize',16);

% print([figdir filesep 'NoModelScatterplots'],'-dtiffn','-r300'); close;

%% 2.2 Figure 1, brainframes
matdir1 = '/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe';
figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
addpath('/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe')
load([matdir filesep 'default_mouse.mat'],'input_struct');

reggroups = zeros(213,1); %Chunk of code to define region_groups
amy = 1:11; cer = 12:23; sub = 24:26; hip = 27:37; hyp = 38:57; 
ncx = 58:95; med = 96:120; mid = 121:141; olf = 142:149; pal = 150:157;
pon = 158:170; str = 171:178; tha = 179:213;
reggroups(amy) = 1; reggroups(cer) = 2; reggroups(sub) = 3; 
reggroups(hip) = 4; reggroups(hyp) = 5; reggroups(ncx) = 6;
reggroups(med) = 7; reggroups(mid) = 8; reggroups(olf) = 9;
reggroups(pal) = 10; reggroups(pon) = 11; reggroups(str) = 12;
reggroups(tha) = 13;
reggroups = [reggroups;reggroups];
cmap = lines(length(unique(reggroups))); %Creating colormap

%% 2.2.1 IbaHippInj end pathology
datinput = ibahippinj_end;
datinput(isnan(datinput)) = 0;
imglabel = 'IbaHippInj_End';
imgview = [0 0 1];
input_struct_ibahipp = brainframe_inputs_mouse(matdir1,'region_groups',reggroups,...
                                             'cmap',cmap,...
                                             'xfac',10,...
                                             'pointsize',1,...
                                             'voxUreg',1,...
                                             'data',datinput,...
                                             'norm_method','mean',...
                                             'center',[1 2],...
                                             'bgcolor','w',...
                                             'con_rescale',20,...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figdir,...
                                             'savenclose',0);
brainframe(input_struct_ibahipp);

%% 2.2.2 Hurtado end pathology
% datinput = taudata_all.Hurtado(:,end);
datinput(isnan(datinput)) = 0;
imglabel = 'Hurtado_End';
imgview = [0 0 -1];
input_struct_hurtado = brainframe_inputs_mouse(matdir1,'region_groups',reggroups,...
                                             'cmap',cmap,...
                                             'xfac',6,...
                                             'pointsize',15,...
                                             'voxUreg',1,...
                                             'data',datinput,...
                                             'norm_method','mean',...
                                             'center',[1 1.5],...
                                             'sphere',0,...
                                             'sphere_npts',15,...
                                             'bgcolor','w',...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figdir,...
                                             'savenclose',0);
brainframe(input_struct_hurtado);
view([0 0 1])

%% 4. Figure 2, model
studyname = 'IbaHippInj';
rescale = 'log';
tvec = [0,mousedata_struct.(studyname).time_stamps];
bvec = 0:0.05:1;
t = 5;
tinds = find(ismember(t*bvec,tvec));

if strcmp(rescale,'log')
    pathology_sn = log(mousedata_struct.(studyname).data+1);
end
seed_sn = mousedata_struct.(studyname).seed;
pathwseed = [seed_sn, pathology_sn];

rate_sn = (pathwseed(:,2:end) - pathwseed(:,1:end-1))./...
    (tvec_data(2:end) - tvec_data(1:end-1));

seed_sn_ccf = DataToCCF(seed_sn,studyname,matdir);
model_ant_mat = zeros(length(seed_sn_ccf),length(tvec));
model_ret_mat = model_ant_mat;
model_nd_mat = model_ant_mat;

Rvals_ant = zeros(1,length(bvec)); Rvals_ant_rate = Rvals_ant;
Rvals_ret = zeros(1,length(bvec)); Rvals_ret_rate = Rvals_ret;
Rvals_nd = zeros(1,length(tvec_data)); Rvals_nd_rate = Rvals_nd;
    Rvals_ant(i) = corr(Xi_ant(:,end),pathology_sn(:,end),'rows','complete');
    Rvals_ant_rate(i) = corr(Xi_ant_rate(:,end),rate_sn(:,end),'rows','complete');
    model_ant_mat(:,:,i) = Xi_ant;

Xi_ant = spreadModel(seed_sn,tvec_model,L_ant,bvec);
if strcmp(rescale,'log')
    Xi_ant = log(Xi_ant + 1);
end
Xi_ant_rate = (Xi_ant(:,tinds(2:end)) - Xi_ant(:,tinds(1:end-1))) ./ ...
    (tvec_model(tinds(2:end)) - tvec_model(tinds(1:end-1)));

Xi_ret = spreadModel(seed_sn,tvec_model,L_ret,bvec);
if strcmp(rescale,'log')
    Xi_ret = log(Xi_ret + 1);
end
Xi_ret_rate = (Xi_ret(:,tinds(2:end)) - Xi_ret(:,tinds(1:end-1))) ./ ...
    (tvec_model(tinds(2:end)) - tvec_model(tinds(1:end-1)));

Xi_nd = spreadModel(seed_sn,tvec_model,(L_ret+L_ant)/2,bvec);
if strcmp(rescale,'log')
    Xi_nd = log(Xi_nd + 1);
end
Xi_nd_rate = (Xi_nd(:,tinds(2:end)) - Xi_nd(:,tinds(1:end-1))) ./ ...
    (tvec_model(tinds(2:end)) - tvec_model(tinds(1:end-1)));




% log(X_t + 1) for constructing R-beta*t curves, both pathology and model
figure; hold on;
plot(bvec*tvec_model(end),Rvals_ret_rate); plot(bvec*tvec_model(end),Rvals_ant_rate); plot(bvec*tvec_model(end),Rvals_nd_rate); 
legend({'ret','ant','nd'}); ylabel('R'); xlabel('\beta t'); set(gca,'FontSize',16); title('Rate R')

figure; hold on;
plot(bvec*tvec_model(end),Rvals_ret); plot(bvec*tvec_model(end),Rvals_ant); plot(bvec*tvec_model(end),Rvals_nd); 
legend({'ret','ant','nd'}); ylabel('R'); xlabel('\beta t'); set(gca,'FontSize',16); title('Pathology R')


%% 3. Conn. from seed
ibahippinj_seed = logical(seed_all.IbaHippInj);
ibahippinj_end = taudata_all.IbaHippInj(:,2);
ibahippinj_end(logical(ibahippinj_seed)) = NaN;
ibahippinj_end = ibahippinj_end(~isnan(ibahippinj_end));
Cout_ibaseed = C_ret(ibahippinj_seed,:).';
Cin_ibaseed = C_ret(:,ibahippinj_seed);
Cout_ibaseed = Cout_ibaseed(~isnan(ibahippinj_end));
Cin_ibaseed = Cin_ibaseed(~isnan(ibahippinj_end));

figure; 
subplot(1,2,1);
scatter(Cout_ibaseed,ibahippinj_end,'bo','filled'); lsline;
xlabel('Conn. from CA3 Seed'); ylabel('IbaHippInj End')
legend(sprintf('R = %.2f',corr(Cout_ibaseed,ibahippinj_end)));
title('Conn. from Seed vs. IbaHippInj');
set(gca,'FontSize',16);

subplot(1,2,2);
scatter(Cin_ibaseed,ibahippinj_end,'ro','filled'); lsline;
xlabel('Conn. to CA3 Seed'); ylabel('IbaHippInj End')
legend(sprintf('R = %.2f',corr(Cin_ibaseed,ibahippinj_end)));
title('Conn. to Seed vs. IbaHippInj');
set(gca,'FontSize',16);

%% 4. S vs. Aggregation
load([matdir filesep 'aggregation_bias_struct.mat']);
svals_pde = aggregation_bias_struct.PDE.svals;
aggrate_pde = aggregation_bias_struct.PDE.gbratio;
logisticfxn = @(params,x) params(4) + (params(3)*ones(1,length(x))) ./ ...
    (1 + exp(-params(1)*(log10(x) - params(2))));
resfxn = @(params) sum(abs(svals_pde - logisticfxn(params,aggrate_pde)));
% ub = [Inf;Inf]; lb = [0;-Inf];
logistic_opt = fmincon(resfxn,[1;1;1;0]);

svals_dnt_adl = [];
for i = 1:length(aggregation_bias_struct.DNT.adl_names)
    svals = aggregation_bias_struct.DNT.svals_adl{i};
    svals = svals(~isnan(svals));
    svals_dnt_adl = [svals_dnt_adl,svals];
end
aggrate_dnt_adl = zeros(1,length(svals_dnt_adl));
for i = 1:length(aggrate_dnt_adl)
    resfxn_i = @(x) abs(svals_dnt_adl(i) - logisticfxn(logistic_opt,x));
    x0 = (aggrate_pde(end) - aggrate_pde(1))/2; 
    lb = aggrate_pde(1); ub = aggrate_pde(end);
    aggrate_dnt_adl(i) = fmincon(resfxn_i,x0,[],[],[],[],lb,ub);
end

svals_dnt_nadl = [];
for i = 1:length(aggregation_bias_struct.DNT.nadl_names)
    svals = aggregation_bias_struct.DNT.svals_nadl{i};
    svals = svals(~isnan(svals));
    svals_dnt_nadl = [svals_dnt_nadl,svals];
end
aggrate_dnt_nadl = zeros(1,length(svals_dnt_nadl));
for i = 1:length(aggrate_dnt_nadl)
    resfxn_i = @(x) abs(svals_dnt_nadl(i) - logisticfxn(logistic_opt,x));
    x0 = (aggrate_pde(end) - aggrate_pde(1))/2; 
    lb = aggrate_pde(1); ub = aggrate_pde(end);
    aggrate_dnt_nadl(i) = fmincon(resfxn_i,x0,[],[],[],[],lb,ub);
end

figure; hold on;
% scatter(aggrate_pde,svals_pde,'filled');
cmap = lines(2);
xrange = linspace(min(aggrate_pde), max(aggrate_pde));
plot(logisticfxn(logistic_opt,xrange),xrange,'k--','LineWidth',2);
scatter(svals_dnt_adl,aggrate_dnt_adl,50,'o','filled',...
    'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor',cmap(1,:));
scatter(svals_dnt_nadl,aggrate_dnt_nadl,30,'s','filled',...
    'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor',cmap(2,:));
legend({'PDE','ADL','NADL'},'Location','northwest');
xlabel('Agg. Rate'); ylabel('s');
set(gca,'FontSize',16,'yscale','log');

% figure('Units','inches','Position',[0 0 10 9]); hold on;
% cmap = hsv(size(biasmat,2));
% smat = 0.5 - 0.5*biasmat;
% legcell = cell(1,length(fraclist));
% for i = 1:length(fraclist)
%     legcell{i} = sprintf('f = %.1f',fraclist(i));
% end
% for i = 1:size(biasmat,2)
%     gbratio = gammalist(8:24)/beta;
%     plot(gbratio,smat(8:24,i,end),'x-','Color',cmap(i,:),'LineWidth',3);
% end
% xlabel('\gamma/\beta'); ylabel('s'); title('Steady State Bias vs. Aggregation');
% legend(legcell,'Location','southeast');
% set(gca,'FontSize',20,'xscale','log');

%% S1. Functions
    function L = genLplcns(mat)
    
        Dr = sum(mat,2);
        Dc = sum(mat,1);
        small = find(Dr < 0.05 * mean(Dr));
        Dr(small(:)) = 0.05 * mean(Dr);
        small = find(Dc < 0.05 * mean(Dc));
        Dc(small(:)) = 0.05 * mean(Dc);
        Dr = diag(Dr);
        Dc = diag(Dc);
        
        L = eye(size(mat)) - ((Dr^-(1/2)) * mat * (Dc^-(1/2)));
    end

