%% 0. Loading data
rng(0); clear; clc; close all;
matdir = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
figdir = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Figures/TauDirectionality';
output_dir = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Results_Tables_TauDir';
load([matdir filesep 'Connectomes.mat'],'Connectomes');
load([matdir filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct')

%% 1. Setting global properties of model
% Connectomes
C_ret = Connectomes.default; C_ant = C_ret.'; C_nd = (C_ret + C_ant)/2;
C_ret = (C_ret - diag(diag(C_ret)));
C_ant = (C_ant - diag(diag(C_ant)));
C_nd = (C_nd - diag(diag(C_nd)));

% Study names
studynames = fieldnames(mousedata_struct);
studynames(ismember(studynames,'IbaP301S')) = []; % Remove this study 

% Non-default NexIS parameters
bootstrapping = 0;
w_dir = 1;
volcorrect = 1;
use_dataspace = 1;
param_init = [NaN,0,1,0.5];

% Connectome properties
L_ret = sum(C_ret) - C_ret;
L_ant = sum(C_ant) - C_ant;

[v_ret,d_ret] = eig(L_ret); d_ret = abs(diag(d_ret)); v_ret = abs(v_ret);
[dretsort,sortinds] = sort(d_ret); vretsort = v_ret(:,sortinds);

[v_ant,d_ant] = eig(L_ant); d_ant = abs(diag(d_ant)); v_ant = abs(v_ant);
[dantsort,sortinds] = sort(d_ant); vantsort = v_ant(:,sortinds);

%% 2. NexIS:global w/directionality modeling
% fit longitudinally alpha/beta/s (if fit_s)
% fix gamma and alpha for per-timepoint, use LinR

%% 2.1 Longitudinal models
saveoutputs = 1;
filename_out = 'outputs_all_baseline';
outputs_all = struct;
modelnames = {'fit_s','ret','ant','nd'};
for i = 1:length(studynames)
    tablename = [filename_out '_' studynames{i}];
    sumtable = [];
    % if ~isnan(mousedata_struct.(studynames{i}).seed)
    %     gammaval = 1/(sum(mousedata_struct.(studynames{i}).seed));
    % else
    %     gammaval = 1;
    % end
    ub = [Inf,Inf,Inf,1]; ubs = repmat(ub,4,1); ubs(:,end) = [1,1,0,0.5].';
    lb = [0,0,0,0]; lbs = repmat(lb,4,1); lbs(:,end) = [0,1,0,0.5].';
    for j = 1:length(modelnames)
        fprintf('Study %d of %d, Model %s\n',i,length(studynames),modelnames{j})
        outputs = NexIS_global('study',studynames{i},...
                                'w_dir',w_dir,...
                                'volcorrect',volcorrect,...
                                'param_init',param_init,...
                                'ub',ubs(j,:),...
                                'lb',lbs(j,:),...
                                'use_dataspace',use_dataspace,...
                                'bootstrapping',bootstrapping);
        outputs_all.(studynames{i}).(modelnames{j}) = outputs;
        sumtable_i = Output2Table(outputs,0,'null','null');
        sumtable_i.Properties.RowNames{1} = modelnames{j};
        sumtable = [sumtable; sumtable_i];
    end
    if saveoutputs
        writetable(sumtable,[output_dir filesep tablename '.csv'],'WriteRowNames',true)
    end
end
if saveoutputs
    save([output_dir filesep filename_out '.mat'],'outputs_all');
end

%% 2.2 Figures per 2.1
preload = 1;
filename_out = 'outputs_all_revvol';
if preload
    load([output_dir filesep filename_out '.mat'],'outputs_all');
end
CompareDirPlots_deltaR(outputs_all,0);
CompareDirPlots_s(outputs_all,0);
CorrComparePlot(outputs_all,0);

%% 2.3.1 Per-timepoint models, Lin R cost function, fix gamma and alpha
saveoutputs = 1;
outputs_all_tpt = struct;
filename_out = 'outputs_all_tpt_fixgammaalpha_baseline';
modelnames = {'fit_s','ret','ant','nd'};
costfun = 'linr';
% excl_tpts = [[2,3];[1,3];[1,2]];
for i = 1:length(studynames)
    tablename = [filename_out '_' studynames{i}];
    sumtable = [];
    for j = 1:length(modelnames)
        fprintf('Study %d of %d, Model %s\n',i,length(studynames),modelnames{j})
        params_opt = outputs_all.(studynames{i}).(modelnames{j}).nexis_global.Full.param_fit;
        gammaval = params_opt(1); alphaval = params_opt(2);
        ub = [gammaval,alphaval,Inf,1]; ubs = repmat(ub,4,1); ubs(:,end) = [1,1,0,0.5].';
        lb = [gammaval,alphaval,0,0]; lbs = repmat(lb,4,1); lbs(:,end) = [0,1,0,0.5].';
        if strcmp(studynames{i},'Hurtado')
            excl_tpts = [[2,3];[1,3];[1,2]];
        else
            excl_tpts = [2; 1];
        end
        for k = 1:size(excl_tpts,1)
            fprintf('Timepoint %d of %d\n',k,size(excl_tpts,1))
            excl_tpt = excl_tpts(k,:);
            % tpt_str = ['t_' num2str(setdiff(1:3,excl_tpt))];
            if strcmp(studynames{i},'Hurtado')
                tpt_str = ['t_' num2str(setdiff(1:3,excl_tpt))];
            else
                tpt_str = ['t_' num2str(setdiff(1:2,excl_tpt))];
            end
            outputs = NexIS_global('study',studynames{i},...
                                    'w_dir',w_dir,...
                                    'volcorrect',volcorrect,...
                                    'param_init',param_init,...
                                    'ub',ubs(j,:),...
                                    'lb',lbs(j,:),...
                                    'use_dataspace',use_dataspace,...
                                    'bootstrapping',bootstrapping,...
                                    'costfun',costfun,...
                                    'excltpts_costfun',excl_tpt);
            outputs_all_tpt.(studynames{i}).(modelnames{j}).(tpt_str) = outputs;
            sumtable_i = Output2Table(outputs,0,'null','null');
            sumtable_i.Properties.RowNames{1} = [modelnames{j} ', ' tpt_str];
            sumtable_i.Properties.VariableNames{16} = 'R';
            sumtable = [sumtable; sumtable_i];
        end
    end
    if saveoutputs
        writetable(sumtable,[output_dir filesep tablename '.csv'],'WriteRowNames',true)
    end
end
if saveoutputs
    save([output_dir filesep filename_out '.mat'],'outputs_all_tpt');
end

%% 2.4 Figures per 2.3
preload = 1;
filename_out = 'outputs_all_tpt_fixgammaalpha';
if preload
    load([output_dir filesep filename_out '.mat'],'outputs_all_tpt');
end
% CompareDirPlots_deltaR(outputs_all_tpt,1);
% [~,sadl,snadl] = CompareDirPlots_s(outputs_all_tpt,1);
for i = 1:2
    % PerTimepointPlot_sbeta(outputs_all_tpt,i-1);
    DirectionalityVsTimePlot(outputs_all_tpt,i-1)
end
% CorrComparePlot(outputs_all_tpt,1);

%% 4. Figure 2, model
% studyname = 'IbaHippInj';
% rescale = 'log';
% tvec = [0,mousedata_struct.(studyname).time_stamps];
% bvec = 0:0.05:1;
% t = 5;
% tinds = find(ismember(t*bvec,tvec));
% 
% if strcmp(rescale,'log')
%     pathology_sn = log(mousedata_struct.(studyname).data+1);
% end
% seed_sn = mousedata_struct.(studyname).seed;
% pathwseed = [seed_sn, pathology_sn];
% 
% rate_sn = (pathwseed(:,2:end) - pathwseed(:,1:end-1))./...
%     (tvec_data(2:end) - tvec_data(1:end-1));
% 
% seed_sn_ccf = DataToCCF(seed_sn,studyname,matdir);
% model_ant_mat = zeros(length(seed_sn_ccf),length(tvec));
% model_ret_mat = model_ant_mat;
% model_nd_mat = model_ant_mat;
% 
% Rvals_ant = zeros(1,length(bvec)); Rvals_ant_rate = Rvals_ant;
% Rvals_ret = zeros(1,length(bvec)); Rvals_ret_rate = Rvals_ret;
% Rvals_nd = zeros(1,length(tvec_data)); Rvals_nd_rate = Rvals_nd;
%     Rvals_ant(i) = corr(Xi_ant(:,end),pathology_sn(:,end),'rows','complete');
%     Rvals_ant_rate(i) = corr(Xi_ant_rate(:,end),rate_sn(:,end),'rows','complete');
%     model_ant_mat(:,:,i) = Xi_ant;
% 
% Xi_ant = spreadModel(seed_sn,tvec_model,L_ant,bvec);
% if strcmp(rescale,'log')
%     Xi_ant = log(Xi_ant + 1);
% end
% Xi_ant_rate = (Xi_ant(:,tinds(2:end)) - Xi_ant(:,tinds(1:end-1))) ./ ...
%     (tvec_model(tinds(2:end)) - tvec_model(tinds(1:end-1)));
% 
% Xi_ret = spreadModel(seed_sn,tvec_model,L_ret,bvec);
% if strcmp(rescale,'log')
%     Xi_ret = log(Xi_ret + 1);
% end
% Xi_ret_rate = (Xi_ret(:,tinds(2:end)) - Xi_ret(:,tinds(1:end-1))) ./ ...
%     (tvec_model(tinds(2:end)) - tvec_model(tinds(1:end-1)));
% 
% Xi_nd = spreadModel(seed_sn,tvec_model,(L_ret+L_ant)/2,bvec);
% if strcmp(rescale,'log')
%     Xi_nd = log(Xi_nd + 1);
% end
% Xi_nd_rate = (Xi_nd(:,tinds(2:end)) - Xi_nd(:,tinds(1:end-1))) ./ ...
%     (tvec_model(tinds(2:end)) - tvec_model(tinds(1:end-1)));
% 
% 
% 
% 
% % log(X_t + 1) for constructing R-beta*t curves, both pathology and model
% figure; hold on;
% plot(bvec*tvec_model(end),Rvals_ret_rate); plot(bvec*tvec_model(end),Rvals_ant_rate); plot(bvec*tvec_model(end),Rvals_nd_rate); 
% legend({'ret','ant','nd'}); ylabel('R'); xlabel('\beta t'); set(gca,'FontSize',16); title('Rate R')
% 
% figure; hold on;
% plot(bvec*tvec_model(end),Rvals_ret); plot(bvec*tvec_model(end),Rvals_ant); plot(bvec*tvec_model(end),Rvals_nd); 
% legend({'ret','ant','nd'}); ylabel('R'); xlabel('\beta t'); set(gca,'FontSize',16); title('Pathology R')


% %% 3. Conn. from seed
% ibahippinj_seed = logical(seed_all.IbaHippInj);
% ibahippinj_end = taudata_all.IbaHippInj(:,2);
% ibahippinj_end(logical(ibahippinj_seed)) = NaN;
% ibahippinj_end = ibahippinj_end(~isnan(ibahippinj_end));
% Cout_ibaseed = C_ret(ibahippinj_seed,:).';
% Cin_ibaseed = C_ret(:,ibahippinj_seed);
% Cout_ibaseed = Cout_ibaseed(~isnan(ibahippinj_end));
% Cin_ibaseed = Cin_ibaseed(~isnan(ibahippinj_end));
% 
% figure; 
% subplot(1,2,1);
% scatter(Cout_ibaseed,ibahippinj_end,'bo','filled'); lsline;
% xlabel('Conn. from CA3 Seed'); ylabel('IbaHippInj End')
% legend(sprintf('R = %.2f',corr(Cout_ibaseed,ibahippinj_end)));
% title('Conn. from Seed vs. IbaHippInj');
% set(gca,'FontSize',16);
% 
% subplot(1,2,2);
% scatter(Cin_ibaseed,ibahippinj_end,'ro','filled'); lsline;
% xlabel('Conn. to CA3 Seed'); ylabel('IbaHippInj End')
% legend(sprintf('R = %.2f',corr(Cin_ibaseed,ibahippinj_end)));
% title('Conn. to Seed vs. IbaHippInj');
% set(gca,'FontSize',16);
% 
% %% 4. S vs. Aggregation
% load([matdir filesep 'aggregation_bias_struct.mat']);
% svals_pde = aggregation_bias_struct.PDE.svals;
% aggrate_pde = aggregation_bias_struct.PDE.gbratio;
% logisticfxn = @(params,x) params(4) + (params(3)*ones(1,length(x))) ./ ...
%     (1 + exp(-params(1)*(log10(x) - params(2))));
% resfxn = @(params) sum(abs(svals_pde - logisticfxn(params,aggrate_pde)));
% % ub = [Inf;Inf]; lb = [0;-Inf];
% logistic_opt = fmincon(resfxn,[1;1;1;0]);
% 
% svals_dnt_adl = [];
% for i = 1:length(aggregation_bias_struct.DNT.adl_names)
%     svals = aggregation_bias_struct.DNT.svals_adl{i};
%     svals = svals(~isnan(svals));
%     svals_dnt_adl = [svals_dnt_adl,svals];
% end
% aggrate_dnt_adl = zeros(1,length(svals_dnt_adl));
% for i = 1:length(aggrate_dnt_adl)
%     resfxn_i = @(x) abs(svals_dnt_adl(i) - logisticfxn(logistic_opt,x));
%     x0 = (aggrate_pde(end) - aggrate_pde(1))/2; 
%     lb = aggrate_pde(1); ub = aggrate_pde(end);
%     aggrate_dnt_adl(i) = fmincon(resfxn_i,x0,[],[],[],[],lb,ub);
% end
% 
% svals_dnt_nadl = [];
% for i = 1:length(aggregation_bias_struct.DNT.nadl_names)
%     svals = aggregation_bias_struct.DNT.svals_nadl{i};
%     svals = svals(~isnan(svals));
%     svals_dnt_nadl = [svals_dnt_nadl,svals];
% end
% aggrate_dnt_nadl = zeros(1,length(svals_dnt_nadl));
% for i = 1:length(aggrate_dnt_nadl)
%     resfxn_i = @(x) abs(svals_dnt_nadl(i) - logisticfxn(logistic_opt,x));
%     x0 = (aggrate_pde(end) - aggrate_pde(1))/2; 
%     lb = aggrate_pde(1); ub = aggrate_pde(end);
%     aggrate_dnt_nadl(i) = fmincon(resfxn_i,x0,[],[],[],[],lb,ub);
% end
% 
% figure; hold on;
% % scatter(aggrate_pde,svals_pde,'filled');
% cmap = lines(2);
% xrange = linspace(min(aggrate_pde), max(aggrate_pde));
% plot(logisticfxn(logistic_opt,xrange),xrange,'k--','LineWidth',2);
% scatter(svals_dnt_adl,aggrate_dnt_adl,50,'o','filled',...
%     'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor',cmap(1,:));
% scatter(svals_dnt_nadl,aggrate_dnt_nadl,30,'s','filled',...
%     'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor',cmap(2,:));
% legend({'PDE','ADL','NADL'},'Location','northwest');
% xlabel('Agg. Rate'); ylabel('s');
% set(gca,'FontSize',16,'yscale','log');
% 
% % figure('Units','inches','Position',[0 0 10 9]); hold on;
% % cmap = hsv(size(biasmat,2));
% % smat = 0.5 - 0.5*biasmat;
% % legcell = cell(1,length(fraclist));
% % for i = 1:length(fraclist)
% %     legcell{i} = sprintf('f = %.1f',fraclist(i));
% % end
% % for i = 1:size(biasmat,2)
% %     gbratio = gammalist(8:24)/beta;
% %     plot(gbratio,smat(8:24,i,end),'x-','Color',cmap(i,:),'LineWidth',3);
% % end
% % xlabel('\gamma/\beta'); ylabel('s'); title('Steady State Bias vs. Aggregation');
% % legend(legcell,'Location','southeast');
% % set(gca,'FontSize',20,'xscale','log');
% 
% % %% S1. Functions
% %     function L = genLplcns(mat)
% % 
% %         Dr = sum(mat,2);
% %         Dc = sum(mat,1);
% %         small = find(Dr < 0.05 * mean(Dr));
% %         Dr(small(:)) = 0.05 * mean(Dr);
% %         small = find(Dc < 0.05 * mean(Dc));
% %         Dc(small(:)) = 0.05 * mean(Dc);
% %         Dr = diag(Dr);
% %         Dc = diag(Dc);
% % 
% %         L = eye(size(mat)) - ((Dr^-(1/2)) * mat * (Dc^-(1/2)));
% %     end
% % 
