function [Rmat,ttest_struct] = CompareGraphMetricPlot(outstruct_,C_,...
    timepoint_flag,matdir_,savenclose_,figdir_)
    
rng(0);
studynames_ = fieldnames(outstruct_);
studynames_(ismember(studynames_,'IbaP301S')) = []; %exclude IbaP301S for too few datapoints
fisher_rtoz = @(r) 0.5*(log(1+r) - log(1-r));

% Define connectomes
C_ret = C_ - diag(diag(C_));
C_ant = C_ret.';

% Connectome properties
% Degree
D_in = sum(C_ret).';
D_out = sum(C_ret,2);
D_mean = (D_in + D_out)/2;

% Normalized Laplacian eigendecomposition, u1
L_ret_norm = eye(size(C_ret)) - diag(D_mean)^(-0.5) * C_ret * diag(D_mean)^(-0.5);
L_ant_norm = eye(size(C_ant)) - diag(D_mean)^(-0.5) * C_ant * diag(D_mean)^(-0.5);
L_ret = L_ret_norm;
L_ant = L_ant_norm;

[v_ret,d_ret] = eig(L_ret);
d_ret = real(diag(d_ret));
[~,sortinds] = sort(d_ret);
v_ret = real(v_ret(:,sortinds));
u1_ret = v_ret(:,1); 

[v_ant,d_ant] = eig(L_ant);
d_ant = real(diag(d_ant));
[~,sortinds] = sort(d_ant);
v_ant= real(v_ant(:,sortinds));
u1_ant = v_ant(:,1);

metric_names = {'\textrm{C, from seed}', '\textrm{C, to seed}', ...
    '\textrm{Out-degree}','\textrm{In-degree}','$v_{1}, L_{ret}$',...
    '$v_{1}, L_{ant}$'};
studynames_ind = [];
Rmat = [];
for i = 1:length(studynames_)
    % Put data in CCF space
    studyname_ = studynames_{i};
    data_i = outstruct_.(studyname_).data;
    data_i = DataToCCF(data_i,studyname_,matdir_);

    % Mean connectivity to seed
    seedreg = outstruct_.(studyname_).seed;
    tptnames = outstruct_.(studyname_).time_stamps;
    if isnan(seedreg)
        seedreg = logical(outstruct_.(studyname_).data(:,1));
        data_i = data_i(:,2:end);
        tptnames = tptnames(2:end);
    end
    if isnumeric(timepoint_flag) && ismember(timepoint_flag,(1:length(tptnames)))
        tptnames = timepoint_flag;
    end
    seedreg_ccf = DataToCCF(seedreg,studyname_,matdir_);
    seedreg_ccf(isnan(seedreg_ccf)) = 0; seedreg_ccf = logical(seedreg_ccf);
    if sum(seedreg) == 1
        conn2seed_out = C_ret(seedreg_ccf,:).'; 
        conn2seed_in = C_ret(:,seedreg_ccf);
    else
        conn2seed_out = mean(C_ret(seedreg_ccf,:)).'; 
        conn2seed_in = mean(C_ret(:,seedreg_ccf),2);
    end
    
    % Remove NaN from all
    naninds = isnan(data_i(:,1));
    data_i(naninds,:) = [];
    D_out_i = D_out; D_out_i(naninds) = [];
    D_in_i = D_in; D_in_i(naninds) = [];
    u1_ret_i = u1_ret; u1_ret_i(naninds) = [];
    u1_ant_i = u1_ant; u1_ant_i(naninds) = [];
    conn2seed_out_i = conn2seed_out; conn2seed_out_i(naninds) = [];
    conn2seed_in_i = conn2seed_in; conn2seed_in_i(naninds) = [];

    metric_cell = {conn2seed_out_i,conn2seed_in_i,D_out_i,D_in_i,u1_ret_i,u1_ant_i};
    % plotcolors = {'b','b','g','g','r','r'};
    % plotshapes = {'o','s','+','x','^','v'};
    
    % Calculate correlations
    Rmat_i = NaN(length(tptnames),length(metric_cell));
    for j = 1:length(metric_cell)
        metric_vals_j = metric_cell{j};
        for k = 1:length(tptnames)
            data_ik = data_i(:,k);
            Rmat_i(k,j) = corr(data_ik,metric_vals_j,'rows','complete');
        end
    end
    Rmat = [Rmat; Rmat_i];
    studyname_i_ind = repmat(i,length(tptnames),1);
    studynames_ind = [studynames_ind; studyname_i_ind];
end
ttest_struct = struct;
zmat = fisher_rtoz(Rmat);
pvals1 = NaN(1,size(Rmat,2));
ttests1 = pvals1;
pvals2 = NaN(1,length(pvals1)/2);
ttests2 = pvals2;
for i = 1:length(pvals1)
    [~,pval1_i,~,stats_i] = ttest(zmat(:,i));
    pvals1(i) = pval1_i; ttests1(i) = stats_i.tstat;
end
ind_ttest = 0;
for i = 1:length(pvals2)
    ind_ttest = ind_ttest + 1;
    z1_i = zmat(:,ind_ttest);
    ind_ttest = ind_ttest + 1;
    z2_i = zmat(:,ind_ttest);
    [~,pval2_i,~,stats_i] = ttest(z1_i,z2_i);
    pvals2(i) = pval2_i; ttests2(i) = stats_i.tstat;
end
ttest_struct.One_Sample.tstat = ttests1;
ttest_struct.One_Sample.pvals = pvals1*6;
ttest_struct.Paired.tstat = ttests2;
ttest_struct.Paired.pvals = pvals2*3;

% Create boxplots
figure('Units','inches','Position',[0 0 13 8.5]); hold on;
% xpos_mat = NaN(size(Rmat)); 
% gbox = xpos_mat; g = studynames_ind;
% coffset1 = 0.05; coffset2 = 0.1; coffset3 = 0.7;
% cmap_boxplot = [[coffset1 coffset2 coffset3]; ...
%     [coffset2 coffset1 coffset3]; [coffset1+0.7 coffset3 coffset2]; ...
%     [coffset2+0.7 coffset3 coffset1]; [coffset3 coffset1 coffset2]; ...
%     [coffset3 coffset2 coffset1]];

cmap_boxplot = [[0.05, 0.40, 1];...
                [0.15, 0.3, 1];...
                [0.50, 0.80, 0.15];...
                [0.50, 0.9, 0.05];...
                [0.75, 0.05, 1];...
                [0.65, 0.15, 1]];
% markerscatter = {'o','s','+','x','^','v'};
% xposscatter = @(y) 0.2 * (2*rand - 1) + y;
% for j = 1:size(Rmat,2)
    % gbox(:,j) = j;
%     for i = 1:size(Rmat,1)
%         xpos_mat(i,j) = xposscatter(j);
%     end
% end
% Rvec = Rmat(:); 
% gbox = gbox(:);
violin(Rmat,'facecolor',cmap_boxplot,'medc',[]);
% b = boxplot(Rvec,gbox,'Colors',cmap_boxplot,'Symbol','');
% set(b,{'linew'},{2});
% for j = 1:size(Rmat,2)
%     gscatter(xpos_mat(:,j),Rmat(:,j),g,cmap_boxplot(j,:),markerscatter{j},7,'off');
% end
plot([0.5,length(metric_names)+0.5],[0,0],'k:','LineWidth',1);
hLegend = findobj(gcf, 'Type', 'Legend'); hLegend.Visible = 'off';
% hLegend.String = {'Mean R'}; hLegend.FontSize = 22; hLegend.Box = 'on';
xticks(1:length(metric_cell)); xlim([0.5,length(metric_names)+0.5]); 
xticklabels(metric_names);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
yplotmax = 0.85; yplotmin = -0.4;
% yplotmax = max(Rvec) + 0.05; yplotmin = min(Rvec) - 0.05;
ylim([yplotmin,yplotmax]);
yticks([-0.4,0,0.4,0.8]);
yticklabels({'-0.4','0','0.4','0.8'});
ylabel("R"); 
set(gca,'FontSize',28,'FontName','Times');

if savenclose_
    print([figdir_ filesep 'NoModelViolins'],'-dtiffn','-r300'); close;
end
end