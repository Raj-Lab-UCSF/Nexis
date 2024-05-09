function [D_out,D_in,u1_ret,u1_ant,conn2seed_out,conn2seed_in] = ...
    DegreeEigenvectorSeedPlots(outstruct_,studyname_,C_,matdir_,...
    savenclose_,figdirectory)

% Define connectomes
C_ = C_ - diag(diag(C_));

% Log transform
% C_max = max(nonzeros(C_(:))); C_min = min(nonzeros(C_(:)));
% C_ = log10(C_); C_max = log10(C_max); C_min = log10(C_min);
% C_(C_ == -Inf) = NaN;
% C_ = (C_max - C_) / (C_max - C_min);
% C_(isnan(C_)) = 0;

C_ret = C_;
C_ant = C_ret.';

% Connectome properties
% Degree
D_in = sum(C_ret).';
D_out = sum(C_ret,2);
D_mean = (D_in + D_out)/2;

% Normalized Laplacian eigendecomposition, u2
% L_ret = diag(D_in) - C_ret;
% L_ant = diag(D_out) - C_ant;
% L_ret = genLplcns(C_ret);
% L_ant = genLplcns(C_ant);
% L_ret_norm = eye(size(C_ret)) - diag(D_out)^(-0.5) * C_ret * diag(D_in)^(-0.5);
% L_ant_norm = eye(size(C_ant)) - diag(D_in)^(-0.5) * C_ant * diag(D_out)^(-0.5);
L_ret_norm = eye(size(C_ret)) - diag(D_mean)^(-0.5) * C_ret * diag(D_mean)^(-0.5);
L_ant_norm = eye(size(C_ant)) - diag(D_mean)^(-0.5) * C_ant * diag(D_mean)^(-0.5);
L_ret = L_ret_norm;
L_ant = L_ant_norm;

[v_ret,d_ret] = eig(L_ret);
d_ret = real(diag(d_ret));
[~,sortinds] = sort(d_ret);
v_ret = abs(v_ret(:,sortinds));
u1_ret = v_ret(:,1); 

[v_ant,d_ant] = eig(L_ant);
d_ant = real(diag(d_ant));
[~,sortinds] = sort(d_ant);
v_ant= abs(v_ant(:,sortinds));
u1_ant = v_ant(:,1);

% Mean connectivity to seed
seedreg = outstruct_.(studyname_).seed;
if isnan(seedreg)
    seedreg = logical(outstruct_.(studyname_).data(:,1));
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

% Put data in CCF space
data_end = outstruct_.(studyname_).data(:,end);
data_end = DataToCCF(data_end,studyname_,matdir_);

% Remove NaN from all
naninds = isnan(data_end);
data_end(naninds) = [];
D_out(naninds) = [];
D_in(naninds) = [];
u1_ret(naninds) = [];
u1_ant(naninds) = [];
conn2seed_out(naninds) = [];
conn2seed_in(naninds) = [];

metric_cell = {D_out,D_in,u1_ret,u1_ant,conn2seed_out,conn2seed_in};
metric_name = {'Out-degree','In-degree','u_1, L', 'u_1, L^T',...
    'C_s_e_e_d, outgoing', 'C_s_e_e_d, incoming'};
plotcolors = {'r','b','r','b','r','b'};
plotshapes = {'o','o','s','s','d','d'};

% Plotting
ylim_plot = [0 max(data_end)];
studyname_plot = strrep(studyname_,'_',' ');
figure('Units','inches','Position',[0 0 25 5]); 
tiledlayout(1,length(metric_name),'TileSpacing','compact','Padding','tight');
for i = 1:length(metric_cell)
    nexttile;
    plotdata_i = metric_cell{i};
    xlim_i = [min(plotdata_i), max(plotdata_i)];
    scatter(plotdata_i,data_end,[plotcolors{i} plotshapes{i}],'filled'); 
    l = lsline; l.LineWidth = 2; l.Color = 'k';
    ylim(ylim_plot); yticks([ylim_plot(1), mean(ylim_plot), ylim_plot(2)]);
    if i == 1
        ylabel([studyname_plot ' Pathology']);
    end
    xlim(xlim_i); 
    xticks([xlim_i(1), mean(xlim_i), xlim_i(2)]);
    if ~ismember(i,[3,4])
        xtickformat('%.1f');
    else
        xtickformat('%.2f');
    end
    xlabel(metric_name{i});
    text(0.6,0.1,sprintf('R = %.2f',corr(plotdata_i,data_end)),...
        'FontSize',20,'FontName','Times','Units','normalized');
    set(gca,'FontSize',20,'FontName','Times');
end

if savenclose_
    print([figdirectory filesep 'NoModelScatterplots'],'-dtiffn','-r300'); close;
end

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

end