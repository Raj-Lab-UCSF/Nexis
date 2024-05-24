function [D_out,D_in,u1_ret,u1_ant,conn2seed_out,conn2seed_in] = ...
    DegreeEigenvectorSeedPlots(outstruct_,studyname_,C_,whichplot_,matdir_,...
    savenclose_,figdir_)

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

metric_cell = {conn2seed_out,conn2seed_in,D_out,D_in,u1_ret,u1_ant};
metric_name = {'C, from seed','C, to seed','Out-degree',...
    'In-degree','$v_{1}, L_{ret}$','$v_{1}, L_{ant}$'};
plotcolors = [[0.05, 0.40, 1];...
                [0.15, 0.3, 1];...
                [0.50, 0.80, 0.15];...
                [0.50, 0.9, 0.05];...
                [0.75, 0.05, 1];...
                [0.65, 0.15, 1]];
plotshapes = {'o','s','+','x','^','v'};

% Plotting
if strcmp(whichplot_,'All')
    ylim_plot = [0 max(data_end)];
    studyname_plot = strrep(studyname_,'_',' ');
    figure('Units','inches','Position',[0 0 9 15]); 
    tiledlayout(3,2,'TileSpacing','compact');
    for i = 1:length(metric_cell)
        nexttile; hold on;
        plotdata_i = metric_cell{i};
        xlim_i = [0, max(plotdata_i)];
        scatter(plotdata_i,data_end,50,'MarkerEdgeColor',plotcolors(i,:),...
            'MarkerFaceColor',plotcolors(i,:),'MarkerFaceAlpha',0.5,'Marker',plotshapes{i}); 
        lm = fitlm(plotdata_i,data_end);
        x_lm = linspace(0,1.5*max(plotdata_i),100).'; 
        [y_lm, y_ci] = predict(lm, x_lm);
        plot(x_lm,y_ci(:,1),'k:'); plot(x_lm,y_ci(:,2),'k:'); 
        fill([x_lm; flipud(x_lm)],[y_ci(:,1); flipud(y_ci(:,2))],[1 0 0.25],...
            'EdgeColor','none','FaceAlpha',0.15);
        plot(x_lm,y_lm,'k','LineWidth',2);
        ylim(ylim_plot); yticks([ylim_plot(1), mean(ylim_plot), ylim_plot(2)]);
        yticklabels({'0',num2str(mean(ylim_plot),'%.1f'),num2str(ylim_plot(2),'%.1f')})  
        if i == 3
            ylabel([studyname_plot ' Pathology']);
        end
        xlim(xlim_i); 
        xticks([0, mean(xlim_i), xlim_i(2)]);
        if ismember(i,[1,2])
            xticklabels({'0',num2str(mean(xlim_i)/10000,'%.1f'),num2str(max(xlim_i)/10000,'%.1f')})        
            text(0.95,-0.17,'\times10^{4}','FontSize',20,'FontName','Times','Units','normalized',...
                'Interpreter','tex');
            xlabel(metric_name{i});
        elseif ismember(i,[5,6])
            xticklabels({'0',num2str(mean(xlim_i),'%.2f'),num2str(max(xlim_i),'%.2f')})
            xlabel(metric_name{i},'Interpreter','latex');
        else
            xticklabels({'0',num2str(mean(xlim_i)/100000,'%.1f'),num2str(max(xlim_i)/100000,'%.1f')})        
            text(0.95,-0.17,'\times10^{5}','FontSize',20,'FontName','Times','Units','normalized',...
                'Interpreter','tex');
            xlabel(metric_name{i});
        end
        [corrR,pval] = corr(plotdata_i,data_end);
        pvalstr = [];
        if pval < 0.05
            pvalstr = [pvalstr '*'];
            if pval < 0.01
                pvalstr = [pvalstr '*'];
                if pval < 0.001
                    pvalstr = [pvalstr '*'];
                end
            end
        end
        yoffsets = [0.1,0.1,0.1,0.1,0.9,0.9];
        if ~isempty(pvalstr)        
            text(0.52,yoffsets(i),sprintf('R = %.2f%s',corrR,pvalstr),...
                'FontSize',18,'FontName','Times','Units','normalized',...
                'FontWeight','bold');
        else
            text(0.62,yoffsets(i),sprintf('R = %.2f%s',corrR,pvalstr),...
                'FontSize',18,'FontName','Times','Units','normalized');
        end
        set(gca,'FontSize',18,'FontName','Times','box','on');
    end
elseif strcmp(whichplot_,'C_seed')
    ylim_plot = [0 max(data_end)];
    studyname_plot = strrep(studyname_,'_',' ');
    figure('Units','inches','Position',[0 0 9 4.5]); 
    tiledlayout(1,2,'TileSpacing','compact');
    for i = 1:2
        nexttile; hold on;
        plotdata_i = metric_cell{i};
        xlim_i = [0, max(plotdata_i)];
        scatter(plotdata_i,data_end,50,'MarkerEdgeColor',plotcolors(i,:),...
            'MarkerFaceColor',plotcolors(i,:),'MarkerFaceAlpha',0.5,'Marker',plotshapes{i}); 
        lm = fitlm(plotdata_i,data_end);
        x_lm = linspace(0,1.5*max(plotdata_i),100).'; 
        [y_lm, y_ci] = predict(lm, x_lm);
        plot(x_lm,y_ci(:,1),'k:'); plot(x_lm,y_ci(:,2),'k:'); 
        fill([x_lm; flipud(x_lm)],[y_ci(:,1); flipud(y_ci(:,2))],[1 0 0.25],...
            'EdgeColor','none','FaceAlpha',0.15);
        plot(x_lm,y_lm,'k','LineWidth',2);
        ylim(ylim_plot); yticks([ylim_plot(1), mean(ylim_plot), ylim_plot(2)]);
        yticklabels({'0',num2str(mean(ylim_plot),'%.1f'),num2str(ylim_plot(2),'%.1f')})  
        if i == 1
            ylabel([studyname_plot ' Pathology']);
        end
        xlim(xlim_i); 
        xticks([0, mean(xlim_i), xlim_i(2)]);
        xticklabels({'0',num2str(mean(xlim_i)/10000,'%.1f'),num2str(max(xlim_i)/10000,'%.1f')})        
        text(0.95,-0.18,'\times10^{4}','FontSize',20,'FontName','Times','Units','normalized',...
            'Interpreter','tex');
        xlabel(metric_name{i});
        [corrR,pval] = corr(plotdata_i,data_end);
        pvalstr = [];
        if pval < 0.05
            pvalstr = [pvalstr '*'];
            if pval < 0.01
                pvalstr = [pvalstr '*'];
                if pval < 0.001
                    pvalstr = [pvalstr '*'];
                end
            end
        end
        yoffset = 0.1;
        if ~isempty(pvalstr)        
            text(0.55,yoffset,sprintf('R = %.2f%s',corrR,pvalstr),...
                'FontSize',18,'FontName','Times','Units','normalized',...
                'FontWeight','bold');
        else
            text(0.62,yoffset,sprintf('R = %.2f%s',corrR,pvalstr),...
                'FontSize',18,'FontName','Times','Units','normalized');
        end
        set(gca,'FontSize',20,'FontName','Times','box','on');
    end
    
end

if savenclose_
    print([figdir_ filesep 'NoModelScatterplots_' whichplot_ '_' studyname_],...
        '-dtiffn','-r300'); close;
end

    % function L = genLplcns(mat)
    % 
    %     Dr = sum(mat,2);
    %     Dc = sum(mat,1);
    %     small = find(Dr < 0.05 * mean(Dr));
    %     Dr(small(:)) = 0.05 * mean(Dr);
    %     small = find(Dc < 0.05 * mean(Dc));
    %     Dc(small(:)) = 0.05 * mean(Dc);
    %     Dr = diag(Dr);
    %     Dc = diag(Dc);
    % 
    %     L = eye(size(mat)) - ((Dr^-(1/2)) * mat * (Dc^-(1/2)));
    % end

end