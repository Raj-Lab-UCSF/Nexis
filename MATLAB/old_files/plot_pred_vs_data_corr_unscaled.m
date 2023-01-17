% Plot prediction vs data 

function plot_pred_vs_data_corr_unscaled(predicted,data,time_stamps)

num_plots = size(data,2);

% % Set min and max points for plot axis
%     notnaninds = ~isnan(data(:,1));
%     data = data(notnaninds,:);
%     predicted = predicted(notnaninds,:);
%     D = [predicted,data];
%     plot_min = min(min(D));
%     plot_max = max(max(D));

figure    
% Plot figure
for j = 1:num_plots
    subplot(1,num_plots,j)
    scatter(predicted(:,j),data(:,j),'o','MarkerFaceColor','k','MarkerEdgeColor','k')
%     axis([plot_min plot_max plot_min plot_max])
    hold on
    h2 = lsline;
    h2.LineWidth = 2; h2.Color = 'r';
%     line([plot_min plot_max],[plot_min plot_max],'LineWidth',2,'Color','r')
    xlabel('Predicted')
    ylabel('Data')
    set(gca,'FontSize',14)
    title(['Time = ' num2str(time_stamps(j))]);
end


