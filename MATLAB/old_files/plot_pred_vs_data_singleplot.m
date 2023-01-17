function plot_pred_vs_data_singleplot(predicted,data,time_stamps,datasetname)

num_tpts = length(time_stamps);

notnaninds = ~isnan(data(:,1));
data = data(notnaninds,:);
predicted = predicted(notnaninds,:);
data_vec = data(:); predicted_vec = predicted(:);
D = [predicted_vec,data_vec];
% plot_min = min(min(D));
% plot_max = max(max(D));

figure('Units','inches','Position',[0 0 10 9]); hold on; 
cmap = hsv(num_tpts);
ptshape = {'o','d','s','^','p'};
leg = cell(1,num_tpts + 1);

for j = 1:num_tpts
    scatter(predicted(:,j),data(:,j),50,ptshape{j},'MarkerFaceColor',cmap(j,:),'MarkerEdgeColor',cmap(j,:));
    if time_stamps(j) == 1
        leg{j} = sprintf('Time = %d Month',time_stamps(j));
    else
        leg{j} = sprintf('Time = %d Months',time_stamps(j));
    end
end
p = polyfit(predicted_vec,data_vec,1);
xlin = linspace(0,max(predicted_vec),100);
ylin = polyval(p,xlin);
plot(xlin,ylin,'k--','LineWidth',3);

mdl = fitlm(predicted_vec, data_vec);
r2 = mdl.Rsquared.Adjusted;
leg{end} = sprintf('R^2 = %0.2f',r2); legend(leg);
xlabel('Predicted'); ylabel('Data'); title(['Dataset: ' datasetname]);
xlim([0,1.1*max(predicted_vec)]); ylim([0,1.1*max(data_vec)]);
set(gca,'FontSize',20)

end

