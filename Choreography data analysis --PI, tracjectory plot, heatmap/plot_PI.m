function PI_mean=plot_PI(PI_mean,PI_sem,num,PI_t_mean,PI_t_sem,str,color,name)
%plot the histogram of PI based on each the mean and SEM of each group
%num is just incase we need to plot multiple figures
%str to indicate if you eliminate the edge or not 
if isempty(PI_t_sem)||isempty(PI_t_sem)
    l=1;
else
    l=2;
end
figure(100+num)
subplot(1,l,1);
hold on;
title(str);
hBar = bar(PI_mean);
errorbar(PI_mean,PI_sem,'.k')
hBar.FaceColor='flat';
xlabel('Groups');
ylabel('Preferential Index: v_x/s');
xlim([0 6]);
xticks([0:1:6]);
% yline(0,'r--');
xticklabels(name);

for i=1:length(PI_mean)
    hBar.CData(i,:)=color(i,:);
end 

% %if there is a significant difference
% %use the following command to draw
% %here, 1 and 2 just mean different groups; the # of stars depend on the
% %p-value
% ctr2 = bsxfun(@plus, hBar(1).XData, hBar(1).XOffset');
% plot(ctr2(1:2), [1 1]*index(1,2)*-9, '-k', 'LineWidth',2)
% plot(mean(ctr2(1:2)), index(1,2)*-9.5, '*k','MarkerSize',7)

if l==2
    figure(100+num)
    subplot(1,l,2);
   
    hBar = bar(PI_t_mean);
    for k1 = 1:size(PI_t_mean,2)
        ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: XOffset Is An Undocumented Feature; This Selects The ‘bar’ Centres
        ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    end

    hold on
    errorbar(ctr, ydt,PI_t_sem' , '.k') 
    legend('0-300s','300-600s','600-900s')
    title('Preferential Index Based on Time');
   
    
    xlabel('Groups');
    ylabel('Navigational Index: v_x/s');
    xlim([0 6]);
    xticks([0:1:6]);
    % yline(0,'r--');
    xticklabels(name);
    hold off;
end
end