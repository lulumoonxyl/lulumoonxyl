function PI_mean=plot_PI(PI_mean,PI_sem,PI_t_mean,PI_t_sem,fig_num,str,color,name,ii_num,t,outdir)
%plot the histogram of PI based on each the mean and SEM of each group
%fig_num is just incase we need to plot multiple figures
%str to indicate if you eliminate the edge or not
if isempty(PI_t_sem)||isempty(PI_t_sem)
    l=1;
else
    l=2;
end
figure(fig_num+5)
subplot(1,l,1);
hold on;
str1=regexprep(str,',',' ');
title(str1);
hBar = bar(PI_mean);
errorbar(PI_mean,PI_sem,'.k')
hBar.FaceColor='flat';
xlabel('Groups');
ylabel('Navigational Index: v_x/s');
xlim([0 ii_num+1]);
xticks([0:1:ii_num+1]);
% yline(0,'r--');
name=[' ' name];
xticklabels(name);

for i=1:length(PI_mean)
    hBar.CData(i,:)=color(i,:);
end


if l==2
    time=0:t:900;
    figure(fig_num+5)
    subplot(1,l,2);

    hBar = bar(PI_t_mean,'grouped');
    [ngroups,nbars] = size(PI_t_mean);
    x = nan(nbars,ngroups);
    for i = 1:nbars
        x(i,:) = hBar(i).XEndPoints;
    end
    % Plot the errorbars
    hold on
    errorbar(x',PI_t_mean,PI_t_sem,'k','linestyle','none');
    hold off


    for i=1:length(time)-1
        legends{i,1}=append(num2str(time(i)),'-',num2str(time(i+1)));
    end
    legend(legends);
    title('Navigational Index Based on Time');


    xlabel('Groups');
    ylabel('Navigational Index: v_x/s');
    xlim([0 ii_num+1]);
    xticks([0:1:ii_num+1]);
    % yline(0,'r--');
    xticklabels(name);
    hold off;
end

file_dir=fullfile(outdir,'Figures');
if ~isfolder(file_dir)
    mkdir(file_dir);
end

filename=fullfile(file_dir,append('preferential index',str,'.fig'));
if isfile(filename)
    delete(filename);
end
    saveas(figure(fig_num+5),filename)

end