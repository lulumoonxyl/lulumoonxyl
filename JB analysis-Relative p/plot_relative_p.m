function t_p_mean=plot_relative_p(t_p_mean,t_p_sem,t_abs_p,x,ii,num,name,ii_num,row,col,ns_removed,edge_removed,both_removed)
%num i just for numbering the fig.#
%name is for naming the subplots within a figure
%ns_removed,edge_removed,both_removed: which data set you are using --'y'
%or []
if ~isempty(both_removed)
    str=', eliminate the first 120s and the edge';
    for i=1:length(name)
        name1(i)=append(name(i),str);
    end
elseif ~isempty(ns_removed)
    str=', eliminate the first 120s';
    for i=1:length(name)
        name1(i)=append(name(i),str);
    end
elseif ~isempty(edge_removed)
    str=', eliminate the edge';
    for i=1:length(name)
        name1(i)=append(name(i),str);
    end
else
    name1=name;
end

a=(x(2)-x(1))/2;
b=x(1)+a;
c=x(end)-a;
x1=b:a*2:c;

figure(num)
hold on ;
subplot(2,3,ii)
patch([x1(1:3)';flipud(x1(1:3)')],[t_p_mean(1:3)-t_p_sem(1:3);flipud(t_p_mean(1:3)+t_p_sem(1:3))],[0 0.45 0.74],'EdgeColor',[1 1 1],'FaceAlpha',0.4)
patch([x1(3:7)';flipud(x1(3:7)')],[t_p_mean(3:7)-t_p_sem(3:7);flipud(t_p_mean(3:7)+t_p_sem(3:7))],[0.22,0.22,0.22],'EdgeColor',[1 1 1],'FaceAlpha',0.4)
patch([x1(7:12)';flipud(x1(7:12)')],[t_p_mean(7:12)-t_p_sem(7:12);flipud(t_p_mean(7:12)+t_p_sem(7:12))],[1 0 0 ],'EdgeColor',[1 1 1],'FaceAlpha',0.4)
patch([x1(12:16)';flipud(x1(12:16)')],[t_p_mean(12:16)-t_p_sem(12:16);flipud(t_p_mean(12:16)+t_p_sem(12:16))],[0.22,0.22,0.22],'EdgeColor',[1 1 1],'FaceAlpha',0.4)
patch([x1(16:18)';flipud(x1(16:18)')],[t_p_mean(16:18)-t_p_sem(16:18);flipud(t_p_mean(16:18)+t_p_sem(16:18))],[0 0.45 0.74],'EdgeColor',[1 1 1],'FaceAlpha',0.4)
hold on;
plot(x1',t_p_mean,'k--','LineWidth',1);
xlabel('Heading direction(degree)')
ylabel('Relative Probability of Orientation')
title(name1(ii));
xlim([x1(1) x1(end)]);
hold off

if ii==ii_num
    %this depends on the number of subplots you have
    figure(num)
    ax1=subplot(row,col,1); ax2=subplot(row,col,2);
    ax3=subplot(row,col,3); ax4=subplot(row,col,4);
    ax5=subplot(row,col,5);
    linkaxes([ax1,ax2,ax3,ax4,ax5],'xy');
end

if ~isempty(t_abs_p)
    figure(1+num)
    hold on ;
    subplot(row,col,ii)
    plot(x1',t_abs_p,'r-','LineWidth',2);
    xlabel('Heading direction(degree)')
    ylabel('Relative Probability of Orientation')
    title(name1(ii));
    xlim([x1(1) x1(end)]);
    hold off;
    if ii==ii_num
        %this depends on the number of subplots you have
        figure(1+num)
        ax1=subplot(row,col,1); ax2=subplot(row,col,2);
        ax3=subplot(row,col,3); ax4=subplot(row,col,4);
        ax5=subplot(row,col,5);
        linkaxes([ax1,ax2,ax3,ax4,ax5],'xy');
    end
end

end