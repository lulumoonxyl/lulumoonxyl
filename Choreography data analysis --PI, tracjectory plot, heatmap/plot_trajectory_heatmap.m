function lim=plot_trajectory_heatmap(dat,fields,ii,uname,xmax,ymax,num)
% this function will just simply plot the trajectory and the heatmap
%for the trajectory
%num is for numbering the figure 
x_pos=find(matches(fields,'x'));
y_pos=find(matches(fields,'y'));

% xc=find(matches(fields,'xcentered'));
% yc=find(matches(fields,'ycentered'));

figure(num);
hold on
subplot(2,3,ii)%this depends on how many groups you have
color=['b','b','b','r','k']; % this will just make different group different color--> natural odor groups
scatter1=scatter(dat{x_pos,1},dat{y_pos,1},0.1,color(ii));
scatter1.MarkerEdgeAlpha=0.05;
title(uname(ii));
xlim([0 xmax]);
ylim([0 ymax]);
axis square
hold off
%% plot centered trajectory
% if (isempty(xc)&isempty(yc))==0
% figure(num+1)
% subplot(2,3,ii)
% hold on 
% scatter2=scatter(dat{xc,1},dat{yc,1},0.1,color(ii));
%scatter2.MarkerEdgeAlpha=0.05;
% title(uname(ii));
% axis square
% hold off
% end 

x=0:5:round(xmax);
y=0:5:round(ymax);
%heatmap
%simply count all the data point in each x and y bins
for a=1:length(x)-1
    for b=1:length(y)-1
        idx=find(dat{x_pos,1}>x(a)&dat{x_pos,1}<=x(a+1)&dat{y_pos}>y(b)&dat{y_pos}<y(b+1));
        num_dat(length(y)-b,a)=length(idx);
        clear idx;
    end
end

figure(num+2)
hold on
subplot(2,3,ii)
title(uname(ii));
imagesc(num_dat);
xticklabels([x(10),x(20),x(30),x(40)]);
yticklabels([y(end)-y(10),y(end)-y(20),y(end)-y(30),y(end)-y(40)]);
colorbar

%output is the variable lim which is the limit of the of the colorbar-->for
%caliberatin later
lim(ii,:)=caxis;
hold off
end
