function [ax3,tt]=calculate_tracking_time(dat_grouped,fields_grouped,xmax,ii,ii_num,name,num,ax3,row,col,color)
%num is for numbering the figure
%this function will plot trackingtime vs xbins
x=0:5:round(xmax);
x_pos= find(matches(fields_grouped,'x'));
y_pos=find(matches(fields_grouped,'y'));
et_pos=find(matches(fields_grouped,'et'));

for j=1:length(x)-1
    tt(1,j)=0;
    for i=1:length(dat_grouped{x_pos,1})
        idx=find(dat_grouped{x_pos,1}{i,1}>x(j)&dat_grouped{x_pos,1}{i,1}<=x(j+1));
        if isempty(idx)
            continue
        elseif length(idx)==length(dat_grouped{x_pos,1}{i,1})
            tt(1,j)=tt(1,j)+dat_grouped{x_pos,1}{i,1}(end)-dat_grouped{x_pos,1}{i,1}(1);
        else
            idx1=find(idx==1);
            idx(idx1)=[];
            tt(1,j)=tt(1,j)+sum(dat_grouped{et_pos,1}{i,1}(idx)-dat_grouped{et_pos,1}{i,1}(idx-1));
        end
    end
end

figure(num+3)
if ii==1
    %this depends on how many groups you have
    tiledlayout(col,row);
end
hold on
ax3(ii)=nexttile;
plot(x(1:end-1),tt,'-','LineWidth',2,'Color',color(ii));
title(name(ii))
ylabel('Total Tracking Time (s)')
xlabel('X position (mm)');
xlim([0 x(end-1)]);
axis square;
    if ii==ii_num %change it based on your looping
        linkaxes(ax3,'xy')
    end
    hold off
end