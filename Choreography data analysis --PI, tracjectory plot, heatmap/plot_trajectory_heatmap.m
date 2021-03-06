function lim=plot_trajectory_heatmap(dat,fields,ii,ii_num,name,color,xmax,ymax,num,col,row,trajectory,centered_plot,heatmap,lim)
% this function will just simply plot the trajectory and the heatmap
%for the trajectory
%num is for numbering the figure
%tell the function to plot the centered plot or not:'y' or []
%trajectory:'y' or []
%heatmap:'y' or []
%ii_num:number of the final loop, 5 in this case
%ax1,ax2,lim are the one that are generated by the plots
x_pos=find(matches(fields,'x'));
y_pos=find(matches(fields,'y'));
%% Plot trajectory
if ~isempty(trajectory)
    
    
    figure(num);
    subplot(row,col,ii)
    hold on
    
    scatter1=scatter(dat{x_pos,1},dat{y_pos,1},0.1,color(ii));
    scatter1.MarkerEdgeAlpha=0.05;
    title(name(ii));
    axis square
    hold on
    if ii==ii_num %change it based on your looping
        figure(num)
        ax1=subplot(row,col,1); ax2=subplot(row,col,2);
        ax3=subplot(row,col,3); ax4=subplot(row,col,4);
        ax5=subplot(row,col,5);
        linkaxes([ax1,ax2,ax3,ax4,ax5],'xy');
        clear ax1 ax2 ax3 ax4 ax5
    end  
end
%% Plot centered trajectories
if ~isempty(centered_plot)
    xc=find(matches(fields,'xcentered'));
    yc=find(matches(fields,'ycentered'));
    
    figure(num+1)
    subplot(row,col,ii)
    hold on
    
    scatter2=scatter(dat{xc,1},dat{yc,1},0.1,color(ii));
    scatter2.MarkerEdgeAlpha=0.05;
    title(name(ii));
    axis square
    hold off
    if ii==ii_num %change it based on your looping
        figure(num+1)
        ax1=subplot(row,col,1); ax2=subplot(row,col,2);
        ax3=subplot(row,col,3); ax4=subplot(row,col,4);
        ax5=subplot(row,col,5);
        linkaxes([ax1,ax2,ax3,ax4,ax5],'xy');
        clear ax1 ax2 ax3 ax4 ax5
    end
    
end
%% Plot heatmap
%simply count all the data point in each x and y bins
if ~isempty(heatmap)
    x=0:5:round(xmax);
    y=0:5:round(ymax);
    
    for a=1:length(x)-1
        for b=1:length(y)-1
            idx=find(dat{x_pos,1}>x(a)&dat{x_pos,1}<=x(a+1)&dat{y_pos}>y(b)&dat{y_pos}<y(b+1));
            num_dat(length(y)-b,a)=length(idx);
            clear idx;
        end
    end
    
    figure(num+2)
    hold on
    subplot(row,col,ii)
    
    imagesc(num_dat);
    xticklabels([x(10),x(20),x(30),x(40)]);
    yticklabels([y(end)-y(10),y(end)-y(20),y(end)-y(30),y(end)-y(40)]);
    title(name(ii));
    colorbar
    lim(ii,:)=caxis;
    hold off
    if ii==5
        figure(num+2)
        axis_max=max(lim(:,2));
        axis_min=min(lim(:,1));
        
        for i=1:ii_num
            
            subplot(row,col,i)
            caxis([axis_min axis_max])
        end
    end
    
else
    lim(ii,:)=[];
end
end
