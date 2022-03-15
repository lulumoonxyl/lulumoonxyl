function tt=calculate_tracking_time(tt,dat_grouped,xmax,ii,ii_num,name,fig_num,row,col,color,str,outdir)

%this function will plot trackingtime vs xbins
x=0:5:round(xmax);

for j=1:length(x)-1
    tt(ii,j)=0;
    for i=1:length(dat_grouped.x)
        idx=find(dat_grouped.x{i,1}>x(j)&dat_grouped.x{i,1}<=x(j+1));
        if isempty(idx)
            continue
        elseif length(idx)==length(dat_grouped.x{i,1})
            tt(ii,j)=tt(ii,j)+dat_grouped.x{i,1}(end)-dat_grouped.x{i,1}(1);
        else
            idx1=find(idx==1);
            idx(idx1)=[];
            tt(ii,j)=tt(ii,j)+sum(dat_grouped.et{i,1}(idx)-dat_grouped.et{i,1}(idx-1));
        end
    end
end
if ii==ii_num

    for i=1:ii_num
        figure(fig_num+3)
        hold on
        plot(x(1:end-1),tt(i,:),'-','LineWidth',2,'Color',color(i,:));
        ylabel('Total Tracking Time (s)')
        xlabel('X position (mm)');
        xlim([0 x(end-1)]);
        str1=regexprep(str,',',' ');
        title(str1);
        axis square;
        legend(name);
        hold off
    end


    file_dir=fullfile(outdir,'Figures');
    if ~isfolder(file_dir)
        mkdir(file_dir);
    end

    filename=fullfile(file_dir,append('tracking time vs x position',str,'.fig'));
    if isfile(filename)
        delete(filename);
    end
        saveas(figure(fig_num+3),filename)
    
end
end

