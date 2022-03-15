function fig_num=plot_PI_timeseries(PI_timeseries,name,color,time,fig_num,outdir,str)
t_diff=time(2)-time(1);
t1=(time(1)+(t_diff/2):t_diff:time(end)-(t_diff/2))';


for i=1:length(PI_timeseries)
    idx=find(isnan(PI_timeseries{i,1}));
    if ~isempty(idx)
    PI_timeseries{i,1}(idx)=[];
    PI_timeseries{i,2}(idx)=[];
    if i==1
    t1(idx)=[];
    end
    end
    figure(fig_num+6)
    str1=regexprep(str,',',' ');
    title(str1);
    patch([t1;flipud(t1)],[PI_timeseries{i,1}-PI_timeseries{i,2};flipud(PI_timeseries{i,1}+PI_timeseries{i,2})],...
        color(i,:),'EdgeColor',color(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.2);
    hold on

    p(i,1)=plot(t1,PI_timeseries{i,1},'-','LineWidth',1,'Color',color(i,:));
    xlabel('Time(s)');
    ylabel('Preferential Index');
    
    axis tight
    hold off
    if i==length(PI_timeseries)
        legend(p,name);
    end
        file_dir=fullfile(outdir,'Figures');
        if ~isfolder(file_dir)
            mkdir(file_dir);
        end

        filename=fullfile(file_dir,append('Preferential index timeseries',str,'.fig'));
        if isfile(filename)
            delete(filename);
        end
            saveas(figure(fig_num+6),filename)
        
    end

end
