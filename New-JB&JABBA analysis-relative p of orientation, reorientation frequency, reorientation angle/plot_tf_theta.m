function ii_num=plot_tf_theta(x,data,fig_num,ii_num,name,color,str,deg,xbins,large_or_small,outdir)
%% plot the tf,theta,large turning probability

    x1=x(1)+xbins/2:xbins:x(end)-xbins/2;
    x1=x1';
    if isequal(large_or_small,'large')
    xaxis_label={'Reorientation frequency (s^-^1)',append('Reorientation frequency (>=',num2str(deg),',s^-^1 )'),'Large turning/all turning',...
        'Absolute reorientation angle (deg)'};
    else 
    xaxis_label={'Reorientation frequency (s^-^1)',append('Reorientation frequency (<',num2str(deg),',s^-^1 )'),'Large turning/all turning',...
        'Absolute reorientation angle (deg)'};
    end 
    
    fn=fieldnames(data);
    
    figure(fig_num)
    for l=1:length(fn)
        for i=1:ii_num
            subplot(2,2,l)
            patch([x1;flipud(x1)],[data.(fn{l})(:,2*i-1)-data.(fn{l})(:,2*i);flipud(data.(fn{l})(:,2*i-1)+data.(fn{l})(:,2*i))],...
                color(i,:),'EdgeColor',color(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.2);
            hold on
            p(i,1)=plot(x1,data.(fn{l})(:,2*i-1),'-','LineWidth',1.5,'Color',color(i,:));
            ylabel(xaxis_label(l));
            title(str);
            xlabel('X (mm)');
            xlim([x1(1) x1(end)])
            if i==ii_num
                legend(p,name);
             end
        end
    end
    
     file_dir=fullfile(outdir,'Figures');
    if ~isfolder(file_dir)
        mkdir(file_dir);
    end 
   
    file_name=fullfile(file_dir,append('Reorientation frequency and cngle_',str,'.fig'));
    if isfile(file_name)
        delete file_name
    end 
    saveas(figure(num),file_name);
    clear idx
end