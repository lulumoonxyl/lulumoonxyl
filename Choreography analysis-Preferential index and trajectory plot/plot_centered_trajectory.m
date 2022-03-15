function ii=plot_centered_trajectory(dat,color,name,ii,ii_num,row,col,fig_num,str,outdir)
figure(fig_num+1)
subplot(row,col,ii)
hold on
scatter2=scatter(dat.xcentered,dat.ycentered,0.1,color(ii,:));
scatter2.MarkerEdgeAlpha=0.05;
title(append(name(ii),str));
axis square
hold off
if ii==ii_num %change it based on your looping
    figure(fig_num+1)
    ax1=subplot(row,col,1); ax2=subplot(row,col,2);
    ax3=subplot(row,col,3); ax4=subplot(row,col,4);
    ax5=subplot(row,col,5);
    linkaxes([ax1,ax2,ax3,ax4,ax5],'xy');
    clear ax1 ax2 ax3 ax4 ax5
    file_dir=fullfile(outdir,'Figures');
    if ~isfolder(file_dir)
        mkdir(file_dir);
    end

    filename=fullfile(file_dir,append('centered trjactory',str,'.fig'));
    if isfile(filename)
        delete(filename);
    end
        saveas(figure(fig_num+1),filename)
  
end
end