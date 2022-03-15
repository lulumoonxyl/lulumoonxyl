function ii=plot_trajectory(dat,color,name,ii,ii_num,row,col,fig_num,str,outdir)
figure(fig_num);
subplot(row,col,ii)
hold on
scatter1=scatter(dat.x,dat.y,0.1,color(ii,:));
scatter1.MarkerEdgeAlpha=0.05;
title(append(name(ii),str));
axis square
hold off
if ii==ii_num %make sure the axes of all the graphys will be the same
    figure(fig_num)
    ax1=subplot(row,col,1); ax2=subplot(row,col,2);
    ax3=subplot(row,col,3); ax4=subplot(row,col,4);
    ax5=subplot(row,col,5);
    linkaxes([ax1,ax2,ax3,ax4,ax5],'xy');

    clear ax1 ax2 ax3 ax4 ax5

    file_dir=fullfile(outdir,'Figures');
    if ~isfolder(file_dir)
        mkdir(file_dir);
    end

    filename=fullfile(file_dir,append('trjactory',str,'.fig'));
    if isfile(filename)
        delete(filename);
    end
        saveas(figure(fig_num),filename)
   
end
end
