function lim=plot_heatmap(dat,name,ii,ii_num,row,col,xmax,ymax,fig_num,str,lim,outdir)
%lim is used for caliberation
%this function will count how many times that the larvae appear in
%different x locaitons
x=0:5:round(xmax);
y=0:5:round(ymax);

for a=1:length(x)-1
    for b=1:length(y)-1
        idx=find(dat.x>x(a)&dat.x<=x(a+1)&dat.y>y(b)&dat.y<y(b+1));
        num_dat(length(y)-b,a)=length(idx);
        clear idx;
    end
end

figure(fig_num+2)
hold on
subplot(row,col,ii)

imagesc(num_dat);
xticklabels([x(10),x(20),x(30),x(40)]);
yticklabels([y(end)-y(10),y(end)-y(20),y(end)-y(30),y(end)-y(40)]);
title(append(name(ii),str));
colorbar
lim(ii,:)=caxis;
hold off
if ii==ii_num
    figure(fig_num+2)
    axis_max=max(lim(:,2));
    axis_min=min(lim(:,1));
    for i=1:ii_num
        subplot(row,col,i)
        caxis([axis_min axis_max]);
    end

    file_dir=fullfile(outdir,'Figures');
    if ~isfolder(file_dir)
        mkdir(file_dir);
    end

    filename=fullfile(file_dir,append('heatmap',str,'.fig'));
    if isfile(filename)
        delete(filename);
    end
        saveas(figure(fig_num+2),filename);
    
end

end