function ii=plot_relative_p(dat_JB,ii,num,name,ii_num,row,col,str,bins,outdir)
%num i just for numbering the fig.#
%name is for naming the subplots within a figure
%ns_removed,edge_removed,both_removed: which data set you are using --'y'
%or []
x1=-180+bins/2:bins:180-bins/2;
quadrant=[-135 -45 45 135];
color=[0 0.45 0.74;0.22,0.22,0.22;1 0 0;0.22,0.22,0.22;0 0.45 0.74];
for i=1:length(quadrant)
    [M,I]=min(abs(x1-quadrant(i)));
    range(i)=I;
end
range=[1 range length(x1)];
figure(num)
hold on ;
subplot(row,col,ii)
for i=1:length(range)-1
    patch([x1(range(i):range(i+1))';flipud(x1(range(i):range(i+1))')],[dat_JB.p(range(i):range(i+1),1)-dat_JB.p(range(i):range(i+1),2);flipud(dat_JB.p(range(i):range(i+1),1)+dat_JB.p(range(i):range(i+1),2))],color(i,:),'EdgeColor',[1 1 1],'FaceAlpha',0.4)
end
hold on;
plot(x1',dat_JB.p(:,1),'k--','LineWidth',1);
xlabel('Heading direction(degree)')
ylabel('Relative Probability of Orientation')
title(append(name(ii),str));
xlim([x1(1) x1(end)]);
hold off

if ii==ii_num
    %this depends on the number of subplots you have
    figure(num)
    ax1=subplot(row,col,1); ax2=subplot(row,col,2);
    ax3=subplot(row,col,3); ax4=subplot(row,col,4);
    ax5=subplot(row,col,5);
    linkaxes([ax1,ax2,ax3,ax4,ax5],'xy');

    file_dir=fullfile(outdir,'Figures');
    if ~isfolder(file_dir)
        mkdir(file_dir);
    end 
   
    file_name=fullfile(file_dir,append('relative probability of orientation',str,'.fig'));
    if isfile(file_name)
        delete file_name
    end 
    saveas(figure(num),file_name);
end

end