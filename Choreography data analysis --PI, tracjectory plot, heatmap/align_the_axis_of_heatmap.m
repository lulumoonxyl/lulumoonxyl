function lim=align_the_axis_of_heatmap(lim,num,ii_num)
axis_max=max(lim(:,2));
axis_min=min(lim(:,1));
for i=1:ii_num
    figure(num)
    subplot(2,3,i)
    caxis([axis_min axis_max])
end 
end 