function [dat_remove120s,dat_grouped_remove120s,fields_remove120s,fields_grouped_remove120s,...
    dat_edge_removed,dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed,...
    dat_edge_120s_removed,dat_grouped_edge_120s_removed,fields_edge_120s_removed,fields_grouped_edge_120s_removed,...
    ax1,ax2,ax3,lim,tt,ax4,ax5,ax6,lim_remove120s,tt_remove120s,ax7,ax8,ax9,lim_edge_removed,...
    tt_edge_removed,ax10,ax11,ax12,lim_edge_120s_removed,tt_edge_120s_removed]=plot_trajectory_heatmap_tt_based_on_input(dat,...
    fields,dat_grouped,fields_grouped,normal,edge,ns,both,both_ns,both_edge,...
    trajectory,centered_plot,heatmap,tt_plot,ii,ii_num,row,col,xmax,ymax,color,name,...
    ax1,ax2,ax3,lim,ax4,ax5,ax6,lim_remove120s,ax7,ax8,ax9,lim_edge_removed,ax10,ax11,ax12,...
    lim_edge_120s_removed)
%normal,edge,ns,both: indicates whether you want to remove the edge-->'y' or []
%usually both_ns and ns=120s-->remove the first 2min out of 15 min
%both_edge and edge will be an array indicate the x and y that you want to
%remove-->[xmin,xmax;ymin;ymax]
%both will just be 'y' or []
%trajectory,centered_plot,heatmap,tt_plot:plot these or not-->'y' or []
%ax..is the axes of the plots, lim is the axis of the heatmap
%ii_num:number of groups you have
%dat..,fields..are the corresponding data and fieldnames
%tt is the tracking time vs x
if ~isempty(normal)
num=1;

if ~isempty(trajectory)|~isempty(centered_plot)|~isempty(heatmap)
    [ax1,ax2,lim]=plot_trajectory_heatmap(dat,fields,ii,ii_num,name,color,xmax,ymax,num,col,row,trajectory,centered_plot,heatmap,ax1,ax2,lim)
else
    ax1=[];ax2=[];lim=[];
end
% calculate and plot the total tracking time vs x bins
if ~isempty(tt_plot)
    [ax3,tt]=calculate_tracking_time(dat_grouped,fields_grouped,xmax,ii,ii_num,name,num,ax3,row,col,color);
else
    tt=[];ax3=[];
end
else 
     tt=[];ax3=[]; ax1=[];ax2=[];lim=[];
end 
%% remove the first n seconds (2min=120s)
if ~isempty(ns)
    [dat_remove120s,dat_grouped_remove120s,fields_remove120s,fields_grouped_remove120s]=remove_the_first_ns(dat,dat_grouped,fields,fields_grouped,ns);
    
    str=', eliminate the first 120s';
    for i=1:length(name)
        name1(i)=append(name(i),str);
    end
    % plot the tracjectory and heatmap and the centered trajectory after eliminating the time
    num_remove120s=10;
    if ~isempty(trajectory)|~isempty(centered_plot)|~isempty(heatmap)
        [ax4,ax5,lim_remove120s]=plot_trajectory_heatmap(dat_remove120s,fields_remove120s,ii,ii_num,name1,color,xmax,ymax,num_remove120s,col,row,trajectory,centered_plot,heatmap,ax4,ax5,lim_remove120s)
    else
          ax4=[];ax5=[];lim_remove120s=[];
    end
    
    if ~isempty(tt_plot)
        [ax6,tt_remove120s]=calculate_tracking_time(dat_grouped_remove120s,fields_grouped_remove120s,xmax,ii,ii_num,name1,num_remove120s,ax6,row,col,color);
    else
        tt_remove120s=[]; ax6=[];
    end
    
else
    ax4=[];ax5=[];lim_remove120s=[];
    tt_remove120s=[]; ax6=[];
    dat_remove120s=[];dat_grouped_remove120s=[];fields_remove120s=[];fields_grouped_remove120s=[];
end
%% eliminate the edge
if ~isempty(edge)
    [dat_edge_removed,dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed]=eliminate_edge(dat, dat_grouped,fields,fields_grouped,edge);
    %plot the tracjectory and heatmap and the centered trajectory after eliminating the time and edge
    str=', eliminate the edge';
    for i=1:length(name)
        name2(i)=append(name(i),str);
    end
    % plot the tracjectory and heatmap and the centered trajectory after eliminating the time
    num_edge_removed=20;
    
    if ~isempty(trajectory)|~isempty(centered_plot)|~isempty(heatmap)
        [ax7,ax8,lim_edge_removed]=plot_trajectory_heatmap(dat_edge_removed,fields_edge_removed,ii,ii_num,name2,color,xmax,ymax,num_edge_removed,col,row,trajectory,centered_plot,heatmap,ax7,ax8,lim_edge_removed)
    else
        ax7=[]; ax8=[]; lim_edge_removed=[];
    end
    
    if ~isempty(tt_plot)
        [ax9,tt_edge_removed]=calculate_tracking_time(dat_grouped_edge_removed,fields_grouped_edge_removed,xmax,ii,ii_num,name2,num_edge_removed,ax9,row,col,color);
    else
        tt_edge_removed=[];ax9=[];
    end
else
    tt_edge_removed=[];ax9=[]; ax7=[]; ax8=[]; lim_edge_removed=[];
    dat_edge_removed=[]; dat_grouped_edge_removed=[]; fields_edge_removed=[]; fields_grouped_edge_removed=[];
end
%% eliminate both 120s and edge
if ~isempty(both)
    if ~isempty (edge)
        [dat_edge_120s_removed,dat_grouped_edge_120s_removed,fields_edge_120s_removed,fields_grouped_edge_120s_removed]=remove_the_first_ns(dat_edge_removed, dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed,both_ns);
        
    elseif ~isempty(ns)
        [dat_edge_120s_removed,dat_grouped_edge_120s_removed,fields_edge_120s_removed,fields_grouped_edge_120s_removed]=eliminate_edge(dat_remove120s,dat_grouped_remove120s,fields_remove120s,fields_grouped_remove120s,both_edge);
    else
         [dat_edge_removed,dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed]=eliminate_edge(dat, dat_grouped,fields,fields_grouped,both_edge);
         [dat_edge_120s_removed,dat_grouped_edge_120s_removed,fields_edge_120s_removed,fields_grouped_edge_120s_removed]=remove_the_first_ns(dat_edge_removed, dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed,both_ns);
    end
    str=', eliminate the first 120s and the edge';
    for i=1:length(name)
        name3(i)=append(name(i),str);
    end
    num_edge_120s_removed=30;
    if ~isempty(trajectory)|~isempty(centered_plot)|~isempty(heatmap)
        [ax10,ax11,lim_edge_120s_removed]=plot_trajectory_heatmap(dat_edge_120s_removed,fields_edge_120s_removed,ii,ii_num,name3,color,xmax,ymax,num_edge_120s_removed,col,row,trajectory,centered_plot,heatmap,ax10,ax11,lim_edge_120s_removed)
    else
        ax10=[]; ax11=[]; lim_edge_120s_removed=[];
    end
    
    if ~isempty(tt_plot)
        [ax12,tt_edge_120s_removed]=calculate_tracking_time(dat_grouped_edge_120s_removed,fields_grouped_edge_120s_removed,xmax,ii,ii_num,name3,num_edge_120s_removed,ax12,row,col,color);
    else
        tt_edge_120s_removed=[];ax12=[];
    end
else
    dat_edge_120s_removed=[];
    dat_grouped_edge_120s_removed=[];
    fields_edge_120s_removed=[];
    fields_grouped_edge_120s_removed=[];
    tt_edge_120s_removed=[];ax12=[];
    ax10=[]; ax11=[]; lim_edge_120s_removed=[];
end
end