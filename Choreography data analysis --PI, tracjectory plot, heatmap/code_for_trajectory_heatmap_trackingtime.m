%% the order of the function cannot be changed bsecause they depend on each other for the input variables
%or you can modify the function
%% this code is for the natural odor
for ii=1:5
     name_of_folder={'D:\choreography-result\CS@CS_replications'};
    [dat,fields,uname]=get_data_from_blob_file(name_of_folder,ii)
    
    %% vertcat all the data
    [dat,xmax,ymax]=vertcat_all_data(dat,fields);
    % grouped data for larvae analysis
    [dat_grouped,fields_grouped]=group_all_data_based_on_AN(dat,fields);
    [dat,dat_grouped,fields,fields_grouped]=generate_centered_data(dat,dat_grouped,fields,fields_grouped);
        %% plot the tracjectory and heatmap and the centered trajectory
    lim=plot_trajectory_heatmap(dat,fields,ii,uname,xmax,ymax,1);
    %
    %     %% to fix the xlim ylim et, edit and use the following function (1 is the figure#,5 is the #of subplot)
    %     % i=fix_figure_title_etc(1,5);
    %     %% calculate and plot the total tracking time vs x bins
    %tt=calculate_tracking_time(dat_grouped,fields,xmax,ii,uname,1);
    
    %% remove the first n seconds (2min=120s)
    %[dat_remove120s,dat_grouped_remove120s,fields_remove120s,fields_grouped_remove120s]=remove_the_first_ns(dat,dat_grouped,fields,fields_grouped,120);
    %% plot the tracjectory and heatmap and the centered trajectory after eliminating the time
    %lim120=plot_trajectory_heatmap(dat_remove120s,fields,ii,uname,xmax,ymax,10);
    %% calculate and plot the total tracking time vs x bins after eliminating the time
    %tt120=calculate_tracking_time(dat_grouped_remove120s,fields_grouped,xmax,ii,uname,10);
    %% eliminate the edge
    %in this case, maynot be even required to use the edge removal function
    %we do not need this for the heatmap, tracjectory or the center
    %trajectory
    
    %[dat_edge_removed,dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed]=eliminate_edge(dat, dat_grouped,fields,fields_grouped);
    %% plot the tracjectory and heatmap and the centered trajectory after eliminating the time and edge
    % edge removal is not suitable for calculating the tracking time and
    % plotting because it just simply cancel out certain x bins
    %     lim120_edge=plot_trajectory_heatmap(dat_edge_removed,fields,ii,uname,xmax,ymax,100);
    %     tt120_edge=calculate_tracking_time(dat_grouped_edge_removed,fields_grouped,xmax,ii,uname,100);
end

