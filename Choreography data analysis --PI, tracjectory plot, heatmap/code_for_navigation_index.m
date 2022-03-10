
for ii=1:5
   
    %name_of_folder={'D:\choreography-result\CS@CS_replications\week4(Yiran)'};
    [dat,fields,uname]=get_data_from_dat_file(name,ii)
    
    %% vertcat all the data
    [dat,xmax,ymax]=vertcat_all_data(dat,fields);
    % grouped data for larvae analysis
    [dat_grouped,fields_grouped]=group_all_data_based_on_AN(dat,fields);
    %% eliminate the edge
    %     [dat_edge_removed,dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed]=eliminate_edge(dat, dat_grouped,fields,fields_grouped);
    %     %% eliminate the first 120s
    %     [dat_remove_120s,dat_grouped_remove_120s,fields_remove_120s,fields_grouped_remove_120s]=remove_the_first_ns(dat,dat_grouped,fields,fields_grouped,120);
    %     %% eliminate both 120s and edge
    %
    %     [dat_edge_120s_removed,dat_grouped_edge_120s_removed,fields_edge_120s_removed,fields_grouped_edge_120s_removed]=remove_the_first_ns(dat_edge_removed, dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed,120);
    %% eliminate the larvae that has tracking time less than 10s
    [fields_grouped,dat_grouped]=eliminate_tt(fields_grouped,dat_grouped,10)
    %% calculate the navigation index
    if ii==1
        dat_array={};
        dat_array_edge_removed={};
        dat_array_remove_120s={};
        dat_array_edge_120s_removed={};
    end
    
    [dat_array,dat_grouped,fields_grouped]=calculate_navigational_index(dat_grouped,fields_grouped,ii,uname,300,dat_array);
   
    if ii==1;
       n=[];
   end 
    n=get_sample_size(fields_grouped,dat_grouped,ii,n)
    %     [dat_array_edge_removed,dat_grouped_edge_removed,fields_grouped_edge_removed]=calculate_navigational_index(dat_grouped_edge_removed,fields_grouped_edge_removed,ii,uname,300,dat_array_edge_removed);
    %     [dat_array_remove_120s,dat_grouped_remove_120s,fields_grouped_remove_120s]=calculate_navigational_index(dat_grouped_remove_120s,fields_grouped_remove_120s,ii,uname,[],dat_array_remove_120s);
    %
    %     [dat_array_edge_120s_removed,dat_grouped_edge_120s_removed,fields_grouped_edge_120s_removed]=calculate_navigational_index(dat_grouped_edge_120s_removed,fields_grouped_edge_120s_removed,ii,uname,300,dat_array_edge_120s_removed);
    %
    %% get the mean and sem of PI
        if ii==1
            PI_mean=[];PI_t_mean=[];PI_t_sem=[];PI_sem=[];
            PI_mean_remove_120s=[];PI_t_mean_remove_120s=[];PI_t_sem_remove_120s=[];PI_sem_remove_120s=[];
            PI_mean_edge_removed=[];PI_t_mean_edge_removed=[];PI_sem_edge_removed=[];PI_t_sem_edge_removed=[];
            PI_mean_edge_120s_removed=[];PI_t_mean_edge_120s_removed=[];PI_t_sem_edge_120s_removed=[];PI_sem_edge_120s_removed=[];
        end
         [PI_mean,PI_t_mean,PI_sem,PI_t_sem]=calculate_SEM_Mean_for_PI(dat_grouped,fields_grouped,ii,PI_mean,PI_t_mean,PI_sem,PI_t_sem);
%          [PI_mean_remove_120s,PI_t_mean_remove_120s,PI_sem_remove_120s,PI_t_sem_remove_120s]=calculate_SEM_Mean_for_PI(dat_grouped_remove_120s,fields_grouped_remove_120s,ii,PI_mean_remove_120s,PI_t_mean_remove_120s,PI_sem_remove_120s,PI_t_sem_remove_120s);
%         [PI_mean_edge_removed,PI_t_mean_edge_removed,PI_sem_edge_removed,PI_t_sem_edge_removed]=calculate_SEM_Mean_for_PI(dat_grouped_edge_removed,fields_grouped_edge_removed,ii,PI_mean_edge_removed,PI_t_mean_edge_removed,PI_sem_edge_removed,PI_t_sem_edge_removed);
%         [PI_mean_edge_120s_removed,PI_t_mean_edge_120s_removed,PI_sem_edge_120s_removed,PI_t_sem_edge_120s_removed]=calculate_SEM_Mean_for_PI(dat_grouped_edge_120s_removed,fields_grouped_edge_120s_removed,ii,PI_mean_edge_120s_removed,PI_t_mean_edge_120s_removed,PI_sem_edge_120s_removed,PI_t_sem_edge_120s_removed);


    clearvars -except n PI_mean PI_t_mean PI_sem PI_t_sem PI_mean_remove_120s PI_t_mean_remove_120s PI_sem_remove_120s PI_t_sem_remove_120s uname ii PI_mean_edge_removed PI_t_mean_edge_removed PI_sem_edge_removed PI_t_sem_edge_removed PI_mean_edge_120s_removed PI_t_mean_edge_120s_removed PI_sem_edge_120s_removed PI_t_sem_edge_120s_removed dat_array
    
end

% %% plot the PI of each group
% %PI_t_mean and PI_t_sem can be []; the number 0 is just for organizing the
% %figures
PI_mean=plot_PI(PI_mean,PI_sem,0,PI_t_mean,PI_t_sem);
% PI_mean_remove_120s=plot_PI(PI_mean_remove_120s,PI_sem_remove_120s,1,[],[]);
% PI_mean_edge_removed=plot_PI(PI_mean_edge_removed,PI_sem_edge_removed,2,PI_t_mean_edge_removed,PI_t_sem_edge_removed);
% PI_mean_edge_120s_removed=plot_PI(PI_mean_edge_120s_removed,PI_sem_edge_120s_removed,3,PI_t_mean_edge_120s_removed,PI_t_sem_edge_120s_removed);
