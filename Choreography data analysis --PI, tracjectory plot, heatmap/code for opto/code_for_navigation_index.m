name_of_folder={'D:\choreography-result\optogenetics\49a@chrimson_exp2','D:\choreography-result\optogenetics\42a@chrimson_exp2','D:\choreography-result\optogenetics\attp2@chrimson_exp2'};
for namelength=1:length(name_of_folder)
    for ii=1:3
        % name_of_folder={'D:\choreography-result\CS@CS_replications\week1';'D:\choreography-result\CS@CS_replications\week2';'D:\choreography-result\CS@CS_replications\week3'}
        [dat,fields,uname]=get_data_from_blob_file_opto(name_of_folder,ii,namelength);
      
        %% vertcat all the data
        [dat,xmax,ymax]=vertcat_all_data(dat,fields);
        % grouped data for larvae analysis
        [dat_grouped,fields_grouped]=group_all_data_based_on_AN(dat,fields);
        %% eliminate the edge
        [dat_edge_removed,dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed]=eliminate_edge(dat, dat_grouped,fields,fields_grouped);
        %% eliminate the first 120s
        [dat_remove_120s,dat_grouped_remove_120s,fields_remove_120s,fields_grouped_remove_120s]=remove_the_first_ns(dat,dat_grouped,fields,fields_grouped,120);
        %% eliminate both 120s and edge
        
        [dat_edge_120s_removed,dat_grouped_edge_120s_removed,fields_edge_120s_removed,fields_grouped_edge_120s_removed]=remove_the_first_ns(dat_edge_removed, dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed,120);
        %% calculate the navigation index
        [dat_grouped,fields_grouped]=calculate_navigational_index(dat_grouped,fields_grouped,ii,uname,300);
        [dat_grouped_edge_removed,fields_grouped_edge_removed]=calculate_navigational_index(dat_grouped_edge_removed,fields_grouped_edge_removed,ii,uname,300);
        [dat_grouped_remove_120s,fields_grouped_remove_120s]=calculate_navigational_index(dat_grouped_remove_120s,fields_grouped_remove_120s,ii,uname,[]);
        
        [dat_grouped_edge_120s_removed,fields_grouped_edge_120s_removed]=calculate_navigational_index(dat_grouped_edge_120s_removed,fields_grouped_edge_120s_removed,ii,uname,300);
        
        %% get the mean and sem of PI
        if ii==1&namelength==1
            PI_mean=[];PI_t_mean=[];PI_t_sem=[];PI_sem=[];
            PI_mean_remove_120s=[];PI_t_mean_remove_120s=[];PI_t_sem_remove_120s=[];PI_sem_remove_120s=[];
            PI_mean_edge_removed=[];PI_t_mean_edge_removed=[];PI_sem_edge_removed=[];PI_t_sem_edge_removed=[];
            PI_mean_edge_120s_removed=[];PI_t_mean_edge_120s_removed=[];PI_t_sem_edge_120s_removed=[];PI_sem_edge_120s_removed=[];
        end
        [PI_mean,PI_t_mean,PI_sem,PI_t_sem]=calculate_SEM_Mean_for_PI(dat_grouped,fields_grouped,ii,PI_mean,PI_t_mean,PI_sem,PI_t_sem,namelength);
        [PI_mean_remove_120s,PI_t_mean_remove_120s,PI_sem_remove_120s,PI_t_sem_remove_120s]=calculate_SEM_Mean_for_PI(dat_grouped_remove_120s,fields_grouped_remove_120s,ii,PI_mean_remove_120s,PI_t_mean_remove_120s,PI_sem_remove_120s,PI_t_sem_remove_120s,namelength);
        [PI_mean_edge_removed,PI_t_mean_edge_removed,PI_sem_edge_removed,PI_t_sem_edge_removed]=calculate_SEM_Mean_for_PI(dat_grouped_edge_removed,fields_grouped_edge_removed,ii,PI_mean_edge_removed,PI_t_mean_edge_removed,PI_sem_edge_removed,PI_t_sem_edge_removed,namelength);
        [PI_mean_edge_120s_removed,PI_t_mean_edge_120s_removed,PI_sem_edge_120s_removed,PI_t_sem_edge_120s_removed]=calculate_SEM_Mean_for_PI(dat_grouped_edge_120s_removed,fields_grouped_edge_120s_removed,ii,PI_mean_edge_120s_removed,PI_t_mean_edge_120s_removed,PI_sem_edge_120s_removed,PI_t_sem_edge_120s_removed,namelength);
        
        clearvars -except PI_mean PI_t_mean PI_sem PI_t_sem PI_mean_remove_120s PI_t_mean_remove_120s PI_sem_remove_120s PI_t_sem_remove_120s uname ii PI_mean_edge_removed PI_t_mean_edge_removed PI_sem_edge_removed PI_t_sem_edge_removed PI_mean_edge_120s_removed PI_t_mean_edge_120s_removed PI_sem_edge_120s_removed PI_t_sem_edge_120s_removed name_of_folder namelength
        
    end
end
%% plot the PI of each group
%PI_t_mean and PI_t_sem can be []; the number 0 is just for organizing the
%figures
PI_mean=plot_PI(PI_mean,PI_sem,0,PI_t_mean,PI_t_sem);
PI_mean_remove_120s=plot_PI(PI_mean_remove_120s,PI_sem_remove_120s,1,[],[]);
PI_mean_edge_removed=plot_PI(PI_mean_edge_removed,PI_sem_edge_removed,2,PI_t_mean_edge_removed,PI_t_sem_edge_removed);
PI_mean_edge_120s_removed=plot_PI(PI_mean_edge_120s_removed,PI_sem_edge_120s_removed,3,PI_t_mean_edge_120s_removed,PI_t_sem_edge_120s_removed);
