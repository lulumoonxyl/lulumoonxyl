for ii=1:5
    %% for getting the JAABA and JB data, make sure to check whether their
    %unames are in the same order
    [dat_JB,fields_JB,uname_JB]=get_JB_data('D:\jb-results\t88\CS@CS\replication\week1',ii);
    %[dat_JAABA,fields_JAABA,uname_JAABA]=get_JAABA_data('D:\JAABA-results\t88\CS@CS\replication\week1',ii);
    %% eliminate the edge
    [dat_edge_removed,fields_edge_removed]=eliminate_edge(dat_JB,fields_JB);
    %[dat_remove_120s,fields_remove_120s]=remove_the_first_ns(dat_JB,fields_JB,120);
    [dat_remove_120s_edge,fields_remove_120s_edge]=remove_the_first_ns(dat_edge_removed,fields_edge_removed,120);
    %% calculate the total trakcing time for each animals
    %     [dat_JB,fields_JB]=calculate_tsum(dat_JB,fields_JB,10);
    %     [dat_JB,fields_JB]=calculate_heading_direction(dat_JB,fields_JB);
    
    %     dat_edge_removed=remove_length_less_than_1(dat_edge_removed,fields_edge_removed)
    %     [dat_edge_removed,fields_edge_removed]=calculate_tsum(dat_edge_removed,fields_edge_removed,10);
    %     [dat_edge_removed,fields_edge_removed]=calculate_heading_direction(dat_edge_removed,fields_edge_removed);
    
%     dat_remove_120s=remove_length_less_than_1(dat_remove_120s,fields_remove_120s);
%     [dat_remove_120s,fields_remove_120s]=calculate_tsum(dat_remove_120s,fields_remove_120s,10);
%     [dat_remove_120s,fields_remove_120s]=calculate_heading_direction(dat_remove_120s,fields_remove_120s);
%     
        dat_remove_120s_edge=remove_length_less_than_1(dat_remove_120s_edge,fields_remove_120s_edge);
        [dat_remove_120s_edge,fields_remove_120s_edge]=calculate_tsum(dat_remove_120s_edge,fields_remove_120s_edge,10);
        [dat_remove_120s_edge,fields_remove_120s_edge]=calculate_heading_direction(dat_remove_120s_edge,fields_remove_120s_edge);
    
    %% calculate the relative p
    %     [t_abs,t_p,x]=calculate_relative_p(dat_JB,fields_JB,20);
    %     [t_p_mean,t_p_sem,t_abs_p]=calculate_mean_of_relative_p(dat_JB,fields_JB,t_p,t_abs)
    %     t_p_mean=plot_relative_p(t_p_mean,t_p_sem,t_abs_p,x,ii,1)
    
%     [t_abs,t_p,x]=calculate_relative_p(dat_edge_removed,fields_edge_removed,20);
%     [t_p_mean,t_p_sem,t_abs_p]=calculate_mean_of_relative_p(dat_edge_removed,fields_edge_removed,t_p,t_abs)
%     
%     t_p_mean=plot_relative_p(t_p_mean,t_p_sem,t_abs_p,x,ii,1);
    
%     [t_abs,t_p,x]=calculate_relative_p(dat_remove_120s,fields_remove_120s,20);
%     [t_p_mean,t_p_sem,t_abs_p]=calculate_mean_of_relative_p(dat_remove_120s,fields_remove_120s,t_p,t_abs);
%     t_p_mean=plot_relative_p(t_p_mean,t_p_sem,t_abs_p,x,ii,1);
    
        [t_abs,t_p,x]=calculate_relative_p(dat_remove_120s_edge,fields_remove_120s_edge,20);
        [t_p_mean,t_p_sem,t_abs_p]=calculate_mean_of_relative_p(dat_remove_120s_edge,fields_remove_120s_edge,t_p,t_abs);
        t_p_mean=plot_relative_p(t_p_mean,t_p_sem,t_abs_p,x,ii,1);
    
    clearvars -except ii
    
end