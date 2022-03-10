for ii=1:5
    %% for getting the JAABA and JB data, make sure to check whether their
    %unames are in the same order
    
    name='D:\organized data\jb_result';
    [dat_JB,fields_JB,uname_JB]=get_JB_data(name,ii);
    %% get x data from xspine
    [dat_JB,fields_JB]=get_xy_center(dat_JB,fields_JB);
    %% eliminate the edge
    edge=[5,222;25,250];ns=120; both='y';
    [dat_edge_removed,fields_edge_removed]=eliminate_edge(dat_JB,fields_JB,edge);
    [dat_remove_120s,fields_remove_120s]=remove_the_first_ns(dat_JB,fields_JB,ns);
    [dat_remove_120s_edge,fields_remove_120s_edge]=remove_both(dat_JB,fields_JB,edge,ns,both);
    %% For normal data (nothing removed)
    % calculate the total trakcing time for each animals
    
    [dat_JB,fields_JB]=calculate_tsum(dat_JB,fields_JB,10);%this function also remove tracking time less than 10s
    %the function calculates the heading direction in a clockwise direction
    %may need to change it if the odor is in the right side
    vector=[-1 0];
    [dat_JB,fields_JB]=calculate_heading_direction(dat_JB,fields_JB,vector);
    name={'10^-^3 GA','10^-^2 GA','10^-^1 GA','10^-^5 EA','H2O'};row=2;col=3;ii_num=5;
    % calculate the relative p
    [t_abs,t_p,x]=calculate_relative_p(dat_JB,fields_JB,20);
    [t_p_mean,t_p_sem,t_abs_p]=calculate_mean_of_relative_p(dat_JB,fields_JB,t_p,t_abs);
    t_p_mean=plot_relative_p(t_p_mean,t_p_sem,[],x,ii,1,name,ii_num,row,col,[],[],[]);
    %% for edge removed data
    if ~isempty(dat_edge_removed)
        dat_edge_removed=remove_length_less_than_1(dat_edge_removed,fields_edge_removed)
        [dat_edge_removed,fields_edge_removed]=calculate_tsum(dat_edge_removed,fields_edge_removed,10);
        [dat_edge_removed,fields_edge_removed]=calculate_heading_direction(dat_edge_removed,fields_edge_removed,vector);
        [t_abs_edge_removed,t_p_edge_removed,x_edge_removed]=calculate_relative_p(dat_edge_removed,fields_edge_removed,20);
        [t_p_mean_edge_removed,t_p_sem_edge_removed,t_abs_p_edge_removed]=calculate_mean_of_relative_p(dat_edge_removed,fields_edge_removed,t_p_edge_removed,t_abs_edge_removed)
        t_p_mean_edge_removed=plot_relative_p(t_p_mean_edge_removed,t_p_sem_edge_removed,[],x_edge_removed,ii,10,name,ii_num,row,col,[],'y',[]);
        
    end
    %% for ns removed data
    if ~isempty(dat_remove_120s)
        dat_remove_120s=remove_length_less_than_1(dat_remove_120s,fields_remove_120s);
        [dat_remove_120s,fields_remove_120s]=calculate_tsum(dat_remove_120s,fields_remove_120s,10);
        [dat_remove_120s,fields_remove_120s]=calculate_heading_direction(dat_remove_120s,fields_remove_120s,vector);
        [t_abs_remove_120s,t_p_remove_120s,x_remove_120s]=calculate_relative_p(dat_remove_120s,fields_remove_120s,20);
        [t_p_mean_remove_120s,t_p_sem_remove_120s,t_abs_p_remove_120s]=calculate_mean_of_relative_p(dat_remove_120s,fields_remove_120s,t_p_remove_120s,t_abs_remove_120s);
        t_p_mean_remove_120s=plot_relative_p(t_p_mean_remove_120s,t_p_sem_remove_120s,[],x_remove_120s,ii,20,name,ii_num,row,col,'y',[],[]);
        
    end
    %% for edge and ns removed data
    if ~isempty(dat_remove_120s_edge)
        dat_remove_120s_edge=remove_length_less_than_1(dat_remove_120s_edge,fields_remove_120s_edge);
        [dat_remove_120s_edge,fields_remove_120s_edge]=calculate_tsum(dat_remove_120s_edge,fields_remove_120s_edge,10);
        [dat_remove_120s_edge,fields_remove_120s_edge]=calculate_heading_direction(dat_remove_120s_edge,fields_remove_120s_edge,vector);
        [t_abs_remove_120s_edge,t_p_remove_120s_edge,x_remove_120s_edge]=calculate_relative_p(dat_remove_120s_edge,fields_remove_120s_edge,20);
        [t_p_mean_remove_120s_edge,t_p_sem_remove_120s_edge,t_abs_p_remove_120s_edge]=calculate_mean_of_relative_p(dat_remove_120s_edge,fields_remove_120s_edge,t_p_remove_120s_edge,t_abs_remove_120s_edge);
        t_p_mean_remove_120s_edge=plot_relative_p(t_p_mean_remove_120s_edge,t_p_sem_remove_120s_edge,[],x_remove_120s_edge,ii,30,name,ii_num,row,col,[],[],'y');
    end
    clearvars -except ii
    
end