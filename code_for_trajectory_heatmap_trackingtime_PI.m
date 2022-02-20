%% this code is for the natural odor
for ii=1:5
    name_of_folder='D:\choreography-result\CS@CS_replications';
    outdir="C:\Users\feihu\OneDrive - McGill University\matlab\Analysis based on larvae\Fall_2021_Replication\All\Choreography analysis";
    
    [dat,fields,uname]=get_data_from_dat_file(name_of_folder,ii);
    % concatenate all the data
    [dat,xmax,ymax]=vertcat_all_data(dat,fields);
    % grouped data for larvae analysis
    [dat_grouped,fields_grouped]=group_all_data_based_on_AN(dat,fields);
    % get the centered x and y
    [dat,dat_grouped,fields,fields_grouped]=generate_centered_data(dat,dat_grouped,fields,fields_grouped);
    %% plot the tracjectory and heatmap and the centered trajectory (No data is removed)
    %if you don't want to plot, make
    %normal,tt_plot,heatmap,trajjectory,centered_plot to be []-->keeping
    %the both,ns,edge,both_ns,both_edge which will enables you to get the
    %corresponding data:eliminate edge and ns, eliminate ns, and eliminate
    %edge
    color=['b','b','b','r','k'];
    name={'10^-^1 GA','10^-^2 GA','10^-^3 GA','10^-^5 EA','H2O'};
    ns=120; both_ns=120;
    ii_num=5;row=2;col=3;
    normal=[];trajectory=[];centered_plot=[];heatmap=[];tt_plot=[];
    edge=[];both_edge=[5,222;25,250];
    both='y';
    %edge=[];both_edge=[];
    PI_plot=[];
    if ii==1
        ax1=[];ax2=[];lim=[];ax3=[];
        ax7=[];ax8=[];lim_edge_removed=[];ax9=[];
        ax4=[];ax5=[];lim_remove120s=[];ax6=[];
        ax10=[];ax11=[];lim_edge_120s_removed=[];ax12=[];
    end
    
    
    [dat_remove120s,dat_grouped_remove120s,fields_remove120s,fields_grouped_remove120s,...
        dat_edge_removed,dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed,...
        dat_edge_120s_removed,dat_grouped_edge_120s_removed,fields_edge_120s_removed,fields_grouped_edge_120s_removed,...
        ax1,ax2,ax3,lim,tt,ax4,ax5,ax6,lim_remove120s,tt_remove120s,ax7,ax8,ax9,lim_edge_removed,...
        tt_edge_removed,ax10,ax11,ax12,lim_edge_120s_removed,tt_edge_120s_removed]=plot_trajectory_heatmap_tt_based_on_input(dat,...
        fields,dat_grouped,fields_grouped,normal,edge,ns,both,both_ns,both_edge,...
        trajectory,centered_plot,heatmap,tt_plot,ii,ii_num,row,col,xmax,ymax,color,name,...
        ax1,ax2,ax3,lim,ax4,ax5,ax6,lim_remove120s,ax7,ax8,ax9,lim_edge_removed,ax10,ax11,ax12,...
        lim_edge_120s_removed);
    %% Calculate the Preferential index
    %eliminate the larvae that has tracking time less than 10s and calculate the navigation index
    if ii==1
        dat_array={};
        dat_array_edge_removed={};
        dat_array_remove120s={};
        dat_array_edge_120s_removed={};
        PI_mean=[];PI_t_mean=[];PI_t_sem=[];PI_sem=[];
        PI_mean_remove120s=[];PI_t_mean_remove120s=[];PI_t_sem_remove120s=[];PI_sem_remove120s=[];
        PI_mean_edge_removed=[];PI_t_mean_edge_removed=[];PI_sem_edge_removed=[];PI_t_sem_edge_removed=[];
        PI_mean_edge_120s_removed=[];PI_t_mean_edge_120s_removed=[];PI_t_sem_edge_120s_removed=[];PI_sem_edge_120s_removed=[];
    end
    
    if ii==5
        sample_size=[];sample_size_edge_removed=[];sample_size_remove120s=[];sample_size_edge_120s_removed=[];
        
    end
    
    [fields_grouped,dat_grouped]=eliminate_tt(fields_grouped,dat_grouped,10);
    [dat_array,dat_grouped,fields_grouped]=calculate_navigational_index(dat_grouped,fields_grouped,ii,300,dat_array);
    [PI_mean,PI_t_mean,PI_sem,PI_t_sem]=calculate_SEM_Mean_for_PI(dat_grouped,fields_grouped,ii,PI_mean,PI_t_mean,PI_sem,PI_t_sem);
    
    for i=1:length(dat_array)
        sample_size(i,1)=length(dat_array{i,1});
    end
    
    if ~isempty(dat_grouped_edge_removed)
        [fields_grouped_edge_removed,dat_grouped_edge_removed]=eliminate_tt(fields_grouped_edge_removed,dat_grouped_edge_removed,10);
        [dat_array_edge_removed,dat_grouped_edge_removed,fields_grouped_edge_removed]=calculate_navigational_index(dat_grouped_edge_removed,fields_grouped_edge_removed,ii,300,dat_array_edge_removed);
        [PI_mean_edge_removed,PI_t_mean_edge_removed,PI_sem_edge_removed,PI_t_sem_edge_removed]=calculate_SEM_Mean_for_PI(dat_grouped_edge_removed,fields_grouped_edge_removed,ii,PI_mean_edge_removed,PI_t_mean_edge_removed,PI_sem_edge_removed,PI_t_sem_edge_removed);
        for i=1:length(dat_array_edge_removed)
            sample_size_edge_removed(i,1)=length(dat_array_edge_removed{i,1});
        end
    end
    
    if ~isempty(dat_grouped_remove120s)
        [fields_grouped_remove120s,dat_grouped_remove120s]=eliminate_tt(fields_grouped_remove120s,dat_grouped_remove120s,10);
        [dat_array_remove120s,dat_grouped_remove120s,fields_grouped_remove120s]=calculate_navigational_index(dat_grouped_remove120s,fields_grouped_remove120s,ii,[],dat_array_remove120s);
        [PI_mean_remove120s,PI_t_mean_remove120s,PI_sem_remove120s,PI_t_sem_remove120s]=calculate_SEM_Mean_for_PI(dat_grouped_remove120s,fields_grouped_remove120s,ii,PI_mean_remove120s,PI_t_mean_remove120s,PI_sem_remove120s,PI_t_sem_remove120s);
        for i=1:length(dat_array_remove120s)
            sample_size_remove120s(i,1)=length(dat_array_remove120s{i,1});
        end
    end
    
    if ~isempty(dat_grouped_edge_120s_removed)
        [fields_grouped_edge_120s_removed,dat_grouped_edge_120s_removed]=eliminate_tt(fields_grouped_edge_120s_removed,dat_grouped_edge_120s_removed,10);
        [dat_array_edge_120s_removed,dat_grouped_edge_120s_removed,fields_grouped_edge_120s_removed]=calculate_navigational_index(dat_grouped_edge_120s_removed,fields_grouped_edge_120s_removed,ii,300,dat_array_edge_120s_removed);
        [PI_mean_edge_120s_removed,PI_t_mean_edge_120s_removed,PI_sem_edge_120s_removed,PI_t_sem_edge_120s_removed]=calculate_SEM_Mean_for_PI(dat_grouped_edge_120s_removed,fields_grouped_edge_120s_removed,ii,PI_mean_edge_120s_removed,PI_t_mean_edge_120s_removed,PI_sem_edge_120s_removed,PI_t_sem_edge_120s_removed);
        for i=1:length(dat_array_edge_120s_removed)
            sample_size_edge_120s_removed(i,1)=length(dat_array_edge_120s_removed{i,1});
        end
    end
    
    if ii==5
        file_dir=fullfile(outdir,'PI_data');
        if ~isdir(file_dir)
            mkdir(file_dir);
        end
        filename=fullfile(file_dir,'data.mat');
        if isfile(filename)
            delete(filename);
        end
        if ~isfile(filename)
            save( filename, 'dat_array', 'dat_array_edge_removed','dat_array_remove120s','dat_array_edge_120s_removed'...
                ,'PI_mean','PI_t_mean','PI_sem','PI_t_sem','PI_mean_remove120s',...
                'PI_t_mean_remove120s', 'PI_sem_remove120s', 'PI_t_sem_remove120s','PI_mean_edge_removed',...
                'PI_t_mean_edge_removed','PI_sem_edge_removed','PI_t_sem_edge_removed','PI_mean_edge_120s_removed',...
                'PI_t_mean_edge_120s_removed','PI_sem_edge_120s_removed','PI_t_sem_edge_120s_removed','name',...
                'sample_size','sample_size_edge_removed','sample_size_remove120s','sample_size_edge_120s_removed')
        end
        %plot the preferential index and save dat_array for the stat, if
        %the PI_plot input is not empty
        if ~isempty(PI_plot)
            str='Total Preferential Index';
            PI_mean=plot_PI(PI_mean,PI_sem,40,PI_t_mean,PI_t_sem,str);
            if ~isempty(PI_mean_remove120s)
                str1='Total Preferential Index,eliminate the first 120s';
                PI_mean_remove_120s=plot_PI(PI_mean_remove_120s,PI_sem_remove_120s,50,[],[],str1);
            end
            if ~isempty(PI_mean_edge_removed)
                str2='Total Preferential Index,eliminate the edge';
                PI_mean_edge_removed=plot_PI(PI_mean_edge_removed,PI_sem_edge_removed,60,PI_t_mean_edge_removed,PI_t_sem_edge_removed,str2);
            end
            if ~isempty(PI_mean_edge_120s_removed)
                str3='Total Preferential Index,eliminate the edge and the first 120s';
                PI_mean_edge_120s_removed=plot_PI(PI_mean_edge_120s_removed,PI_sem_edge_120s_removed,70,PI_t_mean_edge_120s_removed,PI_t_sem_edge_120s_removed,str3);
                
            end
        end
    end
    
    
    clearvars -except n PI_mean PI_t_mean PI_sem PI_t_sem PI_mean_remove120s...
        PI_t_mean_remove120s PI_sem_remove120s PI_t_sem_remove120s uname ii PI_mean_edge_removed...
        PI_t_mean_edge_removed PI_sem_edge_removed PI_t_sem_edge_removed PI_mean_edge_120s_removed...
        PI_t_mean_edge_120s_removed PI_sem_edge_120s_removed PI_t_sem_edge_120s_removed dat_array...
        dat_array dat_array_edge_removed dat_array_remove120s dat_array_edge_120s_removed...
        ax1 ax2 lim ax3 ax7 ax8 lim_edge_removed ax9 ax4 ax5 lim_remove120s ax6 ...
        ax10 ax11 lim_edge_120s_removed ax12 ii_num both heatmap edge ns
    
end
%% to fix the axis of the heatmap
if ~isempty(heatmap)
    if ~isempty(normal)
        lim=align_the_axis_of_heatmap(lim,3,ii_num);
    end
    
    if ~isempty(ns)
        lim_remove120s=align_the_axis_of_heatmap(lim_remove120s,12,ii_num);
    end
    
    if ~isempty(edge)
        lim_edge_removed=align_the_axis_of_heatmap(lim_edge_removed,22,ii_num);
    end
    
    if ~isempty(both)
        lim_edge_120s_removed=align_the_axis_of_heatmap(lim_edge_120s_removed,32,ii_num);
    end
end