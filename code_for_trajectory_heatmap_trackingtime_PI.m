% this code is for the natural odor
ii_num=5;%this is the number of conditions/subfolders
[ns,edge,both_ns,both_edge,plot_input]=get_input_from_user(ii_num);

for ii=1:ii_num
    %% identify the input and output folder
    name_of_folder='D:\choreography-result\choreography_result';
    outdir="C:\Users\feihu\OneDrive - McGill University\matlab\Analysis based on larvae\Fall_2021_Replication\All\Choreography analysis";

    %% load data
    %get raw choreography data (x,y,speed) from the .txt files
    [dat,uname]=get_data_from_dat_file(name_of_folder,ii);
    %concatenate all the data, which will be used for plotting trajectory
    %later
    [dat,xmax,ymax]=vertcat_all_data (dat);
    %grouped data based on their aninmal number (AN)
    dat_grouped=group_all_data_based_on_AN(dat);
    %add two new variables xcentered and ycentered which can be used for
    %the centered trajectory later (all trajectory will have a start point at [0,0])
    [dat,dat_grouped]=generate_centered_data(dat,dat_grouped);

    %% remove both the first 120s and the edge from the data
    if ~isempty(both_edge) && ~isempty(both_ns)
        [dat_edge,dat_grouped_edge]=eliminate_edge(dat,dat_grouped,both_edge);
        [dat_both,dat_grouped_both]=remove_the_first_ns(dat_edge,dat_grouped_edge,both_ns);
        clear dat_edge dat_grouped_edge
    end

    %% eliminate the first 120s of the data
    if ~isempty(ns)
        [dat_ns,dat_grouped_ns]=remove_the_first_ns(dat,dat_grouped,ns);
    end
    %% remove the datapoints at the edge
    if ~isempty(edge)
        [dat_edge,dat_grouped_edge]=eliminate_edge(dat,dat_grouped,edge);
    end
    %% plot
    % 1) the tracjectory
    color=[0.0353,0.2784,0.5020;0.051,0.4157,0.749; 0.4392,0.7333,1.0;0.8392,0.2706,0.3137;0,0,0];
    name={'10^-^1GA','10^-^2GA','10^-^3GA','10^-^5EA','H2O'};
    row=2;col=3;num=1;str='';
    if plot_input(1)=='y'
        ii=plot_trajectory(dat,color,name,ii,ii_num,row,col,num,str,outdir);
    end
    % 2)plot the centered trajectory (all the trajectories start from [0,0])
    if plot_input(2)=='y'
        ii=plot_centered_trajectory(dat,color,name,ii,ii_num,row,col,num,str,outdir);
    end
    % 3)plot heatmap (another representation of trajectory)
    if plot_input(3)=='y'
        if ii==1
            lim=[];
        end
        lim=plot_heatmap(dat,name,ii,ii_num,row,col,xmax,ymax,num,str,lim,outdir);
    end
    % 4)calculate and plot the sum of tracking time at each x bins
    if plot_input(4)=='y'
        if ii==1
            tt=[];
        end
        tt=calculate_tracking_time(tt,dat_grouped,xmax,ii,ii_num,name,num,row,col,color,str,outdir);
    end
    %% Calculate the navigational index
    
    %1)calculate the navigational index:v_x/speed for each timepoint
    if plot_input(5)=='y'
    %eliminate the larvae that has tracking time less than 10s
    dat_grouped=eliminate_tt(dat_grouped,10);
    if ii==1
        PI_mean=[];PI_t_mean=[];PI_t_sem=[];PI_sem=[];PI_timeseries={};
    end
    t=300;
    [dat_grouped,PI_mean,PI_sem,PI_t_mean,PI_t_sem]=calculate_navigational_index_velocity_vs_speed(dat_grouped,ii,t,PI_mean,PI_sem,PI_t_mean,PI_t_sem);
    %2) Caculate the timeseries for navigational index ((number of larvae close to the odor-number of larvae away from the odor)/total number of larvare)
    time=0:10:900;
    [dat_grouped,PI_timeseries]=calculate_preferential_index_vs_time(dat_grouped,ii,time,xmax,PI_timeseries);
    %% plot the navigational index
    if ii==ii_num
        PI_mean=plot_PI(PI_mean,PI_sem,PI_t_mean,PI_t_sem,num,str,color,name,ii_num,t,outdir);
        num=plot_PI_timeseries(PI_timeseries,name,color,time,num,outdir,str);
    end
    end

    %% repeat the plot for different data-eliminate the first 120s, the edge and both
    %plots and data for eliminating the first 120s
    if exist('dat_ns','var')
        str_ns=append(',eliminate the first ',num2str(ns),'s');num_ns=10;
        if plot_input(1)=='y'
            ii=plot_trajectory(dat_ns,color,name,ii,ii_num,row,col,num_ns,str_ns,outdir);
        end

        if plot_input(2)=='y'
            ii=plot_centered_trajectory(dat_ns,color,name,ii,ii_num,row,col,num_ns,str_ns,outdir);
        end

        if plot_input(3)=='y'
            if ii==1
                lim_ns=[];
            end
            lim_ns=plot_heatmap(dat_ns,name,ii,ii_num,row,col,xmax,ymax,num_ns,str_ns,lim_ns,outdir);
        end

        if plot_input(4)=='y'
            if ii==1
                tt_ns=[];
            end
            tt_ns=calculate_tracking_time(tt_ns,dat_grouped_ns,xmax,ii,ii_num,name,num_ns,row,col,color,str_ns,outdir);

        end

        if  plot_input(5)=='y'
            dat_grouped_ns=eliminate_tt(dat_grouped_ns,10);
            if ii==1
                PI_mean_ns=[];PI_t_mean_ns=[];PI_t_sem_ns=[];PI_sem_ns=[];PI_timeseries_ns={};
            end
            t=300;
            [dat_grouped_ns,PI_mean_ns,PI_sem_ns,PI_t_mean_ns,PI_t_sem_ns]=calculate_navigational_index_velocity_vs_speed(dat_grouped_ns,ii,t,PI_mean_ns,PI_sem_ns,PI_t_mean_ns,PI_t_sem_ns);
            time=0:10:900;
            [dat_grouped_ns,PI_timeseries_ns]=calculate_preferential_index_vs_time(dat_grouped_ns,ii,time,xmax,PI_timeseries_ns);
        
        if ii==ii_num
            PI_mean_ns=plot_PI(PI_mean_ns,PI_sem_ns,PI_t_mean_ns,PI_t_sem_ns,num_ns,str_ns,color,name,ii_num,t,outdir);
            num_ns=plot_PI_timeseries(PI_timeseries_ns,name,color,time,num_ns,outdir,str_ns);
        end
        end
    end

%% plots and data for the data with edge eliminated
if exist ('dat_edge','var')
    str_edge=',eliminate the edge';num_edge=20;
    if plot_input(1)=='y'
        ii=plot_trajectory(dat_edge,color,name,ii,ii_num,row,col,num_edge,str_edge,outdir);
    end

    if plot_input(2)=='y'
        ii=plot_centered_trajectory(dat_edge,color,name,ii,ii_num,row,col,num_edge,str_edge,outdir);
    end

    if plot_input(3)=='y'
        if ii==1
            lim_edge=[];
        end
        lim_edge=plot_heatmap(dat_edge,name,ii,ii_num,row,col,xmax,ymax,num_edge,str_edge,lim_edge,outdir);
    end

    if plot_input(4)=='y'
        if ii==1
            tt_edge=[];
        end
        tt_edge=calculate_tracking_time(tt_edge,dat_grouped_edge,xmax,ii,ii_num,name,num_edge,row,col,color,str_edge,outdir);
    end

    if plot_input(5)=='y'
        dat_grouped_edge=eliminate_tt(dat_grouped_edge,10);
        if ii==1
            PI_mean_edge=[];PI_t_mean_edge=[];PI_t_sem_edge=[];PI_sem_edge=[];PI_timeseries_edge={};
        end
        t=300;
        [dat_grouped_edge,PI_mean_edge,PI_sem_edge,PI_t_mean_edge,PI_t_sem_edge]=calculate_navigational_index_velocity_vs_speed(dat_grouped_edge,ii,t,PI_mean_edge,PI_sem_edge,PI_t_mean_edge,PI_t_sem_edge);
        time=0:10:900;
        [dat_grouped_edge,PI_timeseries_edge]=calculate_preferential_index_vs_time(dat_grouped_edge,ii,time,xmax,PI_timeseries_edge);

        if ii==ii_num
            PI_mean_edge=plot_PI(PI_mean_edge,PI_sem_edge,PI_t_mean_edge,PI_t_sem_edge,num_edge,str_edge,color,name,ii_num,t,outdir);
            num_edge=plot_PI_timeseries(PI_timeseries_edge,name,color,time,num_edge,outdir,str_edge);
        end
    end
end

%% eliminate both
if exist('dat_both','var')
    str_both=append(',eliminate both the edge and the first ',num2str(ns),'s');num_both=30;
    if plot_input(1)=='y'
        ii=plot_trajectory(dat_both,color,name,ii,ii_num,row,col,num_both,str_both,outdir);
    end

    if plot_input(2)=='y'
        ii=plot_centered_trajectory(dat_both,color,name,ii,ii_num,row,col,num_both,str_both,outdir);
    end

    if plot_input(3)=='y'
        if ii==1
            lim_both=[];
        end
        lim_both=plot_heatmap(dat_both,name,ii,ii_num,row,col,xmax,ymax,num_both,str_both,lim_both,outdir);
    end

    if plot_input(4)=='y'
        if ii==1
            tt_both=[];
        end
        tt_both=calculate_tracking_time(tt_both,dat_grouped_both,xmax,ii,ii_num,name,num_both,row,col,color,str_both,outdir);
    end

    if plot_input(5)=='y'
        dat_grouped_both=eliminate_tt(dat_grouped_both,10);
        if ii==1
            PI_mean_both=[];PI_t_mean_both=[];PI_t_sem_both=[];PI_sem_both=[];PI_timeseries_both={};
        end
        t=300;
        [dat_grouped_both,PI_mean_both,PI_sem_both,PI_t_mean_both,PI_t_sem_both]=calculate_navigational_index_velocity_vs_speed(dat_grouped_both,ii,t,PI_mean_both,PI_sem_both,PI_t_mean_both,PI_t_sem_both);
        time=0:10:900;
        [dat_grouped_both,PI_timeseries_both]=calculate_preferential_index_vs_time(dat_grouped_both,ii,time,xmax,PI_timeseries_both);


        if ii==ii_num
            PI_mean_both=plot_PI(PI_mean_both,PI_sem_both,PI_t_mean_both,PI_t_sem_both,num_both,str_both,color,name,ii_num,t,outdir);
            num_both=plot_PI_timeseries(PI_timeseries_both,name,color,time,num_both,outdir,str_both);
        end
    end
end

%% save the data file
file_dir=fullfile(outdir,'data');
if ~isdir(file_dir)
    mkdir(file_dir);
end
name_file={'10E-1GA','10E-2GA','10E-3GA','10E-5EA','H2O'};
filename=fullfile(file_dir,append(name_file(ii),'_data.mat'));
if isfile(filename)
    delete(filename);
end
    save( filename, 'dat*')


if ii==5
    %save the preferential index data
    filename=fullfile(file_dir,'_PI_data.mat');

    if isfile(filename)
        delete(filename);
    end
        save( filename, 'PI*')
   
end
clearvars -except PI* ii_num lim* ns edge both_edge both_ns plot_input tt*
end
