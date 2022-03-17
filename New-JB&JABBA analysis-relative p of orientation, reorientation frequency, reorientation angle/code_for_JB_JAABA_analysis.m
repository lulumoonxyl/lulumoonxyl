
ii_num=5;%this depends on the conditions/subfolders you have
[ns,edge,both_ns,both_edge,plot_input]=get_input_from_user(ii_num);
for ii=1:ii_num
    %set output folder
    outdir="C:\Users\tomokolab\OneDrive - McGill University\matlab\Analysis based on larvae\Fall_2021_Replication\All\JB_JAABA_analysis";
    %% load JB data
    name='E:\organized data\jb_result';
    [dat_JB,uname_JB]=get_JB_data(name,ii);
    %get x and y
    dat_JB=get_xy_center(dat_JB);
    %% load JAABA data
    name_JAABA='E:\organized data\jaaba_result';
    [dat_JAABA,uname_JAABA]=get_JAABA_data(name_JAABA,ii);
    %% rearrange the data array for the jaaba
    dat_JAABA=arrange_JAABA_data(dat_JAABA);
    dat_JAABA=align_JB_JAABA_data(dat_JB,dat_JAABA);
    %% remove both the first 120s and the edge from the JB and JAABA data array
    if ~isempty(both_ns)&&~isempty(both_edge)
        [dat_JB_ns,dat_JAABA_ns]=remove_the_first_ns_JAABA_and_JB(dat_JAABA,both_ns,dat_JB);
        [dat_JB_both,dat_JAABA_both]=eliminate_edge_JB_and_JAABA(dat_JB_ns,both_edge,dat_JAABA_ns);
        if isempty(ns)
            clear dat_JB_ns dat_JAABA_ns
        elseif ns~=both_ns
            clear dat_JB_ns dat_JAABA_ns
            %remove only the first n seconds
            [dat_JB_ns,dat_JAABA_ns]=remove_the_first_ns_JAABA_and_JB(dat_JAABA,ns,dat_JB);
        end
    end
    %remove only the edge
    if ~isempty(edge)
        [dat_JB_edge,dat_JAABA_edge]=eliminate_edge_JB_and_JAABA(dat_JB,edge,dat_JAABA);
    end
    
    %% calculate the relative probability of orientation
    % calculate the total trakcing time for each animals and remove the cells that
    %have tracking time less than 5s
    [dat_JB,dat_JAABA]=calculate_tsum(dat_JB,dat_JAABA,5);
    vector=[-1 0]; %the position where odor is at
    dat_JB=calculate_heading_direction(dat_JB,vector);
    %plot a single trajectory to confirm the orientation is calculated
    %correctly
    %         figure(1); i=2;
    %         idx0=dat_JAABA.t0_idx{i,1};idx1=dat_JAABA.t1_idx{i,1};
    %         idx3=[];
    %         for j=1:length(idx1)
    %             idx3=vertcat(idx3,(idx0(j):idx1(j))')
    %         end
    %         hold on
    %         plot(0-dat_JB.x{i,1},dat_JB.y{i,1},'g.');
    %         plot(0-dat_JB.x{i,1}(idx3),dat_JB.y{i,1}(idx3),'k.')
    %         plot(0-dat_JB.x{i,1}(idx0),dat_JB.y{i,1}(idx0),'r.');
    %         plot(0-dat_JB.x{i,1}(idx1),dat_JB.y{i,1}(idx1),'y.')
    %         text(0-dat_JB.x{i,1}(1),dat_JB.y{i,1}(1),'start');
    %         hold off
    bins=20;
    name={'10^-^3GA','10^-^2GA','10^-^1GA','10^-^5EA','H2O'};row=2;col=3;str='';num=1;
    
    
    %% plot the relative probability of orientation
    if isequal(plot_input(1),'y')
        dat_JB=calculate_relative_p(dat_JB,bins);
        ii=plot_relative_p(dat_JB,ii,num,name,ii_num,row,col,str,bins,outdir);
    end
    %% combine some turning events
    dat_JAABA=combine_and_cancel_JAABA_turn(dat_JAABA);
    %get reorientation angle, x position where turning occurs,
    %pre/post-turning angles
    dis=5;%dis is used for pre/post-turning angles (we use the point 'dis' before the turning event to calculate the vector)
    dat_JB=change_in_theta_for_turning_event(dat_JB,dat_JAABA,dis,vector);
    %% calculate the turning frequency and reorientation angle
    if isequal(plot_input(2),'y')
        split_or_not='y';xbins=20;deg=60;large_or_small='large';color=[0.0353,0.2784,0.5020;0.051,0.4157,0.749; 0.4392,0.7333,1.0;0.8392,0.2706,0.3137;0,0,0];
        [result,result_split]=get_turning_frequency_and_theta(dat_JB,deg,split_or_not,xbins,large_or_small);
        if ii==1
            data={};
        end
        data=calculate_mean_tf_theta(result,data,ii);
        %plot the result
        if ii==ii_num
            num=num+1;
            ii_num=plot_tf_theta(result.x,data,num,ii_num,name,color,str,deg,xbins,large_or_small);
        end
        if exist('result_split','var')
            if ii==1
                data_split={};
            end
            data_split=calculate_mean_tf_theta(result_split,data_split,ii);
            if ii==ii_num
                str_q={'q1','q2','q3','q4'};
                for i=1:length(str_q)
                num=num+1;
                data_split_q=get_data_for_each_quadrant(data_split,str_q(i));
                str_plot=append(str_q(i),str);
                ii_num=plot_tf_theta(result.x,data_split_q,num,ii_num,name,color,str_plot,deg,xbins,large_or_small);
                end 
            end
        end
    end


    %% repeat for data with edge removed
    if exist ('dat_JB_edge','var')
        [dat_JB_edge,dat_JAABA_edge]=calculate_tsum(dat_JB_edge,dat_JAABA_edge,5);
        dat_JB_edge=calculate_heading_direction(dat_JB_edge,vector);
        
        bins_edge=20;str_edge=',eliminate the edge';num_edge=10;
        
        
        if isequal(plot_input(1),'y')
            dat_JB_edge=calculate_relative_p(dat_JB_edge,bins_edge);
            ii=plot_relative_p(dat_JB_edge,ii,num_edge,name,ii_num,row,col,str_edge,bins_edge,outdir);
        end
        
        dat_JAABA_edge=combine_and_cancel_JAABA_turn(dat_JAABA_edge);
        dat_JB_edge=change_in_theta_for_turning_event(dat_JB_edge,dat_JAABA_edge,dis,vector);
        
        if isequal(plot_input(2),'y')
            [result_edge,result_split_edge]=get_turning_frequency_and_theta(dat_JB_edge,deg,split_or_not,xbins,large_or_small);
            if ii==1
                data_edge={};
            end
            data_edge=calculate_mean_tf_theta(result_edge,data_edge,ii);
            
            if exist('result_split_edge','var')
                if ii==1
                    data_split_edge={};
                end
                data_split_edge=calculate_mean_tf_theta(result_split_edge,data_split_edge,ii);
                
            end
        end
    end
    %% repeat for data with first n seconds removed
    if exist ('dat_JB_ns','var')
        [dat_JB_ns,dat_JAABA_ns]=calculate_tsum(dat_JB_ns,dat_JAABA_ns,5);
        dat_JB_ns=calculate_heading_direction(dat_JB_ns,vector);
        
        bins_ns=20;str_ns=append(',eliminate the first ',num2str(ns),'s');num_ns=20;
        
        
        if isequal(plot_input(1),'y')
            dat_JB_ns=calculate_relative_p(dat_JB_ns,bins_ns);
            ii=plot_relative_p(dat_JB_ns,ii,num_ns,name,ii_num,row,col,str_ns,bins_ns,outdir);
        end
        
        dat_JAABA_ns=combine_and_cancel_JAABA_turn(dat_JAABA_ns);
        dat_JB_ns=change_in_theta_for_turning_event(dat_JB_ns,dat_JAABA_ns,dis,vector);
        if isequal(plot_input(2),'y')
            
            [result_ns,result_split_ns]=get_turning_frequency_and_theta(dat_JB_ns,deg,split_or_not,xbins,large_or_small);
            if ii==1
                data_ns={};
            end
            data_ns=calculate_mean_tf_theta(result_ns,data_ns,ii);
            
            if exist('result_split_ns','var')
                if ii==1
                    data_split_ns={};
                end
                data_split_ns=calculate_mean_tf_theta(result_split_ns,data_split_ns,ii);
                
            end
        end
    end
    %% repeat for data with both the edge and the first n seconds removed
    if exist ('dat_JB_both','var')
        [dat_JB_both,dat_JAABA_both]=calculate_tsum(dat_JB_both,dat_JAABA_both,5);
        dat_JB_both=calculate_heading_direction(dat_JB_both,vector);
        
        bins_both=20;str_both=append(',eliminate the edge and the first ',num2str(both_ns),'s');num_both=30;
        
        
        if isequal(plot_input(1),'y')
            dat_JB_both=calculate_relative_p(dat_JB_both,bins_both);
            ii=plot_relative_p(dat_JB_both,ii,num_both,name,ii_num,row,col,str_both,bins_both,outdir);
        end
        
        dat_JAABA_both=combine_and_cancel_JAABA_turn(dat_JAABA_both);
        dat_JB_both=change_in_theta_for_turning_event(dat_JB_both,dat_JAABA_both,dis,vector);
        
        if isequal(plot_input(2),'y')
            
            [result_both,result_split_both]=get_turning_frequency_and_theta(dat_JB_both,deg,split_or_not,xbins,large_or_small);
            if ii==1
                data_both={};
            end
            data_both=calculate_mean_tf_theta(result_both,data_both,ii);
            
            if exist('result_split_both','var')
                if ii==1
                    data_split_both={};
                end
                data_split_both=calculate_mean_tf_theta(result_split_both,data_split_both,ii);
                
            end
        end
    end
end