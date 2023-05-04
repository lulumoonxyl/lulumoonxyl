function []=turn_freq_4_quadrants(tracker_num,genotype,condition,filename,output_name,varargin)
%% get the directory for all the data
w=1;

for i=1:length(genotype)
    for j=1:length(condition)
        cond(w,1)=fullfile("/project/6010970/screen/olfactory_output/JB_JAABA",tracker_num,genotype{i,1},condition{j,1},filename);
        name(w,1)=append(genotype{i,1},"@",condition{j,1});
        w=w+1;
    end
end

%% set the output name
outdir=fullfile("/project/6010970/screen/olfactory_output/Properties_turn_events",output_name);
ii_num=length(cond);
%% set the properties for plots
row=3;col=3;
color=[0.75,0.75,0.75;0.627,0.82,1;0.659,0.659,0.659;0.871,0.416,0.451];
hbin=45;% for binning data
input_cond=[];%it will have a same length as name, telling us whether the input is 'ctr' or 'exp'
row_comp=2; col_comp=2; % these two will be used for the plot to show the comparison of ctr vs exp group, this number should be based on how many comparison you have
tmax=900;% tbin is used for binning the turn event's t0s
xmax=225; deg_list=60;
l=3;% how many quadrant, either 3 or 4

for i=1:2:length(varargin)
    if strcmp ('row',varargin{i})
        row=varargin{i+1};
    elseif strcmp ('col',varargin{i})
        col=varargin{i+1};
    elseif strcmp ('color',varargin{i})
        color=varargin{i+1};
    elseif strcmp(varargin{i},'name')
        clear name;
        name=string(varargin{i+1})
    elseif strcmp(varargin{i},'input_cond')
        input_cond=string(varargin{i+1});
    elseif strcmp(varargin{i},'row_comp')
        row_comp=varargin{i+1};
    elseif strcmp(varargin{i},'col_comp')
        col_comp=varargin{i+1};
    elseif strcmp(varargin{i},'xmax')
        xmax=varargin{i+1};
    elseif strcmp(varargin{i},'deg_l')
        deg_list=varargin{i+1};
    elseif strcmp(varargin{i},'tmax')
        tmax=varargin{i+1};
    elseif strcmp(varargin{i},'quadrant')
        l=varargin{i+1};
    end
end

if l==4
    q_legends=["Toward odor"; "Neutral upper"; "Away from odor";"Neutral lower"];% the strings used for the quadrants
elseif l==3
    q_legends=["Toward odor"; "Neutral"; "Away from odor"];
else
    warning('Please put the correct number for quadrant:3 or 4')
end
for ii=1:ii_num
    load(cond(ii));
    %% 1-1) Compute the turning frequency-->count of turning event/time
    % get the data ready for binning
    pre=vertcat(dat.pre_deg{:});
    %     post=vertcat(dat.post_deg{:});
    turn_x=vertcat(dat.turn_x{:});
    %     turn_y=vertcat(dat.turn_y{:});
    reori_deg=vertcat(dat.reorient_deg_abs{:});
    t0s=vertcat(dat.t0s{:});
    if isa(t0s,'cell')
        t0s=vertcat(t0s{:});
    end 
    if length(deg_list)==1
        deg_l=deg_list;%later may be use idx for the large turning event
    else
        deg_l=deg_list(ii);
    end
    disp(append("The cutoff for large turning event of the group ", name(ii)," is ",num2str(deg_l)," deg"));

    idx=find(reori_deg>=deg_l);%later may be use idx for the large turning event
    xbin2=xmax/2;
    hdseries=[-180:hbin:180]';
    %     hd_series1=[-180+hbin/2:hbin:180-hbin/2]';
    %% 1-2) bin the data based on the heading direction before turning --:UNIT WILL BE TURNING PER MIN

    if ii==1
        tf=zeros(ii_num,l); tf_prox=zeros(ii_num,l); tf_dist=zeros(ii_num,l);
        tf_L=zeros(ii_num,l); tf_prox_L=zeros(ii_num,l); tf_dist_L=zeros(ii_num,l);
    end

    t_x1=bin_data_sum_t(dat.orientation,dat.et,hdseries);
    t_x=sum(t_x1,1,'omitnan')';

    %compute the turning freq by dividing the heading directions into 4
    %quadrants: 1) -45-->45 2) 45-->135 3) 135-->180 &&
    %-180-->-135 4)-135-->-45 or combine quadrants 2) and 4) together

    %simply combine two 45 deg quadrants into larger one
    if l==4
        t_q3=assign_num_4quadrants(t_x);
    elseif l==3
        t_q3=assign_num_3quadrants(t_x);
    end
    %-->modify this function based on the quadrants combination you want

    c_x=bin_data_count(pre,hdseries);

    if l==4
        c_q3=assign_num_4quadrants(c_x);
    elseif l==3
        c_q3=assign_num_3quadrants(c_x);
    end

    tf(ii,:)=transpose((c_q3./t_q3).*60);
    %% 1-3) Get the turning frequency for turning events larger than 60 deg
    c_x_L=bin_data_count(pre(idx),hdseries);

    if l==4
        c_q3_L=assign_num_4quadrants(c_x_L);
    elseif l==3
        c_q3_L=assign_num_3quadrants(c_x_L);
    end

    tf_L(ii,:)=transpose((c_q3_L./t_q3).*60);
    clear t_q3 t_x1 c_x* c_q3* t_x

    %% 2-1) bin the data based on x pos of each turning event (proximate or distant)
    t_x=bin_data_sum_t2(dat.orientation,dat.x,hdseries,[0:xbin2:xmax]',dat.et);

    if l==4
        t_q3=assign_num_4quadrants(t_x);
    elseif l==3
        t_q3=assign_num_3quadrants(t_x);
    end

    c_x=bin_data_count2(pre,hdseries,turn_x,[0:xbin2:xmax]');

    if l==4
        c_q3=assign_num_4quadrants(c_x);
    elseif l==3
        c_q3=assign_num_3quadrants(c_x);
    end

    tf_x=(c_q3./t_q3).*60;

    tf_prox(ii,:)=transpose(tf_x(:,2));
    tf_dist(ii,:)=transpose(tf_x(:,1));

    %% 2-2) Get the turning frequency for turning events larger than 60 deg based on the x position (prox and distant)
    c_x2_L=bin_data_count2(pre(idx),hdseries,turn_x(idx),[0:xbin2:xmax]');

    if l==4
        c_q3_L=assign_num_4quadrants(c_x2_L);
    elseif l==3
        c_q3_L=assign_num_3quadrants(c_x2_L);
    end
    tf_x2_L=(c_q3_L./t_q3).*60;

    tf_prox_L(ii,:)=transpose(tf_x2_L(:,2));
    tf_dist_L(ii,:)=transpose(tf_x2_L(:,1));
    clear c_x2_L t_q3 tf_x c_q3* t_x c_x tf_x2_L

    %% 3-1) bin the data based on t0s of each turning event
    if ii==1
        tf_t1=zeros(ii_num,l); tf_t2=zeros(ii_num,l); tf_t3=zeros(ii_num,l);
        tf_t1_L=zeros(ii_num,l); tf_t2_L=zeros(ii_num,l); tf_t3_L=zeros(ii_num,l);
    end

    t=[0:tmax/3:tmax]';% split it in three time sections
    t_t=bin_data_sum_t2(dat.orientation,dat.et,hdseries,t,dat.et);

    if l==4
        t_q3=assign_num_4quadrants(t_t);
    elseif l==3
        t_q3=assign_num_3quadrants(t_t);
    end

    c_t=bin_data_count2(pre,hdseries,t0s,t);

    if l==4
        c_q3=assign_num_4quadrants(c_t);
    elseif l==3
        c_q3=assign_num_3quadrants(c_t);
    end

    tf_t=(c_q3./t_q3).*60;
    tf_t1(ii,:)=transpose(tf_t(:,1)); tf_t2(ii,:)=transpose(tf_t(:,2)); tf_t3(ii,:)=transpose(tf_t(:,3));

    c_t_L=bin_data_count2(pre(idx),hdseries,t0s(idx),t);

    if l==4
        c_q3_L=assign_num_4quadrants(c_t_L);
    elseif l==3
        c_q3_L=assign_num_3quadrants(c_t_L);
    end
    tf_t_L=(c_q3_L./t_q3).*60;
    tf_t1_L(ii,:)=transpose(tf_t_L(:,1)); tf_t2_L(ii,:)=transpose(tf_t_L(:,2)); tf_t3_L(ii,:)=transpose(tf_t_L(:,3));

    clear tf_t c_q3 t_t c_t c_q3_L c_t_L t_q3 tf_t_L
    %% 1) & 2) & 3) use bootstrapping to get 95% CI

    len=200;
    tf_bt=zeros(l,len); tf_bt_prox=zeros(l,len); tf_bt_dist=zeros(l,len);
    tf_bt_L=zeros(l,len);tf_bt_prox_L=zeros(l,len);  tf_bt_dist_L=zeros(l,len);

    tf_bt_t1=zeros(l,len); tf_bt_t2=zeros(l,len); tf_bt_t3=zeros(l,len);
    tf_bt_t1_L=zeros(l,len); tf_bt_t2_L=zeros(l,len); tf_bt_t3_L=zeros(l,len);

    rng('default');
    s = rng;
    [~,bt]=bootstrp(len,[],[1:length(dat.AN)]);
    %get 200 bootstraps samples and compute the turning freq for each of
    %them

    for k=1:len
        data=dat(bt(:,k),:);
        pre_bt=vertcat(data.pre_deg{:});
        reori_deg_bt=vertcat(data.reorient_deg_abs{:});
        turn_x_bt=vertcat(data.turn_x{:});

%         if k==77
%             index=bt(:,k)
%             outdir1=fullfile(outdir,"turn freq","rms turn size",append(num2str(l)," quadrants"));
%             if ~isfolder(outdir1)
%                 mkdir(outdir1);
%             end
% 
%             save(append(outdir1,'index.mat'),'index');
%         end
        t0s_bt=vertcat(data.t0s{:});
        if isa(t0s_bt,'cell')
            t0s_bt=vertcat(t0s_bt{:});
        end
        idx_bt=find(reori_deg_bt>=deg_l);

        % bin the turning events based on their pre-turning heading
        % directions
        t_bt1=bin_data_sum_t(data.orientation,data.et,hdseries);
        t_bt=sum(t_bt1,1,'omitnan')';

        if l==4
            t_bt_q3=assign_num_4quadrants(t_bt);
        elseif l==3
            t_bt_q3=assign_num_3quadrants(t_bt);
        end

        c_bt=bin_data_count(pre_bt,hdseries);


        if l==4
            c_bt_q3=assign_num_4quadrants(c_bt);
        elseif l==3
            c_bt_q3=assign_num_3quadrants(c_bt);
        end
        tf_bt(:,k)=(c_bt_q3./t_bt_q3).*60;

        c_bt_L=bin_data_count(pre_bt(idx_bt),hdseries);


        if l==4
            c_bt_q3_L=assign_num_4quadrants(c_bt_L);
        elseif l==3
            c_bt_q3_L=assign_num_3quadrants(c_bt_L);
        end
        tf_bt_L(:,k)=(c_bt_q3_L./t_bt_q3).*60;
        clear c_bt_q3* t_bt_q3

        % bin the turning events based on their positions and heading
        % directions
        t_bt2=bin_data_sum_t2(data.orientation,data.x,hdseries,[0:xbin2:xmax]',data.et);

        if l==4
            t_bt2_q3=assign_num_4quadrants(t_bt2);
        elseif l==3
            t_bt2_q3=assign_num_3quadrants(t_bt2);
        end

        c_bt2=bin_data_count2(pre_bt,hdseries,turn_x_bt,[0:xbin2:xmax]');

        if l==4
            c_bt2_q3=assign_num_4quadrants(c_bt2);
        elseif l==3
            c_bt2_q3=assign_num_3quadrants(c_bt2);
        end

        c_bt2_L=bin_data_count2(pre_bt(idx_bt),hdseries,turn_x_bt(idx_bt),[0:xbin2:xmax]');


        if l==4
            c_bt2_q3_L=assign_num_4quadrants(c_bt2_L);
        elseif l==3
            c_bt2_q3_L=assign_num_3quadrants(c_bt2_L);
        end

        tf_bt2=(c_bt2_q3./t_bt2_q3).*60;
        tf_bt_prox(:,k)=tf_bt2(:,2);tf_bt_dist(:,k)=tf_bt2(:,1);

        tf_bt2_L=(c_bt2_q3_L./t_bt2_q3).*60;
        tf_bt_prox_L(:,k)=tf_bt2_L(:,2);tf_bt_dist_L(:,k)=tf_bt2_L(:,1);
        clear t_bt2* c_bt2* tf_bt2*

        %for binning the turnning events based on the timing
        t_bt_t=bin_data_sum_t2(dat.orientation,dat.et,hdseries,t,dat.et);
     
        if l==4
            t_bt_t_q3=assign_num_4quadrants(t_bt_t);
        elseif l==3
            t_bt_t_q3=assign_num_3quadrants(t_bt_t);
        end

        c_bt_t=bin_data_count2(pre_bt,hdseries,t0s_bt,t);
        if l==4
            c_bt_t_q3=assign_num_4quadrants(c_bt_t);
        elseif l==3
            c_bt_t_q3=assign_num_3quadrants(c_bt_t);
        end

        tf_bt_t=(c_bt_t_q3./t_bt_t_q3).*60;

        c_bt_t_L=bin_data_count2(pre_bt(idx_bt),hdseries,t0s_bt(idx_bt),t);
   
         if l==4
            c_bt_t_q3_L=assign_num_4quadrants(c_bt_t_L);
        elseif l==3
            c_bt_t_q3_L=assign_num_3quadrants(c_bt_t_L);
        end

        tf_bt_t_L=(c_bt_t_q3_L./t_bt_t_q3).*60;

        tf_bt_t1(:,k)=tf_bt_t(:,1); tf_bt_t2(:,k)=tf_bt_t(:,2); tf_bt_t3(:,k)=tf_bt_t(:,3);
        tf_bt_t1_L(:,k)=tf_bt_t_L(:,1); tf_bt_t2_L(:,k)=tf_bt_t_L(:,2); tf_bt_t3_L(:,k)=tf_bt_t_L(:,3);
        clear tf_bt_t_L c_bt_t*  tf_bt_t t_bt_t*

    end

    % Save data for computing the bootstrap p value
    if ii==1
        tf_bt_all=cell(ii_num,1); tf_bt_prox_all=cell(ii_num,1);tf_bt_dist_all=cell(ii_num,1);
        tf_bt_all_L=cell(ii_num,1); tf_bt_prox_all_L=cell(ii_num,1);tf_bt_dist_all_L=cell(ii_num,1);

        tf_bt_t1_all=cell(ii_num,1); tf_bt_t2_all=cell(ii_num,1);tf_bt_t3_all=cell(ii_num,1);
        tf_bt_t1_all_L=cell(ii_num,1); tf_bt_t2_all_L=cell(ii_num,1);tf_bt_t3_all_L=cell(ii_num,1);
    end

    tf_bt_all{ii,1}=tf_bt; tf_bt_prox_all{ii,1}=tf_bt_prox; tf_bt_dist_all{ii,1}=tf_bt_dist;
    tf_bt_all_L{ii,1}=tf_bt_L; tf_bt_prox_all_L{ii,1}=tf_bt_prox_L; tf_bt_dist_all_L{ii,1}=tf_bt_dist_L;

    tf_bt_t1_all{ii,1}=tf_bt_t1; tf_bt_t2_all{ii,1}=tf_bt_t2; tf_bt_t3_all{ii,1}=tf_bt_t3;
    tf_bt_t1_all_L{ii,1}=tf_bt_t1_L; tf_bt_t2_all_L{ii,1}=tf_bt_t2_L; tf_bt_t3_all_L{ii,1}=tf_bt_t3_L;

    %% Use the data to compute the 95% CI
    % for distant and proximate regions
    tf_diff=tf_bt-tf(ii,:)'; tf_prox_diff=tf_bt_prox-tf_prox(ii,:)'; tf_dist_diff=tf_bt_dist-tf_dist(ii,:)';
    tf_diff_L=tf_bt_L-tf_L(ii,:)'; tf_prox_diff_L=tf_bt_prox_L-tf_prox_L(ii,:)'; tf_dist_diff_L=tf_bt_dist_L-tf_dist_L(ii,:)';

    tf_diff=sort(tf_diff,2); tf_prox_diff=sort(tf_prox_diff,2);tf_dist_diff=sort(tf_dist_diff,2);
    tf_diff_L=sort(tf_diff_L,2); tf_prox_diff_L=sort(tf_prox_diff_L,2);tf_dist_diff_L=sort(tf_dist_diff_L,2);
    % initialize empty arrays for sem from the bootstrapping results
    if ii==1
        tf_pos=zeros(ii_num,l); tf_prox_pos=zeros(ii_num,l); tf_dist_pos=zeros(ii_num,l);
        tf_pos_L=zeros(ii_num,l); tf_prox_pos_L=zeros(ii_num,l); tf_dist_pos_L=zeros(ii_num,l);
        tf_neg=zeros(ii_num,l); tf_prox_neg=zeros(ii_num,l); tf_dist_neg=zeros(ii_num,l);
        tf_neg_L=zeros(ii_num,l); tf_prox_neg_L=zeros(ii_num,l); tf_dist_neg_L=zeros(ii_num,l);

    end
    tf_neg(ii,:)=tf_diff(:,5)'; tf_pos(ii,:)=tf_diff(:,195)';
    tf_neg_L(ii,:)=tf_diff_L(:,5)'; tf_pos_L(ii,:)=tf_diff_L(:,195)';

    tf_dist_neg(ii,:)=tf_dist_diff(:,5)'; tf_dist_pos(ii,:)=tf_dist_diff(:,195)';
    tf_dist_neg_L(ii,:)=tf_dist_diff_L(:,5)'; tf_dist_pos_L(ii,:)=tf_dist_diff_L(:,195)';

    tf_prox_neg(ii,:)=tf_prox_diff(:,5)'; tf_prox_pos(ii,:)=tf_prox_diff(:,195)';
    tf_prox_neg_L(ii,:)=tf_prox_diff_L(:,5)'; tf_prox_pos_L(ii,:)=tf_prox_diff_L(:,195)';
    %% for different timing

    tf_diff_t1=tf_bt_t1-tf_t1(ii,:)'; tf_diff_t2=tf_bt_t2-tf_t2(ii,:)'; tf_diff_t3=tf_bt_t3-tf_t3(ii,:)';
    tf_diff_t1_L=tf_bt_t1_L-tf_t1_L(ii,:)'; tf_diff_t2_L=tf_bt_t2_L-tf_t2_L(ii,:)'; tf_diff_t3_L=tf_bt_t3_L-tf_t3_L(ii,:)';

    tf_diff_t1=sort(tf_diff_t1,2); tf_diff_t2=sort(tf_diff_t2,2);tf_diff_t3=sort(tf_diff_t3,2);
    tf_diff_t1_L=sort(tf_diff_t1_L,2); tf_diff_t2_L=sort(tf_diff_t2_L,2);tf_diff_t3_L=sort(tf_diff_t3_L,2);

    if ii==1
        tf_t1_neg=zeros(ii_num,l); tf_t2_neg=zeros(ii_num,l); tf_t3_neg=zeros(ii_num,l);
        tf_t1_neg_L=zeros(ii_num,l); tf_t2_neg_L=zeros(ii_num,l); tf_t3_neg_L=zeros(ii_num,l);

        tf_t1_pos=zeros(ii_num,l); tf_t2_pos=zeros(ii_num,l); tf_t3_pos=zeros(ii_num,l);
        tf_t1_pos_L=zeros(ii_num,l); tf_t2_pos_L=zeros(ii_num,l); tf_t3_pos_L=zeros(ii_num,l);
    end


    tf_t1_neg(ii,:)=tf_diff_t1(:,5)'; tf_t1_pos(ii,:)=tf_diff_t1(:,195)';
    tf_t1_neg_L(ii,:)=tf_diff_t1_L(:,5)'; tf_t1_pos_L(ii,:)=tf_diff_t1_L(:,195)';

    tf_t2_neg(ii,:)=tf_diff_t2(:,5)'; tf_t2_pos(ii,:)=tf_diff_t2(:,195)';
    tf_t2_neg_L(ii,:)=tf_diff_t2_L(:,5)'; tf_t2_pos_L(ii,:)=tf_diff_t2_L(:,195)';

    tf_t3_neg(ii,:)=tf_diff_t3(:,5)'; tf_t3_pos(ii,:)=tf_diff_t3(:,195)';
    tf_t3_neg_L(ii,:)=tf_diff_t3_L(:,5)'; tf_t3_pos_L(ii,:)=tf_diff_t3_L(:,195)';

    clear tf_diff* tf_prox_diff* tf_dist_diff*
    %% 4) clear variables
    clearvars -except cond row* col* ii* *bin* *max name color outdir tf* output_name input_cond *series* deg* t l q_legends
end

%% 5) plot the bar chart for turn frequency in 3 different quadrants
data_tf={tf',tf_neg',tf_pos'};
data_tf_prox={tf_prox',tf_prox_neg',tf_prox_pos'};
data_tf_dist={tf_dist',tf_dist_neg',tf_dist_pos'};

data_tf_L={tf_L',tf_neg_L',tf_pos_L'};
data_tf_prox_L={tf_prox_L',tf_prox_neg_L',tf_prox_pos_L'};
data_tf_dist_L={tf_dist_L',tf_dist_neg_L',tf_dist_pos_L'};

data_tf_t1={tf_t1',tf_t1_neg',tf_t1_pos'};
data_tf_t2={tf_t2',tf_t2_neg',tf_t2_pos'};
data_tf_t3={tf_t3',tf_t3_neg',tf_t3_pos'};

data_tf_t1_L={tf_t1_L',tf_t1_neg_L',tf_t1_pos_L'};
data_tf_t2_L={tf_t2_L',tf_t2_neg_L',tf_t2_pos_L'};
data_tf_t3_L={tf_t3_L',tf_t3_neg_L',tf_t3_pos_L'};

mlt_subplt(data_tf,[],112,1,2,1,'',"Orientation ","Turn freq (#/min)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color);
mlt_subplt(data_tf_L,[],112,1,2,2,append("tf_bar_",num2str(output_name)),"Orientation ","Turn freq (#/min, >rms turn size)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);

mlt_subplt(data_tf_prox,[],113,1,2,1,'',"Orientation ","Turn freq (#/min, prox)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color);
mlt_subplt(data_tf_prox_L,[],113,1,2,2,append("tf_Proximate_bar_",num2str(output_name)),"Orientation ","Turn freq (#/min, >rms turn size, prox)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);

mlt_subplt(data_tf_dist,[],114,1,2,1,'',"Orientation ","Turn freq (#/min, dist)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color);
mlt_subplt(data_tf_dist_L,[],114,1,2,2,append("tf_Distant_bar_",num2str(output_name)),"Orientation ","Turn freq (#/min, >rms turn size, dist)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);

mlt_subplt(data_tf_t1,[],115,2,3,1,'',"Orientation ","Turn freq (#/min)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color,'title',"0-300s");
mlt_subplt(data_tf_t1_L,[],115,2,3,4,'',"Orientation ","Turn freq (#/min, >rms turn size)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);

mlt_subplt(data_tf_t2,[],115,2,3,2,'',"Orientation ",'',color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color,'title',"300-600s");
mlt_subplt(data_tf_t2_L,[],115,2,3,5,'',"Orientation ",'',color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);

mlt_subplt(data_tf_t3,[],115,2,3,3,'',"Orientation ",'',color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color,'title',"600-900s");
mlt_subplt(data_tf_t3_L,[],115,2,3,6,append("tf across time_bar_",num2str(output_name)),"Orientation ",'',color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);
%% SAVE FIGURES
outdir1=fullfile(outdir,"turn freq","rms turn size",append(num2str(l)," quadrants"));
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;
%% 6-1) stat test for ctr vs exp groups
if ~isempty(input_cond)
    idx1=find(contains(input_cond,"ctr"));
    %split data for bootstrap null hyppothesis testing
    [ctr_tf_all,exp_tf_all,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(tf_bt_all,idx1,name,color,'input_cond',input_cond);
    [ctr_tf_dist_all,exp_tf_dist_all]=split_data_ctr_exp(tf_bt_dist_all,idx1,name,color,'input_cond',input_cond);
    [ctr_tf_prox_all,exp_tf_prox_all]=split_data_ctr_exp(tf_bt_prox_all,idx1,name,color,'input_cond',input_cond);

    [ctr_tf_all_L,exp_tf_all_L]=split_data_ctr_exp(tf_bt_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_tf_dist_all_L,exp_tf_dist_all_L]=split_data_ctr_exp(tf_bt_dist_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_tf_prox_all_L,exp_tf_prox_all_L]=split_data_ctr_exp(tf_bt_prox_all_L,idx1,name,color,'input_cond',input_cond);

    [ctr_t1_all,exp_t1_all]=split_data_ctr_exp(tf_bt_t1_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_all,exp_t2_all]=split_data_ctr_exp(tf_bt_t2_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_all,exp_t3_all]=split_data_ctr_exp(tf_bt_t3_all,idx1,name,color,'input_cond',input_cond);

    [ctr_t1_all_L,exp_t1_all_L]=split_data_ctr_exp(tf_bt_t1_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_all_L,exp_t2_all_L]=split_data_ctr_exp(tf_bt_t2_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_all_L,exp_t3_all_L]=split_data_ctr_exp(tf_bt_t3_all_L,idx1,name,color,'input_cond',input_cond);
    for i=1:length(exp_tf_all)
        % get the color for the ctr and exp groups
        legends=[ctr_name(i,1);exp_name(i,1)];
        color1=[ctr_color(i,:);exp_color(i,:)];
        %% start to compute the bootstrap null hypothesis testing
        bt_diff=exp_tf_all{i,1}-ctr_tf_all{i,1};
        bt_diff_prox=exp_tf_prox_all{i,1}-ctr_tf_prox_all{i,1};
        bt_diff_dist=exp_tf_dist_all{i,1}-ctr_tf_dist_all{i,1};

        bt_diff_L=exp_tf_all_L{i,1}-ctr_tf_all_L{i,1};
        bt_diff_prox_L=exp_tf_prox_all_L{i,1}-ctr_tf_prox_all_L{i,1};
        bt_diff_dist_L=exp_tf_dist_all_L{i,1}-ctr_tf_dist_all_L{i,1};

        bt_diff_t1=exp_t1_all{i,1}-ctr_t1_all{i,1};
        bt_diff_t2=exp_t2_all{i,1}-ctr_t2_all{i,1};
        bt_diff_t3=exp_t3_all{i,1}-ctr_t3_all{i,1};

        bt_diff_t1_L=exp_t1_all_L{i,1}-ctr_t1_all_L{i,1};
        bt_diff_t2_L=exp_t2_all_L{i,1}-ctr_t2_all_L{i,1};
        bt_diff_t3_L=exp_t3_all_L{i,1}-ctr_t3_all_L{i,1};
        %%
        N=width(bt_diff);
        group_p=strings(0);
        data_p=[]; data_p_prox=[]; data_p_dist=[];
        data_p_L=[]; data_p_prox_L=[]; data_p_dist_L=[];

        data_p_t1=[]; data_p_t2=[]; data_p_t3=[];
        data_p_t1_L=[]; data_p_t2_L=[]; data_p_t3_L=[];

        for h=1:l
            group_p=vertcat(group_p,repmat(q_legends(h),N,1));

            data_p=vertcat(data_p,bt_diff(h,:)');
            data_p_prox=vertcat(data_p_prox,bt_diff_prox(h,:)');
            data_p_dist=vertcat(data_p_dist,bt_diff_dist(h,:)');

            data_p_L=vertcat(data_p_L,bt_diff_L(h,:)');
            data_p_prox_L=vertcat(data_p_prox_L,bt_diff_prox_L(h,:)');
            data_p_dist_L=vertcat(data_p_dist_L,bt_diff_dist_L(h,:)');

            data_p_t1=vertcat(data_p_t1,bt_diff_t1(h,:)');
            data_p_t2=vertcat(data_p_t2,bt_diff_t2(h,:)');
            data_p_t3=vertcat(data_p_t3,bt_diff_t3(h,:)');

            data_p_t1_L=vertcat(data_p_t1_L,bt_diff_t1_L(h,:)');
            data_p_t2_L=vertcat(data_p_t2_L,bt_diff_t2_L(h,:)');
            data_p_t3_L=vertcat(data_p_t3_L,bt_diff_t3_L(h,:)');
        end

        mlt_subplt(data_p,group_p,106,row_comp,col_comp,i,append("Bt test_tf_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1)),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_dist,group_p,107,row_comp,col_comp,i,append("Bt test_tf_distant_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_d_i_s_t"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_prox,group_p,108,row_comp,col_comp,i,append("Bt test_tf_proximal_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_p_r_o_x"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);

        mlt_subplt(data_p_L,group_p,109,row_comp,col_comp,i,append("Bt test_tf_larger than rms_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_dist_L,group_p,110,row_comp,col_comp,i,append("Bt test_tf_larger than rms_distant_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_d_i_s_t"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_prox_L,group_p,111,row_comp,col_comp,i,append("Bt test_tf_larger than rms_proximal_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_p_r_o_x"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

        mlt_subplt(data_p_t1,group_p,100,row_comp,col_comp,i,append("Bt test_tf_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_t2,group_p,101,row_comp,col_comp,i,append("Bt test_tf_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_t3,group_p,102,row_comp,col_comp,i,append("Bt test_tf_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);

        mlt_subplt(data_p_t1_L,group_p,103,row_comp,col_comp,i,append("Bt test_tf_larger than rms_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_t2_L,group_p,104,row_comp,col_comp,i,append("Bt test_tf_larger than rms_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_t3_L,group_p,105,row_comp,col_comp,i,append("Bt test_tf_larger than rms_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

    end

end
%% 6-2) Compute the bt hypothesis testing for each condition between different quadrants
str=strings(0);
comp=cell(ii_num,1);comp_prox=cell(ii_num,1);comp_dist=cell(ii_num,1);
comp_t1=cell(ii_num,1);comp_t2=cell(ii_num,1);comp_t3=cell(ii_num,1);
comp_t1_L=cell(ii_num,1);comp_t2_L=cell(ii_num,1);comp_t3_L=cell(ii_num,1);
comp_L=cell(ii_num,1);comp_dist_L=cell(ii_num,1);comp_prox_L=cell(ii_num,1);
for b=1:l-1
    for c=b+1:l
        for a=1:ii_num
            arry=tf_bt_all{a,1};arry_prox=tf_bt_prox_all{a,1};arry_dist=tf_bt_dist_all{a,1};
            arry_t1=tf_bt_t1_all{a,1};arry_t2=tf_bt_t2_all{a,1};arry_t3=tf_bt_t3_all{a,1};

            arry_L=tf_bt_all_L{a,1};arry_prox_L=tf_bt_prox_all_L{a,1};arry_dist_L=tf_bt_dist_all_L{a,1};
            arry_t1_L=tf_bt_t1_all_L{a,1};arry_t2_L=tf_bt_t2_all_L{a,1};arry_t3_L=tf_bt_t3_all_L{a,1};

            comp{a,1}=vertcat(comp{a,1},(arry(c,:)-arry(b,:))');
            comp_prox{a,1}=vertcat(comp_prox{a,1},(arry_prox(c,:)-arry_prox(b,:))');
            comp_dist{a,1}=vertcat(comp_dist{a,1},(arry_prox(c,:)-arry_dist(b,:))');

            comp_t1{a,1}=vertcat(comp_t1{a,1},(arry_t1(c,:)-arry_t1(b,:))');
            comp_t2{a,1}=vertcat(comp_t2{a,1},(arry_t2(c,:)-arry_t2(b,:))');
            comp_t3{a,1}=vertcat(comp_t3{a,1},(arry_t3(c,:)-arry_t3(b,:))');

            comp_L{a,1}=vertcat(comp_L{a,1},(arry_L(c,:)-arry_L(b,:))');
            comp_prox_L{a,1}=vertcat(comp_prox_L{a,1},(arry_prox_L(c,:)-arry_prox_L(b,:))');
            comp_dist_L{a,1}=vertcat(comp_dist_L{a,1},(arry_dist_L(c,:)-arry_dist_L(b,:))');

            comp_t1_L{a,1}=vertcat(comp_t1_L{a,1},(arry_t1_L(c,:)-arry_t1_L(b,:))');
            comp_t2_L{a,1}=vertcat(comp_t2_L{a,1},(arry_t2_L(c,:)-arry_t2_L(b,:))');
            comp_t3_L{a,1}=vertcat(comp_t3_L{a,1},(arry_t3_L(c,:)-arry_t3_L(b,:))');

        end
        str=vertcat(str,repmat(append(q_legends(c),"-",q_legends(b)),width(arry),1));

    end
end
%plot the bt test as a box plot
for a=1:ii_num
    mlt_subplt(comp{a,1},str,200,row,col,a,append("Bt test_tf_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1)),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    mlt_subplt(comp_prox{a,1},str,210 ,row,col,a,append("Bt test_tf_proximal_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_p_r_o_x_i_m_a_l"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    mlt_subplt(comp_dist{a,1},str,220 ,row,col,a,append("Bt test_tf_distant_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_d_i_s_t_a_n_t"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

    mlt_subplt(comp_t1{a,1},str,230 ,row,col,a,append("Bt test_tf_0-300s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    mlt_subplt(comp_t2{a,1},str,240 ,row,col,a,append("Bt test_tf_300-600s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    mlt_subplt(comp_t3{a,1},str,250 ,row,col,a,append("Bt test_tf_600-900s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

    mlt_subplt(comp_L{a,1},str,260 ,row,col,a,append("Bt test_tf_larger than rms_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_p_r_o_x_i_m_a_l_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    mlt_subplt(comp_prox_L{a,1},str,270 ,row,col,a,append("Bt test_tf_larger than rms_proximal_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_p_r_o_x_i_m_a_l"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    mlt_subplt(comp_dist_L{a,1},str,280 ,row,col,a,append("Bt test_tf_larger than rms_distant_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_d_i_s_t_a_n_t"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

    mlt_subplt(comp_t1_L{a,1},str,290 ,row,col,a,append("Bt test_tf_larger than rms_0-300s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    mlt_subplt(comp_t2_L{a,1},str,300 ,row,col,a,append("Bt test_tf_larger than rms_300-600s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    mlt_subplt(comp_t3_L{a,1},str,310 ,row,col,a,append("Bt test_tf_larger than rms_600-900s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

end

%% SAVE FIGURES
outdir1=fullfile(outdir,"turn freq","rms turn size",append(num2str(l)," quadrants"),"stat");
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;
end