function []=turn_freq_t_pos_4_quadrants(tracker_num,genotype,condition,filename,output_name,varargin)
%% this function will compute the turning freq per min at the prox vs distant region across different duration (0-300s, 300-600s, and 600-900s)

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
xmax=225; 
deg_list=60;% this is used for the cutoff for the turn size -->can be an arry (size ii_num*1) or a number 
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
    %% count the number of turning events

    if ii==1

        tf_t1_prox=zeros(ii_num,l);  tf_t2_prox=zeros(ii_num,l);  tf_t3_prox=zeros(ii_num,l);
        tf_t1_dist=zeros(ii_num,l);  tf_t2_dist=zeros(ii_num,l);  tf_t3_dist=zeros(ii_num,l);

        tf_t1_prox_L=zeros(ii_num,l);  tf_t2_prox_L=zeros(ii_num,l);  tf_t3_prox_L=zeros(ii_num,l);
        tf_t1_dist_L=zeros(ii_num,l);  tf_t2_dist_L=zeros(ii_num,l);  tf_t3_dist_L=zeros(ii_num,l);
    end
    t=[0:tmax/3:tmax]';% split it in three time sections
    t3=bin_data_sum_t3(dat.orientation,dat.x,dat.et,hdseries,[0:xbin2:xmax]',t,dat.et);
    c=bin_data_count3(pre,hdseries,turn_x,[0:xbin2:xmax]',t0s,t);

    if l==4
        t3_prox=assign_num_4quadrants(squeeze(t3(:,2,:)));
        t3_dist=assign_num_4quadrants(squeeze(t3(:,1,:)));

        c_prox=assign_num_4quadrants(squeeze(c(:,2,:)));
        c_dist=assign_num_4quadrants(squeeze(c(:,1,:)));
    elseif l==3
        t3_prox=assign_num_3quadrants(squeeze(t3(:,2,:)));
        t3_dist=assign_num_3quadrants(squeeze(t3(:,1,:)));

        c_prox=assign_num_3quadrants(squeeze(c(:,2,:)));
        c_dist=assign_num_3quadrants(squeeze(c(:,1,:)));
    end


    tf_prox=(c_prox./t3_prox).*60;
    tf_dist=(c_dist./t3_dist).*60;
    tf_t1_prox(ii,:)=tf_prox(:,1); tf_t2_prox(ii,:)=tf_prox(:,2); tf_t3_prox(ii,:)=tf_prox(:,3);
    tf_t1_dist(ii,:)=tf_dist(:,1); tf_t2_dist(ii,:)=tf_dist(:,2); tf_t3_dist(ii,:)=tf_dist(:,3);

    clear c_prox c_dist tf_dist tf_prox c
    %% 1-3) Compute the turning freq for turning event larger than 60 deg
    c_L=bin_data_count3(pre(idx),hdseries,turn_x(idx),[0:xbin2:xmax]',t0s(idx),t);
    
    if l==4
        c_prox_L=assign_num_4quadrants(squeeze(c_L(:,2,:)));
        c_dist_L=assign_num_4quadrants(squeeze(c_L(:,1,:)));
    elseif l==3
        c_prox_L=assign_num_3quadrants(squeeze(c_L(:,2,:)));
        c_dist_L=assign_num_3quadrants(squeeze(c_L(:,1,:)));
    end

    tf_prox_L=(c_prox_L./t3_prox).*60;
    tf_dist_L=(c_dist_L./t3_dist).*60;

    tf_t1_prox_L(ii,:)=tf_prox_L(:,1); tf_t2_prox_L(ii,:)=tf_prox_L(:,2); tf_t3_prox_L(ii,:)=tf_prox_L(:,3);
    tf_t1_dist_L(ii,:)=tf_dist_L(:,1); tf_t2_dist_L(ii,:)=tf_dist_L(:,2); tf_t3_dist_L(ii,:)=tf_dist_L(:,3);

    clear c_L c_prox_L c_dist_L tf_prox_L tf_dist_L
    %% 1-4) Use the bootstrapping to get the 95% CI
    len=200;
    tf_bt_t1_prox=zeros(l,len); tf_bt_t2_prox=zeros(l,len); tf_bt_t3_prox=zeros(l,len);
    tf_bt_t1_prox_L=zeros(l,len); tf_bt_t2_prox_L=zeros(l,len); tf_bt_t3_prox_L=zeros(l,len);

    tf_bt_t1_dist=zeros(l,len); tf_bt_t2_dist=zeros(l,len); tf_bt_t3_dist=zeros(l,len);
    tf_bt_t1_dist_L=zeros(l,len); tf_bt_t2_dist_L=zeros(l,len); tf_bt_t3_dist_L=zeros(l,len);

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
        t0_bt=vertcat(data.t0s{:});
        if iscell(t0_bt)

            t0s_bt=vertcat(t0_bt{:});
        else
            t0s_bt=t0_bt;
        end
        idx_bt=find(reori_deg_bt>=deg_l);

        % bin the turning events based on their pre-turning heading
        % directions
        t3_bt=bin_data_sum_t3(data.orientation,data.x,data.et,hdseries,[0:xbin2:xmax]',t,data.et);
        c_bt=bin_data_count3(pre_bt,hdseries,turn_x_bt,[0:xbin2:xmax]',t0s_bt,t);

        if l==4
            t3_bt_prox=assign_num_4quadrants(squeeze(t3_bt(:,2,:)));
            t3_bt_dist=assign_num_4quadrants(squeeze(t3_bt(:,1,:)));

            c_bt_prox=assign_num_4quadrants(squeeze(c_bt(:,2,:)));
            c_bt_dist=assign_num_4quadrants(squeeze(c_bt(:,1,:)));
        elseif l==3
            t3_bt_prox=assign_num_3quadrants(squeeze(t3_bt(:,2,:)));
            t3_bt_dist=assign_num_3quadrants(squeeze(t3_bt(:,1,:)));

            c_bt_prox=assign_num_3quadrants(squeeze(c_bt(:,2,:)));
            c_bt_dist=assign_num_3quadrants(squeeze(c_bt(:,1,:)));
        end
        

        tf_bt_prox=(c_bt_prox./t3_bt_prox).*60;
        tf_bt_dist=(c_bt_dist./t3_bt_dist).*60;
        tf_bt_t1_prox(:,k)=tf_bt_prox(:,1); tf_bt_t2_prox(:,k)=tf_bt_prox(:,2); tf_bt_t3_prox(:,k)=tf_bt_prox(:,3);
        tf_bt_t1_dist(:,k)=tf_bt_dist(:,1); tf_bt_t2_dist(:,k)=tf_bt_dist(:,2); tf_bt_t3_dist(:,k)=tf_bt_dist(:,3);

        c_bt_L=bin_data_count3(pre_bt(idx_bt),hdseries,turn_x_bt(idx_bt),[0:xbin2:xmax]',t0s_bt(idx_bt),t);

        if l==4
            c_bt_prox_L=assign_num_4quadrants(squeeze(c_bt_L(:,2,:)));
            c_bt_dist_L=assign_num_4quadrants(squeeze(c_bt_L(:,1,:)));
        elseif l==3
            c_bt_prox_L=assign_num_3quadrants(squeeze(c_bt_L(:,2,:)));
            c_bt_dist_L=assign_num_3quadrants(squeeze(c_bt_L(:,1,:)));
        end

        tf_bt_prox_L=(c_bt_prox_L./t3_bt_prox).*60;
        tf_bt_dist_L=(c_bt_dist_L./t3_bt_dist).*60;

        tf_bt_t1_prox_L(:,k)=tf_bt_prox_L(:,1); tf_bt_t2_prox_L(:,k)=tf_bt_prox_L(:,2); tf_bt_t3_prox_L(:,k)=tf_bt_prox_L(:,3);
        tf_bt_t1_dist_L(:,k)=tf_bt_dist_L(:,1); tf_bt_t2_dist_L(:,k)=tf_bt_dist_L(:,2); tf_bt_t3_dist_L(:,k)=tf_bt_dist_L(:,3);


        clear c_bt* t3_bt* tf_bt_prox* tf_dist*
    end


    % Save data for computing the bootstrap p value
    if ii==1
        tf_bt_t1_prox_all=cell(ii_num,1); tf_bt_t2_prox_all=cell(ii_num,1);tf_bt_t3_prox_all=cell(ii_num,1);
        tf_bt_t1_dist_all=cell(ii_num,1); tf_bt_t2_dist_all=cell(ii_num,1);tf_bt_t3_dist_all=cell(ii_num,1);

        tf_bt_t1_prox_all_L=cell(ii_num,1); tf_bt_t2_prox_all_L=cell(ii_num,1);tf_bt_t3_prox_all_L=cell(ii_num,1);
        tf_bt_t1_dist_all_L=cell(ii_num,1); tf_bt_t2_dist_all_L=cell(ii_num,1);tf_bt_t3_dist_all_L=cell(ii_num,1);

    end

    tf_bt_t1_prox_all{ii,1}=tf_bt_t1_prox; tf_bt_t2_prox_all{ii,1}=tf_bt_t2_prox; tf_bt_t3_prox_all{ii,1}=tf_bt_t3_prox;
    tf_bt_t1_dist_all{ii,1}=tf_bt_t1_dist; tf_bt_t2_dist_all{ii,1}=tf_bt_t2_dist; tf_bt_t3_dist_all{ii,1}=tf_bt_t3_dist;

    tf_bt_t1_prox_all_L{ii,1}=tf_bt_t1_prox_L; tf_bt_t2_prox_all_L{ii,1}=tf_bt_t2_prox_L; tf_bt_t3_prox_all_L{ii,1}=tf_bt_t3_prox_L;
    tf_bt_t1_dist_all_L{ii,1}=tf_bt_t1_dist_L; tf_bt_t2_dist_all_L{ii,1}=tf_bt_t2_dist_L; tf_bt_t3_dist_all_L{ii,1}=tf_bt_t3_dist_L;
    %% Use the data to compute the 95% CI
    % for distant and proximate regions
    tf_t1_prox_diff=tf_bt_t1_prox-tf_t1_prox(ii,:)';tf_t2_prox_diff=tf_bt_t2_prox-tf_t2_prox(ii,:)';tf_t3_prox_diff=tf_bt_t3_prox-tf_t3_prox(ii,:)';
    tf_t1_dist_diff=tf_bt_t1_dist-tf_t1_dist(ii,:)';tf_t2_dist_diff=tf_bt_t2_dist-tf_t2_dist(ii,:)';tf_t3_dist_diff=tf_bt_t3_dist-tf_t3_dist(ii,:)';

    tf_t1_prox_diff_L=tf_bt_t1_prox_L-tf_t1_prox_L(ii,:)';tf_t2_prox_diff_L=tf_bt_t2_prox_L-tf_t2_prox_L(ii,:)';tf_t3_prox_diff_L=tf_bt_t3_prox_L-tf_t3_prox_L(ii,:)';
    tf_t1_dist_diff_L=tf_bt_t1_dist_L-tf_t1_dist_L(ii,:)';tf_t2_dist_diff_L=tf_bt_t2_dist_L-tf_t2_dist_L(ii,:)';tf_t3_dist_diff_L=tf_bt_t3_dist_L-tf_t3_dist_L(ii,:)';

    %% sort the tf from small to large and get the 95% CI
    tf_t1_prox_diff=sort(tf_t1_prox_diff,2); tf_t2_prox_diff=sort(tf_t2_prox_diff,2); tf_t3_prox_diff=sort(tf_t3_prox_diff,2);
    tf_t1_dist_diff=sort(tf_t1_dist_diff,2); tf_t2_dist_diff=sort(tf_t2_dist_diff,2); tf_t3_dist_diff=sort(tf_t3_dist_diff,2);

    tf_t1_prox_diff_L=sort(tf_t1_prox_diff_L,2); tf_t2_prox_diff_L=sort(tf_t2_prox_diff_L,2); tf_t3_prox_diff_L=sort(tf_t3_prox_diff_L,2);
    tf_t1_dist_diff_L=sort(tf_t1_dist_diff_L,2); tf_t2_dist_diff_L=sort(tf_t2_dist_diff_L,2); tf_t3_dist_diff_L=sort(tf_t3_dist_diff_L,2);
    if ii==1
        tf_t1_prox_pos=zeros(ii_num,l); tf_t2_prox_pos=zeros(ii_num,l); tf_t3_prox_pos=zeros(ii_num,l);
        tf_t1_prox_neg=zeros(ii_num,l); tf_t2_prox_neg=zeros(ii_num,l); tf_t3_prox_neg=zeros(ii_num,l);

        tf_t1_dist_pos=zeros(ii_num,l); tf_t2_dist_pos=zeros(ii_num,l); tf_t3_dist_pos=zeros(ii_num,l);
        tf_t1_dist_neg=zeros(ii_num,l); tf_t2_dist_neg=zeros(ii_num,l); tf_t3_dist_neg=zeros(ii_num,l);

        tf_t1_prox_pos_L=zeros(ii_num,l); tf_t2_prox_pos_L=zeros(ii_num,l); tf_t3_prox_pos_L=zeros(ii_num,l);
        tf_t1_prox_neg_L=zeros(ii_num,l); tf_t2_prox_neg_L=zeros(ii_num,l); tf_t3_prox_neg_L=zeros(ii_num,l);

        tf_t1_dist_pos_L=zeros(ii_num,l); tf_t2_dist_pos_L=zeros(ii_num,l); tf_t3_dist_pos_L=zeros(ii_num,l);
        tf_t1_dist_neg_L=zeros(ii_num,l); tf_t2_dist_neg_L=zeros(ii_num,l); tf_t3_dist_neg_L=zeros(ii_num,l);
    end

    tf_t1_prox_neg(ii,:)=tf_t1_prox_diff(:,5)'; tf_t1_prox_pos(ii,:)=tf_t1_prox_diff(:,195)';
    tf_t2_prox_neg(ii,:)=tf_t2_prox_diff(:,5)'; tf_t2_prox_pos(ii,:)=tf_t2_prox_diff(:,195)';
    tf_t3_prox_neg(ii,:)=tf_t3_prox_diff(:,5)'; tf_t3_prox_pos(ii,:)=tf_t3_prox_diff(:,195)';

    tf_t1_dist_neg(ii,:)=tf_t1_dist_diff(:,5)'; tf_t1_dist_pos(ii,:)=tf_t1_dist_diff(:,195)';
    tf_t2_dist_neg(ii,:)=tf_t2_dist_diff(:,5)'; tf_t2_dist_pos(ii,:)=tf_t2_dist_diff(:,195)';
    tf_t3_dist_neg(ii,:)=tf_t3_dist_diff(:,5)'; tf_t3_dist_pos(ii,:)=tf_t3_dist_diff(:,195)';

    tf_t1_prox_neg_L(ii,:)=tf_t1_prox_diff_L(:,5)'; tf_t1_prox_pos_L(ii,:)=tf_t1_prox_diff_L(:,195)';
    tf_t2_prox_neg_L(ii,:)=tf_t2_prox_diff_L(:,5)'; tf_t2_prox_pos_L(ii,:)=tf_t2_prox_diff_L(:,195)';
    tf_t3_prox_neg_L(ii,:)=tf_t3_prox_diff_L(:,5)'; tf_t3_prox_pos_L(ii,:)=tf_t3_prox_diff_L(:,195)';

    tf_t1_dist_neg_L(ii,:)=tf_t1_dist_diff_L(:,5)'; tf_t1_dist_pos_L(ii,:)=tf_t1_dist_diff_L(:,195)';
    tf_t2_dist_neg_L(ii,:)=tf_t2_dist_diff_L(:,5)'; tf_t2_dist_pos_L(ii,:)=tf_t2_dist_diff_L(:,195)';
    tf_t3_dist_neg_L(ii,:)=tf_t3_dist_diff_L(:,5)'; tf_t3_dist_pos_L(ii,:)=tf_t3_dist_diff_L(:,195)';
    clear tf_*_prox_diff* tf_*_dist_diff*
    %% 2) clear variables
    clearvars -except cond row* col* ii* *bin* *max name color outdir tf* output_name input_cond *series* deg* t l len q_legends
end

%plot the figures
data_tf_t1_prox={tf_t1_prox',tf_t1_prox_neg',tf_t1_prox_pos'};
data_tf_t2_prox={tf_t2_prox',tf_t2_prox_neg',tf_t2_prox_pos'};
data_tf_t3_prox={tf_t3_prox',tf_t3_prox_neg',tf_t3_prox_pos'};

data_tf_t1_dist={tf_t1_dist',tf_t1_dist_neg',tf_t1_dist_pos'};
data_tf_t2_dist={tf_t2_dist',tf_t2_dist_neg',tf_t2_dist_pos'};
data_tf_t3_dist={tf_t3_dist',tf_t3_dist_neg',tf_t3_dist_pos'};

data_tf_t1_prox_L={tf_t1_prox_L',tf_t1_prox_neg_L',tf_t1_prox_pos_L'};
data_tf_t2_prox_L={tf_t2_prox_L',tf_t2_prox_neg_L',tf_t2_prox_pos_L'};
data_tf_t3_prox_L={tf_t3_prox_L',tf_t3_prox_neg_L',tf_t3_prox_pos_L'};

data_tf_t1_dist_L={tf_t1_dist_L',tf_t1_dist_neg_L',tf_t1_dist_pos_L'};
data_tf_t2_dist_L={tf_t2_dist_L',tf_t2_dist_neg_L',tf_t2_dist_pos_L'};
data_tf_t3_dist_L={tf_t3_dist_L',tf_t3_dist_neg_L',tf_t3_dist_pos_L'};

mlt_subplt(data_tf_t1_prox,[],100,2,3,1,'',"","Turn freq (#/min,prox)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color,'title',"0-300s");
mlt_subplt(data_tf_t2_prox,[],100,2,3,2,'',""," ",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color,'title',"300-600s");
mlt_subplt(data_tf_t3_prox,[],100,2,3,3,'',""," ",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color,'title',"600-900s");

mlt_subplt(data_tf_t1_prox_L,[],100,2,3,4,'',"Orientation ","Turn freq (#/min,>rms turn size,prox)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);
mlt_subplt(data_tf_t2_prox_L,[],100,2,3,5,'',"Orientation "," ",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);
mlt_subplt(data_tf_t3_prox_L,[],100,2,3,6,append("tf_Proximate_across time_bar_",num2str(output_name)),"Orientation "," ",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);

mlt_subplt(data_tf_t1_dist,[],101,2,3,1,'',"","Turn freq (#/min,dist)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color,'title',"0-300s");
mlt_subplt(data_tf_t2_dist,[],101,2,3,2,'',""," ",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color,'title',"300-600s");
mlt_subplt(data_tf_t3_dist,[],101,2,3,3,'',""," ",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 7],'bar_color',color,'title',"600-900s");

mlt_subplt(data_tf_t1_dist_L,[],101,2,3,4,'',"Orientation ","Turn freq (#/min,>rms turn size,dist)",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);
mlt_subplt(data_tf_t2_dist_L,[],101,2,3,5,'',"Orientation "," ",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'y lim',[0 3],'bar_color',color);
mlt_subplt(data_tf_t3_dist_L,[],101,2,3,6,append("tf_Distant_across time_bar_",num2str(output_name)),"Orientation "," ",color,q_legends,"bar",'ii_num',ii_num,'bt_CI',1,'legends',name,'ylim',[0 3],'bar_color',color);
%% SAVE FIGURES
outdir1=fullfile(outdir,"turn freq","rms turn size",append(num2str(l)," quadrants"));
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;

%% perfomr the stat test
%% 3-1) this will compute the bt difference between ctr and each exp group for each quadrant

if ~isempty(input_cond)
    idx1=find(contains(input_cond,"ctr"));
    %split data for bootstrap null hyppothesis testing
    [ctr_t1_prox,exp_t1_prox,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(tf_bt_t1_prox_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_prox,exp_t2_prox]=split_data_ctr_exp(tf_bt_t2_prox_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_prox,exp_t3_prox]=split_data_ctr_exp(tf_bt_t3_prox_all,idx1,name,color,'input_cond',input_cond);

    [ctr_t1_dist,exp_t1_dist]=split_data_ctr_exp(tf_bt_t1_dist_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_dist,exp_t2_dist]=split_data_ctr_exp(tf_bt_t2_dist_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_dist,exp_t3_dist]=split_data_ctr_exp(tf_bt_t3_dist_all,idx1,name,color,'input_cond',input_cond);

    [ctr_t1_prox_L,exp_t1_prox_L]=split_data_ctr_exp(tf_bt_t1_prox_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_prox_L,exp_t2_prox_L]=split_data_ctr_exp(tf_bt_t2_prox_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_prox_L,exp_t3_prox_L]=split_data_ctr_exp(tf_bt_t3_prox_all_L,idx1,name,color,'input_cond',input_cond);

    [ctr_t1_dist_L,exp_t1_dist_L]=split_data_ctr_exp(tf_bt_t1_dist_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_dist_L,exp_t2_dist_L]=split_data_ctr_exp(tf_bt_t2_dist_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_dist_L,exp_t3_dist_L]=split_data_ctr_exp(tf_bt_t3_dist_all_L,idx1,name,color,'input_cond',input_cond);
    for i=1:length(ctr_t1_prox)
        % get the color for the ctr and exp groups
        legends=[ctr_name(i,1);exp_name(i,1)];
        color1=[ctr_color(i,:);exp_color(i,:)];
        %% start to compute the bootstrap null hypothesis testing
        t1_prox_diff=exp_t1_prox{i,1}-ctr_t1_prox{i,1};
        t2_prox_diff=exp_t2_prox{i,1}-ctr_t2_prox{i,1};
        t3_prox_diff=exp_t3_prox{i,1}-ctr_t3_prox{i,1};

        t1_dist_diff=exp_t1_dist{i,1}-ctr_t1_dist{i,1};
        t2_dist_diff=exp_t2_dist{i,1}-ctr_t2_dist{i,1};
        t3_dist_diff=exp_t3_dist{i,1}-ctr_t3_dist{i,1};

        t1_prox_diff_L=exp_t1_prox_L{i,1}-ctr_t1_prox_L{i,1};
        t2_prox_diff_L=exp_t2_prox_L{i,1}-ctr_t2_prox_L{i,1};
        t3_prox_diff_L=exp_t3_prox_L{i,1}-ctr_t3_prox_L{i,1};

        t1_dist_diff_L=exp_t1_dist_L{i,1}-ctr_t1_dist_L{i,1};
        t2_dist_diff_L=exp_t2_dist_L{i,1}-ctr_t2_dist_L{i,1};
        t3_dist_diff_L=exp_t3_dist_L{i,1}-ctr_t3_dist_L{i,1};
        %% Compute the bootstraping hypothesis testing
        N=width(t1_prox_diff);
        group_p=strings(0);
        data_p_t1_prox=[]; data_p_t2_prox=[]; data_p_t3_prox=[];
        data_p_t1_dist=[]; data_p_t2_dist=[]; data_p_t3_dist=[];

        data_p_t1_prox_L=[]; data_p_t2_prox_L=[]; data_p_t3_prox_L=[];
        data_p_t1_dist_L=[]; data_p_t2_dist_L=[]; data_p_t3_dist_L=[];
        for h=1:l
            group_p=vertcat(group_p,repmat(q_legends(h),N,1));

            data_p_t1_prox=vertcat(data_p_t1_prox,t1_prox_diff(h,:)');
            data_p_t2_prox=vertcat(data_p_t2_prox,t2_prox_diff(h,:)');
            data_p_t3_prox=vertcat(data_p_t3_prox,t3_prox_diff(h,:)');

            data_p_t1_dist=vertcat(data_p_t1_dist,t1_dist_diff(h,:)');
            data_p_t2_dist=vertcat(data_p_t2_dist,t2_dist_diff(h,:)');
            data_p_t3_dist=vertcat(data_p_t3_dist,t3_dist_diff(h,:)');

            data_p_t1_prox_L=vertcat(data_p_t1_prox_L,t1_prox_diff_L(h,:)');
            data_p_t2_prox_L=vertcat(data_p_t2_prox_L,t2_prox_diff_L(h,:)');
            data_p_t3_prox_L=vertcat(data_p_t3_prox_L,t3_prox_diff_L(h,:)');

            data_p_t1_dist_L=vertcat(data_p_t1_dist_L,t1_dist_diff_L(h,:)');
            data_p_t2_dist_L=vertcat(data_p_t2_dist_L,t2_dist_diff_L(h,:)');
            data_p_t3_dist_L=vertcat(data_p_t3_dist_L,t3_dist_diff_L(h,:)');

        end
        mlt_subplt(data_p_t1_prox,group_p,104,row_comp,col_comp,i,append("Bt test_tf_proximal_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_p_r_o_x_i_m_a_l_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_t2_prox,group_p,105,row_comp,col_comp,i,append("Bt test_tf_proximal_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_p_r_o_x_i_m_a_l_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_t3_prox,group_p,106,row_comp,col_comp,i,append("Bt test_tf_proximal_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_p_r_o_x_i_m_a_l_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);

        mlt_subplt(data_p_t1_dist,group_p,107,row_comp,col_comp,i,append("Bt test_tf_distant_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_d_i_s_t_a_n_t_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_t2_dist,group_p,108,row_comp,col_comp,i,append("Bt test_tf_distant_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_d_i_s_t_a_n_t_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_t3_dist,group_p,108,row_comp,col_comp,i,append("Bt test_tf_distant_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_d_i_s_t_a_n_t_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);

        mlt_subplt(data_p_t1_prox_L,group_p,109,row_comp,col_comp,i,append("Bt test_tf_larger than rms_proximal_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_p_r_o_x_i_m_a_l_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_t2_prox_L,group_p,110,row_comp,col_comp,i,append("Bt test_tf_larger than rms_proximal_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_p_r_o_x_i_m_a_l_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_t3_prox_L,group_p,111,row_comp,col_comp,i,append("Bt test_tf_larger than rms_proximal_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_p_r_o_x_i_m_a_l_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

        mlt_subplt(data_p_t1_dist_L,group_p,112,row_comp,col_comp,i,append("Bt test_tf_larger than rms_distant_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_d_i_s_t_a_n_t_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_t2_dist_L,group_p,113,row_comp,col_comp,i,append("Bt test_tf_larger than rms_distant_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_d_i_s_t_a_n_t_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_t3_dist_L,group_p,114,row_comp,col_comp,i,append("Bt test_tf_larger than rms_distant_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Orientation ","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_l_a_r_g_e_,_d_i_s_t_a_n_t_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

    end
end
%% 3-2) this will perform the bt hypothesis test for different quadrants within the same condition
str=strings(0);
comp_t1_prox=cell(ii_num,1);comp_t2_prox=cell(ii_num,1);comp_t3_prox=cell(ii_num,1);
comp_t1_dist=cell(ii_num,1);comp_t2_dist=cell(ii_num,1);comp_t3_dist=cell(ii_num,1);
comp_t1_prox_L=cell(ii_num,1);comp_t2_prox_L=cell(ii_num,1);comp_t3_prox_L=cell(ii_num,1);
comp_t1_dist_L=cell(ii_num,1);comp_t2_dist_L=cell(ii_num,1);comp_t3_dist_L=cell(ii_num,1);
for b=1:l-1
    for c=b+1:l
        for a=1:ii_num
            arry_t1_prox=tf_bt_t1_prox_all{a,1};arry_t2_prox=tf_bt_t2_prox_all{a,1};arry_t3_prox=tf_bt_t3_prox_all{a,1};
            arry_t1_dist=tf_bt_t1_dist_all{a,1};arry_t2_dist=tf_bt_t2_dist_all{a,1};arry_t3_dist=tf_bt_t3_dist_all{a,1};

            arry_t1_prox_L=tf_bt_t1_prox_all_L{a,1};arry_t2_prox_L=tf_bt_t2_prox_all_L{a,1};arry_t3_prox_L=tf_bt_t3_prox_all_L{a,1};
            arry_t1_dist_L=tf_bt_t1_dist_all_L{a,1};arry_t2_dist_L=tf_bt_t2_dist_all_L{a,1};arry_t3_dist_L=tf_bt_t3_dist_all_L{a,1};

            comp_t1_prox{a,1}=vertcat(comp_t1_prox{a,1},(arry_t1_prox(c,:)-arry_t1_prox(b,:))');
            comp_t2_prox{a,1}=vertcat(comp_t2_prox{a,1},(arry_t2_prox(c,:)-arry_t2_prox(b,:))');
            comp_t3_prox{a,1}=vertcat(comp_t3_prox{a,1},(arry_t3_prox(c,:)-arry_t3_prox(b,:))');

            comp_t1_dist{a,1}=vertcat(comp_t1_dist{a,1},(arry_t1_dist(c,:)-arry_t1_dist(b,:))');
            comp_t2_dist{a,1}=vertcat(comp_t2_dist{a,1},(arry_t2_dist(c,:)-arry_t2_dist(b,:))');
            comp_t3_dist{a,1}=vertcat(comp_t3_dist{a,1},(arry_t3_dist(c,:)-arry_t3_dist(b,:))');

            comp_t1_prox_L{a,1}=vertcat(comp_t1_prox_L{a,1},(arry_t1_prox_L(c,:)-arry_t1_prox_L(b,:))');
            comp_t2_prox_L{a,1}=vertcat(comp_t2_prox_L{a,1},(arry_t2_prox_L(c,:)-arry_t2_prox_L(b,:))');
            comp_t3_prox_L{a,1}=vertcat(comp_t3_prox_L{a,1},(arry_t3_prox_L(c,:)-arry_t3_prox_L(b,:))');

            comp_t1_dist_L{a,1}=vertcat(comp_t1_dist_L{a,1},(arry_t1_dist_L(c,:)-arry_t1_dist_L(b,:))');
            comp_t2_dist_L{a,1}=vertcat(comp_t2_dist_L{a,1},(arry_t2_dist_L(c,:)-arry_t2_dist_L(b,:))');
            comp_t3_dist_L{a,1}=vertcat(comp_t3_dist_L{a,1},(arry_t3_dist_L(c,:)-arry_t3_dist_L(b,:))');

            
        end
        str=vertcat(str,repmat(append(q_legends(c),"-",q_legends(b)),width(arry_t1_prox),1));
    end
end
%plot the bt test as a box plot
for a=1:ii_num

     mlt_subplt(comp_t1_prox{a,1},str,114,row,col,a,append("Bt test_tf_proximal_0-300s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_p_r_o_x_i_m_a_l_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
     mlt_subplt(comp_t2_prox{a,1},str,124,row,col,a,append("Bt test_tf_proximal_300-600s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_p_r_o_x_i_m_a_l_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
     mlt_subplt(comp_t3_prox{a,1},str,134,row,col,a,append("Bt test_tf_proximal_600-900s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_p_r_o_x_i_m_a_l_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    
     mlt_subplt(comp_t1_dist{a,1},str,144,row,col,a,append("Bt test_tf_distant_0-300s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_d_i_s_t_a_n_t_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
     mlt_subplt(comp_t2_dist{a,1},str,154,row,col,a,append("Bt test_tf_distant_300-600s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_d_i_s_t_a_n_t_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
     mlt_subplt(comp_t3_dist{a,1},str,164,row,col,a,append("Bt test_tf_distant_600-900s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_d_i_s_t_a_n_t_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    
     mlt_subplt(comp_t1_prox_L{a,1},str,174,row,col,a,append("Bt test_tf_larger than rms_proximal_0-300s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_p_r_o_x_i_m_a_l_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
     mlt_subplt(comp_t2_prox_L{a,1},str,184,row,col,a,append("Bt test_tf_larger than rms_proximal_300-600s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_p_r_o_x_i_m_a_l_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
     mlt_subplt(comp_t3_prox_L{a,1},str,194,row,col,a,append("Bt test_tf_larger than rms_proximal_600-900s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_p_r_o_x_i_m_a_l_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
    
     mlt_subplt(comp_t1_dist_L{a,1},str,204,row,col,a,append("Bt test_tf_larger than rms_distant_0-300s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_d_i_s_t_a_n_t_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
     mlt_subplt(comp_t2_dist_L{a,1},str,214,row,col,a,append("Bt test_tf_larger than rms_distant_300-600s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_d_i_s_t_a_n_t_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
     mlt_subplt(comp_t3_dist_L{a,1},str,224,row,col,a,append("Bt test_tf_larger than rms_distant_600-900s_",num2str(output_name),"_3q boxplot"),"Orientation ","Difference",[],append(name(a,1),"_l_a_r_g_e_,_d_i_s_t_a_n_t_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

end 
%% SAVE FIGURES
outdir1=fullfile(outdir,"turn freq","rms turn size",append(num2str(l)," quadrants"),"stat");
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;
clear
end
