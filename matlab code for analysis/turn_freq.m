function []=turn_freq(tracker_num,genotype,condition,filename,output_name,varargin)
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
hbin=10;% for binning data
input_cond=[];%it will have a same length as name, telling us whether the input is 'ctr' or 'exp'
row_comp=2; col_comp=2; % these two will be used for the plot to show the comparison of ctr vs exp group, this number should be based on how many comparison you have
tmax=900;
xmax=225; deg_list=60;
%deg_list is the cutoff for defining large turning event, it can be a number
%used for all groups or it can be an array that has a length of ii_num*1
for i=1:2:length(varargin)
    if strcmp ('row',varargin{i})
        row=varargin{i+1};
    elseif strcmp ('col',varargin{i})
        col=varargin{i+1};
    elseif strcmp ('hbin',varargin{i})
        hbin=varargin{i+1};
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
    end
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
    idx=find(reori_deg>=deg_l);
    xbin2=xmax/2;
    hdseries=[-180:hbin:180]';
    hd_series1=[-180+hbin/2:hbin:180-hbin/2]';

    %% 1-2) bin the data based on the heading direction before turning --:UNIT WILL BE TURNING PER MIN
    if ii==1
        tf_x=cell(ii_num,1); tf_prox=cell(ii_num,1); tf_dist=cell(ii_num,1);
        tf_x_L=cell(ii_num,1); tf_prox_L=cell(ii_num,1); tf_dist_L=cell(ii_num,1);
    end

    t_x1=bin_data_sum_t(dat.orientation,dat.et,hdseries);
    t_x=sum(t_x1,1,'omitnan')';

    c_x=bin_data_count(pre,hdseries);
    tf_x{ii,1}(:,1)=(c_x./t_x).*60;
    clear t_x1 c_x
    %% 1-3) bin the data based on x pos of each turning event (proximate or distant)
    t_x2=bin_data_sum_t2(dat.orientation,dat.x,hdseries,[0:xbin2:xmax]',dat.et);
    c_x2=bin_data_count2(pre,hdseries,turn_x,[0:xbin2:xmax]');
    tf=(c_x2./t_x2).*60;

    tf_prox{ii,1}(:,1)=tf(:,2);
    tf_dist{ii,1}(:,1)=tf(:,1);
    clear tf sem c_x2
    %% 2-1) Get the turning frequency for turning events larger than 60 deg
    c_x_L=bin_data_count(pre(idx),hdseries);
    tf_x_L{ii,1}(:,1)=(c_x_L./t_x).*60;
    clear c_x_L
    %% 2-2) Get the turning frequency for turning events larger than 60 deg based on the x position (prox and distant)
    c_x2_L=bin_data_count2(pre(idx),hdseries,turn_x(idx),[0:xbin2:xmax]');
    tf_L=(c_x2_L./t_x2).*60;

    tf_prox_L{ii,1}(:,1)=tf_L(:,2);
    tf_dist_L{ii,1}(:,1)=tf_L(:,1);
    clear tf_L c_x2_L

    %% 3-1) bin the data based on t0s of each turning event
    t=[0:tmax/3:tmax]';% split it in three time sections
    t_t=bin_data_sum_t2(dat.orientation,dat.et,hdseries,t,dat.et);
    c_t=bin_data_count2(pre,hdseries,t0s,t);
    tf_t=(c_t./t_t).*60;

    c_t_L=bin_data_count2(pre(idx),hdseries,t0s(idx),t);
    tf_t_L=(c_t_L./t_t).*60;
    if ii==1
        tf_t1=cell(ii_num,1); tf_t2=cell(ii_num,1); tf_t3=cell(ii_num,1);
        tf_t1_L=cell(ii_num,1); tf_t2_L=cell(ii_num,1); tf_t3_L=cell(ii_num,1);
    end

    tf_t1{ii,1}(:,1)=tf_t(:,1); tf_t2{ii,1}(:,1)=tf_t(:,2); tf_t3{ii,1}(:,1)=tf_t(:,3);
    tf_t1_L{ii,1}(:,1)=tf_t_L(:,1); tf_t2_L{ii,1}(:,1)=tf_t_L(:,2); tf_t3_L{ii,1}(:,1)=tf_t_L(:,3);

    clear tf_t tf_t_L
    %% 1) & 2) use bootstrapping to get 95% CI
    l=length(hd_series1);
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

        t0s_bt=vertcat(data.t0s{:});
        if isa(t0s_bt,'cell')
            t0s_bt=vertcat(t0s_bt{:});
        end
        idx_bt=find(reori_deg_bt>=deg_l);

        % bin the turning events based on their pre-turning heading
        % directions
        t_bt1=bin_data_sum_t(data.orientation,data.et,hdseries);
        t_bt=sum(t_bt1,1,'omitnan')';
        c_bt=bin_data_count(pre_bt,hdseries);
        tf_bt(:,k)=(c_bt./t_bt).*60;

        % bin the turning events based on their positions and heading
        % directions
        t_bt2=bin_data_sum_t2(data.orientation,data.x,hdseries,[0:xbin2:xmax]',data.et);
        c_bt2=bin_data_count2(pre_bt,hdseries,turn_x_bt,[0:xbin2:xmax]');
        tf_bt2=(c_bt2./t_bt2).*60;
        tf_bt_prox(:,k)=tf_bt2(:,2);tf_bt_dist(:,k)=tf_bt2(:,1);

        % for turning event larger than 60 deg
        c_bt_L=bin_data_count(pre_bt(idx_bt),hdseries);
        tf_bt_L(:,k)=(c_bt_L./t_bt).*60;
        c_bt2_L=bin_data_count2(pre_bt(idx_bt),hdseries,turn_x_bt(idx_bt),[0:xbin2:xmax]');
        tf_bt2_L=(c_bt2_L./t_bt2).*60;
        tf_bt_prox_L(:,k)=tf_bt2_L(:,2);tf_bt_dist_L(:,k)=tf_bt2_L(:,1);

        %for binning the turnning events based on the timing
        t_bt_t=bin_data_sum_t2(dat.orientation,dat.et,hdseries,t,dat.et);
        c_bt_t=bin_data_count2(pre_bt,hdseries,t0s_bt,t);
        c_bt_t_L=bin_data_count2(pre_bt(idx_bt),hdseries,t0s_bt(idx_bt),t);
        tf_bt_t=(c_bt_t./t_bt_t).*60;
        %Compute the turn freq for large turning events
        tf_bt_t_L=(c_bt_t_L./t_bt_t).*60;
        %save them based on the time
        tf_bt_t1(:,k)=tf_bt_t(:,1); tf_bt_t2(:,k)=tf_bt_t(:,2); tf_bt_t3(:,k)=tf_bt_t(:,3);
        tf_bt_t1_L(:,k)=tf_bt_t_L(:,1); tf_bt_t2_L(:,k)=tf_bt_t_L(:,2); tf_bt_t3_L(:,k)=tf_bt_t_L(:,3);

        clear tf_bt2 tf_bt2_L c_bt c_bt_L idx_bt pre_bt turn_x_bt reori_deg_bt t0s_bt t0_bt
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
    tf_diff=tf_bt-tf_x{ii,1}; tf_prox_diff=tf_bt_prox-tf_prox{ii,1}; tf_dist_diff=tf_bt_dist-tf_dist{ii,1};
    tf_diff_L=tf_bt_L-tf_x_L{ii,1}; tf_prox_diff_L=tf_bt_prox_L-tf_prox_L{ii,1}; tf_dist_diff_L=tf_bt_dist_L-tf_dist_L{ii,1};

    tf_diff=sort(tf_diff,2); tf_prox_diff=sort(tf_prox_diff,2);tf_dist_diff=sort(tf_dist_diff,2);
    tf_diff_L=sort(tf_diff_L,2); tf_prox_diff_L=sort(tf_prox_diff_L,2);tf_dist_diff_L=sort(tf_dist_diff_L,2);

    tf_x{ii,1}(:,3)=tf_diff(:,5); tf_x{ii,1}(:,2)=tf_diff(:,195);
    tf_x_L{ii,1}(:,3)=tf_diff_L(:,5); tf_x_L{ii,1}(:,2)=tf_diff_L(:,195);

    tf_dist{ii,1}(:,3)=tf_dist_diff(:,5); tf_dist{ii,1}(:,2)=tf_dist_diff(:,195);
    tf_dist_L{ii,1}(:,3)=tf_dist_diff_L(:,5); tf_dist_L{ii,1}(:,2)=tf_dist_diff_L(:,195);

    tf_prox{ii,1}(:,3)=tf_prox_diff(:,5); tf_prox{ii,1}(:,2)=tf_prox_diff(:,195);
    tf_prox_L{ii,1}(:,3)=tf_prox_diff_L(:,5); tf_prox_L{ii,1}(:,2)=tf_prox_diff_L(:,195);
    %% for different timing

    tf_diff_t1=tf_bt_t1-tf_t1{ii,1}; tf_diff_t2=tf_bt_t2-tf_t2{ii,1}; tf_diff_t3=tf_bt_t3-tf_t3{ii,1};
    tf_diff_t1_L=tf_bt_t1_L-tf_t1_L{ii,1}; tf_diff_t2_L=tf_bt_t2_L-tf_t2_L{ii,1}; tf_diff_t3_L=tf_bt_t3_L-tf_t3_L{ii,1};

    tf_diff_t1=sort(tf_diff_t1,2); tf_diff_t2=sort(tf_diff_t2,2);tf_diff_t3=sort(tf_diff_t3,2);
    tf_diff_t1_L=sort(tf_diff_t1_L,2); tf_diff_t2_L=sort(tf_diff_t2_L,2);tf_diff_t3_L=sort(tf_diff_t3_L,2);

    tf_t1{ii,1}(:,3)=tf_diff_t1(:,5); tf_t1{ii,1}(:,2)=tf_diff_t1(:,195);
    tf_t1_L{ii,1}(:,3)=tf_diff_t1_L(:,5); tf_t1_L{ii,1}(:,2)=tf_diff_t1_L(:,195);

    tf_t2{ii,1}(:,3)=tf_diff_t2(:,5); tf_t2{ii,1}(:,2)=tf_diff_t2(:,195);
    tf_t2_L{ii,1}(:,3)=tf_diff_t2_L(:,5); tf_t2_L{ii,1}(:,2)=tf_diff_t2_L(:,195);

    tf_t3{ii,1}(:,3)=tf_diff_t3(:,5); tf_t3{ii,1}(:,2)=tf_diff_t3(:,195);
    tf_t3_L{ii,1}(:,3)=tf_diff_t3_L(:,5); tf_t3_L{ii,1}(:,2)=tf_diff_t3_L(:,195);

    clear tf_diff* tf_prox_diff* tf_dist_diff*
    %% 4) clear variables
    clearvars -except cond row* col* ii* *bin* *max name color outdir tf* output_name input_cond *series* deg* t

end
%% Get the ctr and exp groups and plot them together
if ~isempty(input_cond)
    idx1=find(contains(input_cond,"ctr"));
    [ctr_tfx,exp_tfx,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(tf_x,idx1,name,color,'input_cond',input_cond);
    [ctr_tf_prox,exp_tf_prox]=split_data_ctr_exp(tf_prox,idx1,name,color,'input_cond',input_cond);
    [ctr_tf_dist,exp_tf_dist]=split_data_ctr_exp(tf_dist,idx1,name,color,'input_cond',input_cond);

    [ctr_tf_prox_L,exp_tf_prox_L]=split_data_ctr_exp(tf_prox_L,idx1,name,color,'input_cond',input_cond);
    [ctr_tf_dist_L,exp_tf_dist_L]=split_data_ctr_exp(tf_dist_L,idx1,name,color,'input_cond',input_cond);
    [ctr_tfx_L,exp_tfx_L]=split_data_ctr_exp(tf_x_L,idx1,name,color,'input_cond',input_cond);

    [ctr_t1,exp_t1]=split_data_ctr_exp(tf_t1,idx1,name,color,'input_cond',input_cond);
    [ctr_t2,exp_t2]=split_data_ctr_exp(tf_t2,idx1,name,color,'input_cond',input_cond);
    [ctr_t3,exp_t3]=split_data_ctr_exp(tf_t3,idx1,name,color,'input_cond',input_cond);
    % for different timing
    [ctr_t1_L,exp_t1_L]=split_data_ctr_exp(tf_t1_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_L,exp_t2_L]=split_data_ctr_exp(tf_t2_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_L,exp_t3_L]=split_data_ctr_exp(tf_t3_L,idx1,name,color,'input_cond',input_cond);
    %split data for bootstrap null hyppothesis testing
    [ctr_tf_all,exp_tf_all]=split_data_ctr_exp(tf_bt_all,idx1,name,color,'input_cond',input_cond);
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

    %% plot the polar plot
    for i=1:length(exp_tfx)
        % get the color for the ctr and exp groups
        legends=[ctr_name(i,1);exp_name(i,1)];
        color1=[ctr_color(i,:);exp_color(i,:)];
        % GET THE DATA
        ctr_tfx_polar=ctr_tfx{i,1};
        exp_tfx_polar=exp_tfx{i,1};
        ctr_tfx_L_polar=ctr_tfx_L{i,1};
        exp_tfx_L_polar=exp_tfx_L{i,1};

        ctr_tf_prox_polar=ctr_tf_prox{i,1};
        exp_tf_prox_polar=exp_tf_prox{i,1};
        ctr_tf_prox_L_polar=ctr_tf_prox_L{i,1};
        exp_tf_prox_L_polar=exp_tf_prox_L{i,1};

        ctr_tf_dist_polar=ctr_tf_dist{i,1};
        exp_tf_dist_polar=exp_tf_dist{i,1};
        ctr_tf_dist_L_polar=ctr_tf_dist_L{i,1};
        exp_tf_dist_L_polar=exp_tf_dist_L{i,1};

        ctr_t1_polar=ctr_t1{i,1};
        exp_t1_polar=exp_t1{i,1};
        ctr_t1_L_polar=ctr_t1_L{i,1};
        exp_t1_L_polar=exp_t1_L{i,1};


        ctr_t2_polar=ctr_t2{i,1};
        exp_t2_polar=exp_t2{i,1};
        ctr_t2_L_polar=ctr_t2_L{i,1};
        exp_t2_L_polar=exp_t2_L{i,1};

        ctr_t3_polar=ctr_t3{i,1};
        exp_t3_polar=exp_t3{i,1};
        ctr_t3_L_polar=ctr_t3_L{i,1};
        exp_t3_L_polar=exp_t3_L{i,1};


        %% Turning frequency for ctr and exp groups
        hbin1=hbin/2;
        hd_series2=[-180+hbin1/2:hbin1:180-hbin1/2]';
        l1=length(hd_series2); l2=length(ctr_tfx_polar);
        em_arry=zeros(l1-l2,1);

        for j=1:width(ctr_tfx_polar)
            data_tfx_polar{1,1}(:,j)=reshape([ctr_tfx_polar(:,j)';em_arry'],1,[])';
            data_tfx_polar{2,1}(:,j)=reshape([em_arry';exp_tfx_polar(:,j)'],1,[])';
            data_tfx_L_polar{1,1}(:,j)=reshape([ctr_tfx_L_polar(:,j)';em_arry'],1,[])';
            data_tfx_L_polar{2,1}(:,j)=reshape([em_arry';exp_tfx_L_polar(:,j)'],1,[])';

            data_tf_prox_polar{1,1}(:,j)=reshape([ctr_tf_prox_polar(:,j)';em_arry'],1,[])';
            data_tf_prox_polar{2,1}(:,j)=reshape([em_arry';exp_tf_prox_polar(:,j)'],1,[])';
            data_tf_prox_L_polar{1,1}(:,j)=reshape([ctr_tf_prox_L_polar(:,j)';em_arry'],1,[])';
            data_tf_prox_L_polar{2,1}(:,j)=reshape([em_arry';exp_tf_prox_L_polar(:,j)'],1,[])';

            data_tf_dist_polar{1,1}(:,j)=reshape([ctr_tf_dist_polar(:,j)';em_arry'],1,[])';
            data_tf_dist_polar{2,1}(:,j)=reshape([em_arry';exp_tf_dist_polar(:,j)'],1,[])';
            data_tf_dist_L_polar{1,1}(:,j)=reshape([ctr_tf_dist_L_polar(:,j)';em_arry'],1,[])';
            data_tf_dist_L_polar{2,1}(:,j)=reshape([em_arry';exp_tf_dist_L_polar(:,j)'],1,[])';


            data_t1_polar{1,1}(:,j)=reshape([ctr_t1_polar(:,j)';em_arry'],1,[])';
            data_t1_polar{2,1}(:,j)=reshape([em_arry';exp_t1_polar(:,j)'],1,[])';
            data_t1_L_polar{1,1}(:,j)=reshape([ctr_t1_L_polar(:,j)';em_arry'],1,[])';
            data_t1_L_polar{2,1}(:,j)=reshape([em_arry';exp_t1_L_polar(:,j)'],1,[])';

            data_t2_polar{1,1}(:,j)=reshape([ctr_t2_polar(:,j)';em_arry'],1,[])';
            data_t2_polar{2,1}(:,j)=reshape([em_arry';exp_t2_polar(:,j)'],1,[])';
            data_t2_L_polar{1,1}(:,j)=reshape([ctr_t2_L_polar(:,j)';em_arry'],1,[])';
            data_t2_L_polar{2,1}(:,j)=reshape([em_arry';exp_t2_L_polar(:,j)'],1,[])';

            data_t3_polar{1,1}(:,j)=reshape([ctr_t3_polar(:,j)';em_arry'],1,[])';
            data_t3_polar{2,1}(:,j)=reshape([em_arry';exp_t3_polar(:,j)'],1,[])';
            data_t3_L_polar{1,1}(:,j)=reshape([ctr_t3_L_polar(:,j)';em_arry'],1,[])';
            data_t3_L_polar{2,1}(:,j)=reshape([em_arry';exp_t3_L_polar(:,j)'],1,[])';
        end
        mlt_subplt(data_tfx_polar,[],100,row_comp,col_comp,i,append("tf_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1),"  vs. ",exp_name(i,1)),"overlapped porlar plot",'bins',hbin1,'max',8,'bt_CI',1);
        mlt_subplt(data_tfx_L_polar,[],101,row_comp,col_comp,i,append("tf larger than rms_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1),"  vs. ",exp_name(i,1)),"overlapped polar plot",'bins',hbin1,'max',2.5,'bt_CI',1);

        mlt_subplt(data_tf_prox_polar,[],102,row_comp,col_comp,i,append("tf_Proximate region_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1),"  vs. ",exp_name(i,1),"_p_r_o_x"),"overlapped polar plot",'bins',hbin1,'max',8,'bt_CI',1);
        mlt_subplt(data_tf_prox_L_polar,[],103,row_comp,col_comp,i,append("tf larger than rms_Proximate region_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1),"  vs. ",exp_name(i,1),"_p_r_o_x_i_m_a_t_e"),"overlapped polar plot",'bins',hbin1,'max',2.5,'bt_CI',1);

        mlt_subplt(data_tf_dist_polar,[],104,row_comp,col_comp,i,append("tf_Distant region_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1),"  vs. ",exp_name(i,1),"_d_i_s_t"),"overlapped polar plot",'bins',hbin1,'max',8,'bt_CI',1);
        mlt_subplt(data_tf_dist_L_polar,[],105,row_comp,col_comp,i,append("tf larger than rms_Distant region_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1),"  vs. ",exp_name(i,1),"_d_i_s_t_a_n_t"),"overlapped polar plot",'bins',hbin1,'max',2.5,'bt_CI',1);

        mlt_subplt(data_t1_polar,[],112,row_comp,col_comp,i,append("tf_0-300s_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1),"  vs. ",exp_name(i,1),"_0_-_3_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',8,'bt_CI',1);
        mlt_subplt(data_t1_L_polar,[],113,row_comp,col_comp,i,append("tf larger than rms_0-300s_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1),"  vs. ",exp_name(i,1),"_0_-_3_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',2.5,'bt_CI',1);

        mlt_subplt(data_t2_polar,[],114,row_comp,col_comp,i,append("tf_300-600s_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_3_0_0_-_6_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',8,'bt_CI',1);
        mlt_subplt(data_t2_L_polar,[],115,row_comp,col_comp,i,append("tf larger than rms_300-600s_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_3_0_0_-_6_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',2.5,'bt_CI',1);

        mlt_subplt(data_t3_polar,[],116,row_comp,col_comp,i,append("tf_600-900s_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_6_0_0_-_9_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',8,'bt_CI',1);
        mlt_subplt(data_t3_L_polar,[],117,row_comp,col_comp,i,append("tf larger than rms _600-900s_porlar plot_ ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_6_0_0_-_9_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',2.5,'bt_CI',1);

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

        N=width(bt_diff);
        group_p=strings(0);
        data_p=[]; data_p_prox=[]; data_p_dist=[];
        data_p_L=[]; data_p_prox_L=[]; data_p_dist_L=[];

        data_p_t1=[]; data_p_t2=[]; data_p_t3=[];
        data_p_t1_L=[]; data_p_t2_L=[]; data_p_t3_L=[];
        for h=1:length(hd_series1)
            group_p=vertcat(group_p,repmat(num2str(hd_series1(h)),N,1));

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

        mlt_subplt(data_p,group_p,106,row_comp,col_comp,i,append("Bt_test_tf_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1)),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_dist,group_p,107,row_comp,col_comp,i,append("Bt_test_tf_distant_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_d_i_s_t_a_n_t"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_prox,group_p,108,row_comp,col_comp,i,append("Bt_test_tf_proximal_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_p_r_o_x_i_m_a_t_e"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);

        mlt_subplt(data_p_L,group_p,109,row_comp,col_comp,i,append("Bt_test_tf larger than rms_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_dist_L,group_p,110,row_comp,col_comp,i,append("Bt_test_tf larger than rms_distant_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_d_i_s_t_a_n_t"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_prox_L,group_p,111,row_comp,col_comp,i,append("Bt_test_tf larger than rms_proximal_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_p_r_o_x_i_m_a_t_e"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

        mlt_subplt(data_p_t1,group_p,118,row_comp,col_comp,i,append("Bt_test_tf_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_t2,group_p,119,row_comp,col_comp,i,append("Bt_test_tf_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);
        mlt_subplt(data_p_t3,group_p,120,row_comp,col_comp,i,append("Bt_test_tf_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-2.5 2.5]);

        mlt_subplt(data_p_t1_L,group_p,121,row_comp,col_comp,i,append("Bt_test_tf larger than rms_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_t2_L,group_p,122,row_comp,col_comp,i,append("Bt_test_tf larger than rms_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);
        mlt_subplt(data_p_t3_L,group_p,123,row_comp,col_comp,i,append("Bt_test_tf larger than rms_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-1.5 1.5]);

    end
end
%% SAVE FIGURES
% need to modify this if you want to save it in different folder name
outdir1=fullfile(outdir,"turn freq","rms turn size");
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;

end