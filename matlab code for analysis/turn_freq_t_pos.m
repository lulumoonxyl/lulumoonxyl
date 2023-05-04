function []=turn_freq_t_pos(tracker_num,genotype,condition,filename,output_name,varargin)
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
tmax=900;% tbin is used for binning the turn event's t0s
xmax=225; deg_list=60;
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
    %% 1-1) Compute the turning frequency-->count of turning
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
    %% get the cutoff for large turning events
    if length(deg_list)==1
        deg_l=deg_list;%later may be use idx for the large turning event
    else
        deg_l=deg_list(ii);
    end
    disp(append("The cutoff for large turning event of the group ", name(ii)," is ",num2str(deg_l)," deg"));

    idx=find(reori_deg>=deg_l);%later may be use idx for the large turning event
    xbin2=xmax/2;
    hdseries=[-180:hbin:180]';
    hd_series1=[-180+hbin/2:hbin:180-hbin/2]';
    %% count the number of turning events
    if ii==1
        tf_t1_prox=cell(ii_num,1);  tf_t2_prox=cell(ii_num,1);  tf_t3_prox=cell(ii_num,1);
        tf_t1_dist=cell(ii_num,1);  tf_t2_dist=cell(ii_num,1);  tf_t3_dist=cell(ii_num,1);

        tf_t1_prox_L=cell(ii_num,1);  tf_t2_prox_L=cell(ii_num,1);  tf_t3_prox_L=cell(ii_num,1);
        tf_t1_dist_L=cell(ii_num,1);  tf_t2_dist_L=cell(ii_num,1);  tf_t3_dist_L=cell(ii_num,1);
    end
    t=[0:tmax/3:tmax]';% split it in three time sections
    t3=bin_data_sum_t3(dat.orientation,dat.x,dat.et,hdseries,[0:xbin2:xmax]',t,dat.et);
    tt=sum(t3,1,'omitnan');
    c=bin_data_count3(pre,hdseries,turn_x,[0:xbin2:xmax]',t0s,t);
    tf=(c./tt).*60;
    tf_t1_prox{ii,1}(:,1)=tf(:,2,1); tf_t2_prox{ii,1}(:,1)=tf(:,2,2); tf_t3_prox{ii,1}(:,1)=tf(:,2,3);
    tf_t1_dist{ii,1}(:,1)=tf(:,1,1); tf_t2_dist{ii,1}(:,1)=tf(:,1,2); tf_t3_dist{ii,1}(:,1)=tf(:,1,3);
    %compute the turn freq for turning event>60
    c_L=bin_data_count3(pre(idx),hdseries,turn_x(idx),[0:xbin2:xmax]',t0s(idx),t);
    tf_L=(c_L./tt).*60;
    tf_t1_prox_L{ii,1}(:,1)=tf_L(:,2,1); tf_t2_prox_L{ii,1}(:,1)=tf_L(:,2,2); tf_t3_prox_L{ii,1}(:,1)=tf_L(:,2,3);
    tf_t1_dist_L{ii,1}(:,1)=tf_L(:,1,1); tf_t2_dist_L{ii,1}(:,1)=tf_L(:,1,2); tf_t3_dist_L{ii,1}(:,1)=tf_L(:,1,3);
    clear tf tf_L
    %% 2) use bootstrapping to get 95% CI
    l=length(hd_series1);
    len=200;
    tf_bt_t1_prox=zeros(l,len); tf_bt_t2_prox=zeros(l,len); tf_bt_t3_prox=zeros(l,len);
    tf_bt_t1_dist=zeros(l,len); tf_bt_t2_dist=zeros(l,len); tf_bt_t3_dist=zeros(l,len);

    tf_bt_t1_prox_L=zeros(l,len); tf_bt_t2_prox_L=zeros(l,len); tf_bt_t3_prox_L=zeros(l,len);
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
        t0s_bt=vertcat(data.t0s{:});
        if isa(t0s_bt,'cell')
            t0s_bt=vertcat(t0s_bt{:});
        end
        idx_bt=find(reori_deg_bt>=deg_l);

        t_bt3=bin_data_sum_t3(data.orientation,data.x,data.et,hdseries,[0:xbin2:xmax]',t,data.et);
        c_bt=bin_data_count3(pre_bt,hdseries,turn_x_bt,[0:xbin2:xmax]',t0s_bt,t);
        c_bt_L=bin_data_count3(pre_bt(idx_bt),hdseries,turn_x_bt(idx_bt),[0:xbin2:xmax]',t0s_bt(idx_bt),t);

        t_bt=sum(t_bt3,1,'omitnan');
        tf_bt=(c_bt./t_bt).*60;
        tf_bt_L=(c_bt_L./t_bt).*60;

        tf_bt_t1_prox(:,k)=tf_bt(:,2,1); tf_bt_t2_prox(:,k)=tf_bt(:,2,2); tf_bt_t3_prox(:,k)=tf_bt(:,2,3);
        tf_bt_t1_dist(:,k)=tf_bt(:,1,1); tf_bt_t2_dist(:,k)=tf_bt(:,1,2); tf_bt_t3_dist(:,k)=tf_bt(:,1,3);

        tf_bt_t1_prox_L(:,k)=tf_bt_L(:,2,1); tf_bt_t2_prox_L(:,k)=tf_bt_L(:,2,2); tf_bt_t3_prox_L(:,k)=tf_bt_L(:,2,3);
        tf_bt_t1_dist_L(:,k)=tf_bt_L(:,1,1); tf_bt_t2_dist_L(:,k)=tf_bt_L(:,1,2); tf_bt_t3_dist_L(:,k)=tf_bt_L(:,1,3);
        clear tf_bt tf_bt_L c_bt c_bt_L idx_bt pre_bt turn_x_bt reori_deg_bt t0s_bt t0_bt
    end
    if ii==1
        tf_bt_t1_prox_all=cell(ii_num,1); tf_bt_t2_prox_all=cell(ii_num,1);tf_bt_t3_prox_all=cell(ii_num,1);
        tf_bt_t1_dist_all=cell(ii_num,1); tf_bt_t2_dist_all=cell(ii_num,1);tf_bt_t3_dist_all=cell(ii_num,1);

        tf_bt_t1_prox_all_L=cell(ii_num,1); tf_bt_t2_prox_all_L=cell(ii_num,1);tf_bt_t3_prox_all_L=cell(ii_num,1);
        tf_bt_t1_dist_all_L=cell(ii_num,1); tf_bt_t2_dist_all_L=cell(ii_num,1);tf_bt_t3_dist_all_L=cell(ii_num,1);
    end
    %save data for boxplot later
    tf_bt_t1_prox_all{ii,1}=tf_bt_t1_prox; tf_bt_t2_prox_all{ii,1}=tf_bt_t2_prox; tf_bt_t3_prox_all{ii,1}=tf_bt_t3_prox;
    tf_bt_t1_dist_all{ii,1}=tf_bt_t1_dist; tf_bt_t2_dist_all{ii,1}=tf_bt_t2_dist; tf_bt_t3_dist_all{ii,1}=tf_bt_t3_dist;

    tf_bt_t1_prox_all_L{ii,1}=tf_bt_t1_prox_L; tf_bt_t2_prox_all_L{ii,1}=tf_bt_t2_prox_L; tf_bt_t3_prox_all_L{ii,1}=tf_bt_t3_prox_L;
    tf_bt_t1_dist_all_L{ii,1}=tf_bt_t1_dist_L; tf_bt_t2_dist_all_L{ii,1}=tf_bt_t2_dist_L; tf_bt_t3_dist_all_L{ii,1}=tf_bt_t3_dist_L;
    %% Use the data to compute the 95% CI
    % for distant and proximate regions
    tf_diff_t1_prox=tf_bt_t1_prox-tf_t1_prox{ii,1}; tf_diff_t2_prox=tf_bt_t2_prox-tf_t2_prox{ii,1}; tf_diff_t3_prox=tf_bt_t3_prox-tf_t3_prox{ii,1};
    tf_diff_t1_dist=tf_bt_t1_dist-tf_t1_dist{ii,1}; tf_diff_t2_dist=tf_bt_t2_dist-tf_t2_dist{ii,1}; tf_diff_t3_dist=tf_bt_t3_dist-tf_t3_dist{ii,1};

    tf_diff_t1_prox_L=tf_bt_t1_prox_L-tf_t1_prox_L{ii,1}; tf_diff_t2_prox_L=tf_bt_t2_prox_L-tf_t2_prox_L{ii,1}; tf_diff_t3_prox_L=tf_bt_t3_prox_L-tf_t3_prox_L{ii,1};
    tf_diff_t1_dist_L=tf_bt_t1_dist_L-tf_t1_dist_L{ii,1}; tf_diff_t2_dist_L=tf_bt_t2_dist_L-tf_t2_dist_L{ii,1}; tf_diff_t3_dist_L=tf_bt_t3_dist_L-tf_t3_dist_L{ii,1};

    tf_diff_t1_prox=sort(tf_diff_t1_prox,2); tf_diff_t2_prox=sort(tf_diff_t2_prox,2); tf_diff_t3_prox=sort(tf_diff_t3_prox,2);
    tf_diff_t1_dist=sort(tf_diff_t1_dist,2); tf_diff_t2_dist=sort(tf_diff_t2_dist,2); tf_diff_t3_dist=sort(tf_diff_t3_dist,2);

    tf_diff_t1_prox_L=sort(tf_diff_t1_prox_L,2); tf_diff_t2_prox_L=sort(tf_diff_t2_prox_L,2); tf_diff_t3_prox_L=sort(tf_diff_t3_prox_L,2);
    tf_diff_t1_dist_L=sort(tf_diff_t1_dist_L,2); tf_diff_t2_dist_L=sort(tf_diff_t2_dist_L,2); tf_diff_t3_dist_L=sort(tf_diff_t3_dist_L,2);

    tf_t1_prox{ii,1}(:,3)=tf_diff_t1_prox(:,5); tf_t1_prox{ii,1}(:,2)=tf_diff_t1_prox(:,195);
    tf_t1_prox_L{ii,1}(:,3)=tf_diff_t1_prox_L(:,5); tf_t1_prox_L{ii,1}(:,2)=tf_diff_t1_prox_L(:,195);

    tf_t2_prox{ii,1}(:,3)=tf_diff_t2_prox(:,5); tf_t2_prox{ii,1}(:,2)=tf_diff_t2_prox(:,195);
    tf_t2_prox_L{ii,1}(:,3)=tf_diff_t2_prox_L(:,5); tf_t2_prox_L{ii,1}(:,2)=tf_diff_t2_prox_L(:,195);

    tf_t3_prox{ii,1}(:,3)=tf_diff_t3_prox(:,5); tf_t3_prox{ii,1}(:,2)=tf_diff_t3_prox(:,195);
    tf_t3_prox_L{ii,1}(:,3)=tf_diff_t3_prox_L(:,5); tf_t3_prox_L{ii,1}(:,2)=tf_diff_t3_prox_L(:,195);


    tf_t1_dist{ii,1}(:,3)=tf_diff_t1_dist(:,5); tf_t1_dist{ii,1}(:,2)=tf_diff_t1_dist(:,195);
    tf_t1_dist_L{ii,1}(:,3)=tf_diff_t1_dist_L(:,5); tf_t1_dist_L{ii,1}(:,2)=tf_diff_t1_dist_L(:,195);

    tf_t2_dist{ii,1}(:,3)=tf_diff_t2_dist(:,5); tf_t2_dist{ii,1}(:,2)=tf_diff_t2_dist(:,195);
    tf_t2_dist_L{ii,1}(:,3)=tf_diff_t2_dist_L(:,5); tf_t2_dist_L{ii,1}(:,2)=tf_diff_t2_dist_L(:,195);

    tf_t3_dist{ii,1}(:,3)=tf_diff_t3_dist(:,5); tf_t3_dist{ii,1}(:,2)=tf_diff_t3_dist(:,195);
    tf_t3_dist_L{ii,1}(:,3)=tf_diff_t3_dist_L(:,5); tf_t3_dist_L{ii,1}(:,2)=tf_diff_t3_dist_L(:,195);

    clear tf_diff*
    %% 3) clear variables
    clearvars -except cond row* col* ii* *bin* *max name color outdir tf* output_name input_cond *series* deg* t

end

if ~isempty(input_cond)
    idx1=find(contains(input_cond,"ctr"));
    [ctr_t1_prox,exp_t1_prox,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(tf_t1_prox,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_prox,exp_t2_prox]=split_data_ctr_exp(tf_t2_prox,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_prox,exp_t3_prox]=split_data_ctr_exp(tf_t3_prox,idx1,name,color,'input_cond',input_cond);

    [ctr_t1_dist,exp_t1_dist]=split_data_ctr_exp(tf_t1_dist,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_dist,exp_t2_dist]=split_data_ctr_exp(tf_t2_dist,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_dist,exp_t3_dist]=split_data_ctr_exp(tf_t3_dist,idx1,name,color,'input_cond',input_cond);

    [ctr_t1_prox_L,exp_t1_prox_L]=split_data_ctr_exp(tf_t1_prox_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_prox_L,exp_t2_prox_L]=split_data_ctr_exp(tf_t2_prox_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_prox_L,exp_t3_prox_L]=split_data_ctr_exp(tf_t3_prox_L,idx1,name,color,'input_cond',input_cond);

    [ctr_t1_dist_L,exp_t1_dist_L]=split_data_ctr_exp(tf_t1_dist_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_dist_L,exp_t2_dist_L]=split_data_ctr_exp(tf_t2_dist_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_dist_L,exp_t3_dist_L]=split_data_ctr_exp(tf_t3_dist_L,idx1,name,color,'input_cond',input_cond);
    
    %split data for bootstrap null hyppothesis testing
    [ctr_t1_prox_all,exp_t1_prox_all]=split_data_ctr_exp(tf_bt_t1_prox_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t1_dist_all,exp_t1_dist_all]=split_data_ctr_exp(tf_bt_t1_dist_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t1_prox_all_L,exp_t1_prox_all_L]=split_data_ctr_exp(tf_bt_t1_prox_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t1_dist_all_L,exp_t1_dist_all_L]=split_data_ctr_exp(tf_bt_t1_dist_all_L,idx1,name,color,'input_cond',input_cond);

    [ctr_t2_prox_all,exp_t2_prox_all]=split_data_ctr_exp(tf_bt_t2_prox_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_dist_all,exp_t2_dist_all]=split_data_ctr_exp(tf_bt_t2_dist_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_prox_all_L,exp_t2_prox_all_L]=split_data_ctr_exp(tf_bt_t2_prox_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t2_dist_all_L,exp_t2_dist_all_L]=split_data_ctr_exp(tf_bt_t2_dist_all_L,idx1,name,color,'input_cond',input_cond);

    [ctr_t3_prox_all,exp_t3_prox_all]=split_data_ctr_exp(tf_bt_t3_prox_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_dist_all,exp_t3_dist_all]=split_data_ctr_exp(tf_bt_t3_dist_all,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_prox_all_L,exp_t3_prox_all_L]=split_data_ctr_exp(tf_bt_t3_prox_all_L,idx1,name,color,'input_cond',input_cond);
    [ctr_t3_dist_all_L,exp_t3_dist_all_L]=split_data_ctr_exp(tf_bt_t3_dist_all_L,idx1,name,color,'input_cond',input_cond);
    %% plot the polar plot
    for i=1:length(exp_t1_prox)
        % get the color for the ctr and exp groups
        legends=[ctr_name(i,1);exp_name(i,1)];
        color1=[ctr_color(i,:);exp_color(i,:)];
        % GET THE DATA
        ctr_t1_prox_polar=ctr_t1_prox{i,1};
        exp_t1_prox_polar=exp_t1_prox{i,1};
        ctr_t1_prox_L_polar=ctr_t1_prox_L{i,1};
        exp_t1_prox_L_polar=exp_t1_prox_L{i,1};

        ctr_t1_dist_polar=ctr_t1_dist{i,1};
        exp_t1_dist_polar=exp_t1_dist{i,1};
        ctr_t1_dist_L_polar=ctr_t1_dist_L{i,1};
        exp_t1_dist_L_polar=exp_t1_dist_L{i,1};

        ctr_t2_prox_polar=ctr_t2_prox{i,1};
        exp_t2_prox_polar=exp_t2_prox{i,1};
        ctr_t2_prox_L_polar=ctr_t2_prox_L{i,1};
        exp_t2_prox_L_polar=exp_t2_prox_L{i,1};

        ctr_t2_dist_polar=ctr_t2_dist{i,1};
        exp_t2_dist_polar=exp_t2_dist{i,1};
        ctr_t2_dist_L_polar=ctr_t2_dist_L{i,1};
        exp_t2_dist_L_polar=exp_t2_dist_L{i,1};

        ctr_t3_prox_polar=ctr_t3_prox{i,1};
        exp_t3_prox_polar=exp_t3_prox{i,1};
        ctr_t3_prox_L_polar=ctr_t3_prox_L{i,1};
        exp_t3_prox_L_polar=exp_t3_prox_L{i,1};

        ctr_t3_dist_polar=ctr_t3_dist{i,1};
        exp_t3_dist_polar=exp_t3_dist{i,1};
        ctr_t3_dist_L_polar=ctr_t3_dist_L{i,1};
        exp_t3_dist_L_polar=exp_t3_dist_L{i,1};

        %% Turning frequency for ctr and exp groups
        hbin1=hbin/2;
        hd_series2=[-180+hbin1/2:hbin1:180-hbin1/2]';
        l1=length(hd_series2); l2=length(ctr_t1_prox_polar);
        em_arry=zeros(l1-l2,1);
        for j=1:width(ctr_t1_prox_polar)
            data_t1_prox_polar{1,1}(:,j)=reshape([ctr_t1_prox_polar(:,j)';em_arry'],1,[])';
            data_t1_prox_polar{2,1}(:,j)=reshape([em_arry';exp_t1_prox_polar(:,j)'],1,[])';
            data_t1_prox_L_polar{1,1}(:,j)=reshape([ctr_t1_prox_L_polar(:,j)';em_arry'],1,[])';
            data_t1_prox_L_polar{2,1}(:,j)=reshape([em_arry';exp_t1_prox_L_polar(:,j)'],1,[])';

            data_t1_dist_polar{1,1}(:,j)=reshape([ctr_t1_dist_polar(:,j)';em_arry'],1,[])';
            data_t1_dist_polar{2,1}(:,j)=reshape([em_arry';exp_t1_dist_polar(:,j)'],1,[])';
            data_t1_dist_L_polar{1,1}(:,j)=reshape([ctr_t1_dist_L_polar(:,j)';em_arry'],1,[])';
            data_t1_dist_L_polar{2,1}(:,j)=reshape([em_arry';exp_t1_dist_L_polar(:,j)'],1,[])';

            data_t2_prox_polar{1,1}(:,j)=reshape([ctr_t2_prox_polar(:,j)';em_arry'],1,[])';
            data_t2_prox_polar{2,1}(:,j)=reshape([em_arry';exp_t2_prox_polar(:,j)'],1,[])';
            data_t2_prox_L_polar{1,1}(:,j)=reshape([ctr_t2_prox_L_polar(:,j)';em_arry'],1,[])';
            data_t2_prox_L_polar{2,1}(:,j)=reshape([em_arry';exp_t2_prox_L_polar(:,j)'],1,[])';

            data_t2_dist_polar{1,1}(:,j)=reshape([ctr_t2_dist_polar(:,j)';em_arry'],1,[])';
            data_t2_dist_polar{2,1}(:,j)=reshape([em_arry';exp_t2_dist_polar(:,j)'],1,[])';
            data_t2_dist_L_polar{1,1}(:,j)=reshape([ctr_t2_dist_L_polar(:,j)';em_arry'],1,[])';
            data_t2_dist_L_polar{2,1}(:,j)=reshape([em_arry';exp_t2_dist_L_polar(:,j)'],1,[])';

            data_t3_prox_polar{1,1}(:,j)=reshape([ctr_t3_prox_polar(:,j)';em_arry'],1,[])';
            data_t3_prox_polar{2,1}(:,j)=reshape([em_arry';exp_t3_prox_polar(:,j)'],1,[])';
            data_t3_prox_L_polar{1,1}(:,j)=reshape([ctr_t3_prox_L_polar(:,j)';em_arry'],1,[])';
            data_t3_prox_L_polar{2,1}(:,j)=reshape([em_arry';exp_t3_prox_L_polar(:,j)'],1,[])';

            data_t3_dist_polar{1,1}(:,j)=reshape([ctr_t3_dist_polar(:,j)';em_arry'],1,[])';
            data_t3_dist_polar{2,1}(:,j)=reshape([em_arry';exp_t3_dist_polar(:,j)'],1,[])';
            data_t3_dist_L_polar{1,1}(:,j)=reshape([ctr_t3_dist_L_polar(:,j)';em_arry'],1,[])';
            data_t3_dist_L_polar{2,1}(:,j)=reshape([em_arry';exp_t3_dist_L_polar(:,j)'],1,[])';
        end

        mlt_subplt(data_t1_prox_polar,[],100,row_comp,col_comp,i,append("tf_Proximate_0-300s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_p_r_o_x_ _0_-_3_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.3,'bt_CI',1);
        mlt_subplt(data_t2_prox_polar,[],101,row_comp,col_comp,i,append("tf_Proximate_300-600s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_p_r_o_x_ _3_0_0_-_6_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.3,'bt_CI',1);
        mlt_subplt(data_t3_prox_polar,[],102,row_comp,col_comp,i,append("tf_Proximate_600-900s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_p_r_o_x_ _6_0_0_-_9_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.3,'bt_CI',1);

        mlt_subplt(data_t1_prox_L_polar,[],103,row_comp,col_comp,i,append("tf larger than rms_Proximate_0-300s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_p_r_o_x_ _0_-_3_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.12,'bt_CI',1);
        mlt_subplt(data_t2_prox_L_polar,[],104,row_comp,col_comp,i,append("tf larger than rms_Proximate_300-600s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_p_r_o_x_ _3_0_0_-_6_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.12,'bt_CI',1);
        mlt_subplt(data_t3_prox_L_polar,[],105,row_comp,col_comp,i,append("tf larger than rms_Proximate_300-900s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_p_r_o_x_ _6_0_0_-_9_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.12,'bt_CI',1);

        mlt_subplt(data_t1_dist_polar,[],106,row_comp,col_comp,i,append("tf_Distant_0-300s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_d_i_s_t_ _0_-_3_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.3,'bt_CI',1);
        mlt_subplt(data_t2_dist_polar,[],107,row_comp,col_comp,i,append("tf_Distant_300-600s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_d_i_s_t_ _3_0_0_-_6_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.3,'bt_CI',1);
        mlt_subplt(data_t3_dist_polar,[],108,row_comp,col_comp,i,append("tf_Distant_600-900s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_d_i_s_t_ _6_0_0_-_9_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.3,'bt_CI',1);

        mlt_subplt(data_t1_dist_L_polar,[],109,row_comp,col_comp,i,append("tf larger than rms_Distant_0-300s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_d_i_s_t_ _0_-_3_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.12,'bt_CI',1);
        mlt_subplt(data_t2_dist_L_polar,[],110,row_comp,col_comp,i,append("tf larger than rms_Distant_300-600s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_d_i_s_t_ _3_0_0_-_6_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.12,'bt_CI',1);
        mlt_subplt(data_t3_dist_L_polar,[],111,row_comp,col_comp,i,append("tf larger than rms_Distant_600-900s_polar plot_",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),"_l_a_r_g_e_,_d_i_s_t_ _6_0_0_-_9_0_0_s"),"overlapped polar plot",'bins',hbin1,'max',0.12,'bt_CI',1);
        %% compute the bootstrap null hypothesis testing
        bt_diff_t1_prox=exp_t1_prox_all{i,1}-ctr_t1_prox_all{i,1};
        bt_diff_t2_prox=exp_t2_prox_all{i,1}-ctr_t2_prox_all{i,1};
        bt_diff_t3_prox=exp_t3_prox_all{i,1}-ctr_t3_prox_all{i,1};

        bt_diff_t1_prox_L=exp_t1_prox_all_L{i,1}-ctr_t1_prox_all_L{i,1};
        bt_diff_t2_prox_L=exp_t2_prox_all_L{i,1}-ctr_t2_prox_all_L{i,1};
        bt_diff_t3_prox_L=exp_t3_prox_all_L{i,1}-ctr_t3_prox_all_L{i,1};

        bt_diff_t1_dist=exp_t1_dist_all{i,1}-ctr_t1_dist_all{i,1};
        bt_diff_t2_dist=exp_t2_dist_all{i,1}-ctr_t2_dist_all{i,1};
        bt_diff_t3_dist=exp_t3_dist_all{i,1}-ctr_t3_dist_all{i,1};

        bt_diff_t1_dist_L=exp_t1_dist_all_L{i,1}-ctr_t1_dist_all_L{i,1};
        bt_diff_t2_dist_L=exp_t2_dist_all_L{i,1}-ctr_t2_dist_all_L{i,1};
        bt_diff_t3_dist_L=exp_t3_dist_all_L{i,1}-ctr_t3_dist_all_L{i,1};

        N=width(bt_diff_t1_prox);
        group_p=strings(0);
        data_p_t1_prox=[]; data_p_t2_prox=[]; data_p_t3_prox=[];
        data_p_t1_dist=[]; data_p_t2_dist=[]; data_p_t3_dist=[];

        data_p_t1_prox_L=[]; data_p_t2_prox_L=[]; data_p_t3_prox_L=[];
        data_p_t1_dist_L=[]; data_p_t2_dist_L=[]; data_p_t3_dist_L=[];

        for h=1:length(hd_series1)
            group_p=vertcat(group_p,repmat(num2str(hd_series1(h)),N,1));

            data_p_t1_prox=vertcat(data_p_t1_prox,bt_diff_t1_prox(h,:)');
            data_p_t2_prox=vertcat(data_p_t2_prox,bt_diff_t2_prox(h,:)');
            data_p_t3_prox=vertcat(data_p_t3_prox,bt_diff_t3_prox(h,:)');

            data_p_t1_prox_L=vertcat(data_p_t1_prox_L,bt_diff_t1_prox_L(h,:)');
            data_p_t2_prox_L=vertcat(data_p_t2_prox_L,bt_diff_t2_prox_L(h,:)');
            data_p_t3_prox_L=vertcat(data_p_t3_prox_L,bt_diff_t3_prox_L(h,:)');

            data_p_t1_dist=vertcat(data_p_t1_dist,bt_diff_t1_dist(h,:)');
            data_p_t2_dist=vertcat(data_p_t2_dist,bt_diff_t2_dist(h,:)');
            data_p_t3_dist=vertcat(data_p_t3_dist,bt_diff_t3_dist(h,:)');

            data_p_t1_dist_L=vertcat(data_p_t1_dist_L,bt_diff_t1_dist_L(h,:)');
            data_p_t2_dist_L=vertcat(data_p_t2_dist_L,bt_diff_t2_dist_L(h,:)');
            data_p_t3_dist_L=vertcat(data_p_t3_dist_L,bt_diff_t3_dist_L(h,:)');
        end
        mlt_subplt(data_p_t1_prox,group_p,112,row_comp,col_comp,i,append("Bt_test_tf_Proximate_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_p_r_o_x_ _0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);
        mlt_subplt(data_p_t2_prox,group_p,113,row_comp,col_comp,i,append("Bt_test_tf_Proximate_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_p_r_o_x_ _3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);
        mlt_subplt(data_p_t3_prox,group_p,114,row_comp,col_comp,i,append("Bt_test_tf_Proximate_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_p_r_o_x_ _6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);

        mlt_subplt(data_p_t1_prox_L,group_p,115,row_comp,col_comp,i,append("Bt_test_tf_larger_than_rms_Proximate_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_p_r_o_x_ _0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);
        mlt_subplt(data_p_t2_prox_L,group_p,116,row_comp,col_comp,i,append("Bt_test_tf_larger_than_rms_Proximate_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_p_r_o_x_ _3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);
        mlt_subplt(data_p_t3_prox_L,group_p,117,row_comp,col_comp,i,append("Bt_test_tf_larger_than_rms_Proximate_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_p_r_o_x_ _6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);

        mlt_subplt(data_p_t1_dist,group_p,118,row_comp,col_comp,i,append("Bt_test_tf_Distant_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_d_i_s_t_ _0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);
        mlt_subplt(data_p_t2_dist,group_p,119,row_comp,col_comp,i,append("Bt_test_tf_Distant_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_d_i_s_t_ _3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);
        mlt_subplt(data_p_t3_dist,group_p,120,row_comp,col_comp,i,append("Bt_test_tf_Distant_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_d_i_s_t_ _6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);

        mlt_subplt(data_p_t1_dist_L,group_p,121,row_comp,col_comp,i,append("Bt_test_tf_larger_than_rms_Distant_0-300s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_d_i_s_t_ _0_-_3_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);
        mlt_subplt(data_p_t2_dist_L,group_p,122,row_comp,col_comp,i,append("Bt_test_tf_larger_than_rms_Distant_300-600s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_d_i_s_t_ _3_0_0_-_6_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);
        mlt_subplt(data_p_t3_dist_L,group_p,123,row_comp,col_comp,i,append("Bt_test_tf_larger_than_rms_Distant_600-900s_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1)," vs. ",exp_name(i,1),"_d_i_s_t_ _6_0_0_-_9_0_0_s"),"boxplot",'yline',0,'ylim',[-0.1 0.1]);

    end
end
%% SAVE FIGURES
outdir1=fullfile(outdir,"turn freq","rms turn size");
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;