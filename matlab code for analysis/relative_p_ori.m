function []=relative_p_ori(tracker_num,genotype,condition,filename,output_name,varargin)
%condition should be a string array, genotype can also be a string array
%output_name is the output folder name
w=1;
for i=1:length(genotype)
    for j=1:length(condition)
        cond(w,1)=fullfile("/project/6010970/screen/olfactory_output/JB_JAABA",tracker_num,genotype{i,1},condition{j,1},filename);
        name(w,1)=append(genotype{i,1},"@",condition{j,1});
        w=w+1;
    end
end
%% set the output name
outdir=fullfile("/project/6010970/screen/olfactory_output/relative_prob_of_orientation",output_name);
ii_num=length(cond);
%% set the properties for plots
row=3;col=3;
color=[0.75,0.75,0.75;0.627,0.82,1;0.659,0.659,0.659;0.871,0.416,0.451];
hbin=10;% for binning data
input_cond=[];%it will have a same length as name, telling us whether the input is 'ctr' or 'exp'
row_comp=2; col_comp=2;
xmax=215;
% these two will be used for the plot to show the comparison of ctr vs exp group, this number should be based on how many comparison you have
for i=1:2:length(varargin)
    if strcmp ('row',varargin{i})
        row=varargin{i+1};
    elseif strcmp ('col',varargin{i})
        col=varargin{i+1};
    elseif strcmp ('hbin',varargin{i})
        hbin=varargin{i+1};
    elseif strcmp ('xmax',varargin{i})
        xmax=varargin{i+1};
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
        
    end
end

for ii=1:ii_num
    %% load data-->"data_grad.mat"
    % load("C:\Users\feihu\OneDrive - McGill University\matlab\Analysis based on larvae\data.mat");
    load(cond(ii));
     
    %% 1-1) compute the total time spend in each bin of orientation
    hd_series=[-180:hbin:180]';
    hd_series1=[-180+hbin/2:hbin:180-hbin/2]';
    
    if ii==1
        p_all=cell(ii_num,1);
    end
    
    %t0_idx and t1_idx will be used to eliminate the time when turning
    %bin dat.orientation with [-180:bin:180] for each larvae
    %p is the probability of time that each larvae spend in a HD bin;p_mean is the weighted mean from p
    
    % difference betwee p_mean and p_all: p_mean is for each animal,
    % therefore it has sem; p_all use t/all time
    [t,p_all{ii,1},ttl_t,tsum]=bin_data_sum_t(dat.orientation,dat.et,hd_series);
    
    %% 1-2)Plot the relative p of orientation p_all
    t_all=sum(t,'omitnan')'; 
    %p_all{ii,1}=t_all./sum(tsum);
    t_all1{ii,1}=t_all;
    mlt_subplt(t_all,hd_series1,100,row,col,ii,append('Time spent in orientation@bin=',num2str(hbin), ',bar chart'),'Orientation (deg)','Time(s)',color,name(ii),'bar','ii_num',ii_num);
    
    
    %% 2-1) Compute the time spend in half of the region and hd_series bin-->p_mean_half is the mean p of orientation for each obj
    if ii==1
        p_prox=cell(ii_num,1);p_dist=cell(ii_num,1);
    end 
    [t_half,~,p_half,t_hd_half,tsum_half,ttl_t_half]=bin_data_sum_t2(dat.orientation,dat.x,hd_series,[0:xmax/2:xmax]',dat.et);
    %p_half=t_half./sum(t_hd_half,1);
    p_prox{ii,1}(:,1)=p_half(:,2); 
    p_dist{ii,1}(:,1)=p_half(:,1); 
    clear p_half
    %% use bootstrapping to get 95% CI
    l=length(hd_series1);
    len=200;
    p_bt_prox=zeros(l,len);
    p_bt_dist=zeros(l,len);
    p_bt=zeros(l,len);
    rng('default');
    s = rng;
    [~,bt]=bootstrp(len,[],[1:length(dat.AN)]);
    for k=1:len
        data=dat(bt(:,k),:);
        [~,p_bt(:,k)]=bin_data_sum_t(data.orientation,data.et,hd_series);
        [~,~,p_bt_half]=bin_data_sum_t2(data.orientation,data.x,hd_series,[0:xmax/2:xmax]',data.et);
        p_bt_prox(:,k)=p_bt_half(:,2);
        p_bt_dist(:,k)=p_bt_half(:,1);
        clear data p_bt_half
    end 
    %% save the bootstrap data for compute the PI
    if ii==1
        p_bt_prox_all=cell(ii_num,1);
        p_bt_dist_all=cell(ii_num,1);
        p_bt_all=cell(ii_num,1);
    end 
    p_bt_prox_all{ii,1}=p_bt_prox;
    p_bt_dist_all{ii,1}=p_bt_dist;
    p_bt_all{ii,1}=p_bt;
    %% compute the difference between each bootstrap sample and the actual one
    p_prox_diff=p_bt_prox-p_prox{ii,1};
    p_dist_diff=p_bt_dist-p_dist{ii,1};
    p_diff=p_bt-p_all{ii,1};

    p_diff=sort(p_diff,2);
    p_prox_diff=sort(p_prox_diff,2);
    p_dist_diff=sort(p_dist_diff,2);

    p_all{ii,1}(:,3)=p_diff(:,5);
    p_all{ii,1}(:,2)=p_diff(:,195);

    p_prox{ii,1}(:,3)=p_prox_diff(:,5);
    p_prox{ii,1}(:,2)=p_prox_diff(:,195);
    
    p_dist{ii,1}(:,3)=p_dist_diff(:,5);
    p_dist{ii,1}(:,2)=p_dist_diff(:,195);
    clear p_bt p_dist_diff p_prox_diff
    %% clear variables
    clearvars -except row* col* input_cond color name outdir output_name *bin ii* *max t_all* p* t* hd_* cond
end

%% 1-3) plot the relative p spend in each orientation bin, one line for each condition
mlt_subplt(p_all,hd_series1,103,1,3,1,'',"Orientation(deg)","Relative probability",color,name,"line",'ii_num',ii_num,'title','All position','ylim',[0.01 0.08],'xlim',[-180+hbin/2 180-hbin/2],'bt_CI',1);
mlt_subplt(p_prox,hd_series1,103,1,3,2,'',"Orientation(deg)","Relative probability",color,name,"line",'ii_num',ii_num,'title','Proximate','ylim',[0.01 0.08],'xlim',[-180+hbin/2 180-hbin/2],'bt_CI',1);
mlt_subplt(p_dist,hd_series1,103,1,3,3,append('Relative p of orientation@bin=',num2str(hbin)),"Orientation(deg)","Relative probability",color,name,"line",'ii_num',ii_num,'title','Distant','ylim',[0.01 0.06],'xlim',[-180+hbin/2 180-hbin/2],'bt_CI',1);

%% 2) compare the ctr and exp in the same bar chart
if ~isempty(input_cond)
    idx=find(contains(input_cond,"ctr"));
    %% plot the ctr and exp groups in the same subplot plot
    %     [ctr,exp,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(t_all1,idx,name,color,'input_cond',input_cond);
    [ctr_p,exp_p,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(p_all,idx,name,color,'input_cond',input_cond);
    [ctr_prox,exp_prox]=split_data_ctr_exp(p_prox,idx,name,color,'input_cond',input_cond); 
    [ctr_dist,exp_dist]=split_data_ctr_exp(p_dist,idx,name,color,'input_cond',input_cond);
    %split data for bootstrap null hyppothesis testing
    [ctr_p_all,exp_p_all]=split_data_ctr_exp(p_bt_all,idx,name,color,'input_cond',input_cond); 
    [ctr_p_dist_all,exp_p_dist_all]=split_data_ctr_exp(p_bt_dist_all,idx,name,color,'input_cond',input_cond);
    [ctr_p_prox_all,exp_p_prox_all]=split_data_ctr_exp(p_bt_prox_all,idx,name,color,'input_cond',input_cond);
    for i=1:length(exp_p)
        data={ctr_p{i,1};exp_p{i,1}};
        data_prox={ctr_prox{i,1};exp_prox{i,1}};
        data_dist={ctr_dist{i,1};exp_dist{i,1}};
        legends=[ctr_name(i,1);exp_name(i,1)];
        color1=[ctr_color(i,:);exp_color(i,:)];
        mlt_subplt(data,hd_series1,104,row_comp,col_comp,i,append("Probability of orientation bar plot_",num2str(output_name),"_ctr vs exp"),"Orientation(deg)","Relative probability",color1,append(ctr_name(i,1),"-",exp_name(i,1)),"overlapped bar",'legend',legends,'xlim',[-180 180],'ylim',[0.01 0.06],'bt_CI',1);
        mlt_subplt(data_prox,hd_series1,106,row_comp,col_comp,i,append("Probability of orientation bar plot_proximate_",num2str(output_name),"_ctr vs exp"),"Orientation(deg)","Relative probability",color1,append(ctr_name(i,1),"-",exp_name(i,1),"_p_r_o_x"),"overlapped bar",'legend',legends,'xlim',[-180 180],'ylim',[0.01 0.06],'bt_CI',1);
        mlt_subplt(data_dist,hd_series1,107,row_comp,col_comp,i,append("Probability of orientation bar plot_distant_",num2str(output_name),"_ctr vs exp"),"Orientation(deg)","Relative probability",color1,append(ctr_name(i,1),"-",exp_name(i,1),"_d_i_s_t"),"overlapped bar",'legend',legends,'xlim',[-180 180],'ylim',[0.01 0.06],'bt_CI',1);
        
        %% relative probability of orientation, polar plot between ctr and exp
        ctr_p_polar=ctr_p{i,1};
        exp_p_polar=exp_p{i,1};
        ctr_prox_polar=ctr_prox{i,1};
        exp_prox_polar=exp_prox{i,1};
        ctr_dist_polar=ctr_dist{i,1};
        exp_dist_polar=exp_dist{i,1};
        hbin1=hbin/2;
        hd_series2=[-180+hbin1/2:hbin1:180-hbin1/2]';
        l1=length(hd_series2); l2=length(ctr_p_polar);
        em_arry=zeros(l1-l2,1);
        for k=1:width(ctr_p_polar)
            data_p1{1,1}(:,k)=reshape([ctr_p_polar(:,k)';em_arry'],1,[])';
            data_p1{2,1}(:,k)=reshape([em_arry';exp_p_polar(:,k)'],1,[])';
    
            data_prox1{1,1}(:,k)=reshape([ctr_prox_polar(:,k)';em_arry'],1,[])';
            data_prox1{2,1}(:,k)=reshape([em_arry';exp_prox_polar(:,k)'],1,[])';
            data_dist1{1,1}(:,k)=reshape([ctr_dist_polar(:,k)';em_arry'],1,[])';
            data_dist1{2,1}(:,k)=reshape([em_arry';exp_dist_polar(:,k)'],1,[])';
        end
        mlt_subplt(data_p1,[],105,row_comp,col_comp,i,append("Relative probability of orientation_",num2str(output_name),"_ctr vs exp polar plot")," "," ",color1,append(ctr_name(i,1),"-",exp_name(i,1)),"overlapped polar plot",'bins',hbin1,'max',0.05,'bt_CI',1);
        mlt_subplt(data_prox1,[],108,row_comp,col_comp,i,append("Relative probability of orientation_proximate_",num2str(output_name),"_ctr vs exp polar plot")," "," ",color1,append(ctr_name(i,1),"-",exp_name(i,1),"_p_r_o_x"),"overlapped polar plot",'bins',hbin1,'max',0.05,'bt_CI',1);
        mlt_subplt(data_dist1,[],109,row_comp,col_comp,i,append("Relative probability of orientation_distant_",num2str(output_name),"_ctr vs exp polar plot")," "," ",color1,append(ctr_name(i,1),"-",exp_name(i,1),"_d_i_s_t"),"overlapped polar plot",'bins',hbin1,'max',0.05,'bt_CI',1);
        %% start to compute the bootstrap null hypothesis testing
        bt_diff=exp_p_all{i,1}-ctr_p_all{i,1}; 
        bt_diff_prox=exp_p_prox_all{i,1}-ctr_p_prox_all{i,1}; 
        bt_diff_dist=exp_p_dist_all{i,1}-ctr_p_dist_all{i,1}; 

        N=width(bt_diff);
        group_p=strings(0); data_p=[];
        data_p_prox=[];
        data_p_dist=[];
        for h=1:length(hd_series1)
            group_p=vertcat(group_p,repmat(num2str(hd_series1(h)),N,1));
            data_p=vertcat(data_p,bt_diff(h,:)');

            data_p_prox=vertcat(data_p_prox,bt_diff_prox(h,:)');
            data_p_dist=vertcat(data_p_dist,bt_diff_dist(h,:)');
        end 

        mlt_subplt(data_p,group_p,110,row_comp,col_comp,i,append("Rel_p_bootstrap_hypothesis_testing_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1)),"boxplot",'ylim',[-0.02 0.02],'yline',0);
        mlt_subplt(data_p_dist,group_p,111,row_comp,col_comp,i,append("Rel_p_bootstrap_hypothesis_testing_distant_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_d_i_s_t_a_n_t"),"boxplot",'ylim',[-0.02 0.02],'yline',0);
        mlt_subplt(data_p_prox,group_p,112,row_comp,col_comp,i,append("Rel_p_bootstrap_hypothesis_testing_proximal_",num2str(output_name),"_ctr vs exp boxplot"),"Heading direction (deg)","Difference",[],append(ctr_name(i,1),"-",exp_name(i,1),"_p_r_o_x_i_m_a_l"),"boxplot",'ylim',[-0.02 0.02],'yline',0);
      
    end

end
%% SAVE FIGURES
if ~isfolder(outdir)
    mkdir(outdir);
end
save_all_figures(outdir);
close all;
end