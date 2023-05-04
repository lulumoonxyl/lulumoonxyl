function []=accept_reject_turn_q1_3(tracker_num,genotype,condition,filename,output_name,varargin)
%% this fDistantunction will compute the turning freq per min at the prox vs distant region across different duration (0-300s, 300-600s, and 600-900s)

%% get the directory for all the data
w=1;

for i=1:length(genotype)
    for j=1:length(condition)
        cond(w,1)=fullfile("/project/6010970/screen/olfactory_output/JB_JAABA",tracker_num,genotype{i,1},condition{j,1},filename);
        %   cond(w,1)=fullfile("C:\Users\Fei Huang\OneDrive - McGill University\matlab\Analysis based on larvae\JB_JAABA",tracker_num,genotype{i,1},condition{j,1},filename);
        name(w,1)=append(genotype{i,1},"@",condition{j,1});
        w=w+1;
    end
end

%% set the output name
outdir=fullfile("/project/6010970/screen/olfactory_output/Properties_turn_events",output_name);

%outdir=fullfile("C:\Users\Fei Huang\OneDrive - McGill University\matlab\Analysis based on larvae\olfactory_output\Properties_turn_events",output_name);

ii_num=length(cond);
%% set the properties for plots
row=3;col=3;
color=[0.75,0.75,0.75;0.627,0.82,1;0.659,0.659,0.659;0.871,0.416,0.451];

input_cond=[];%it will have a same length as name, telling us whether the input is 'ctr' or 'exp'
row_comp=2; col_comp=2; % these two will be used for the plot to show the comparison of ctr vs exp group, this number should be based on how many comparison you have
tmax=900;% tbin is used for binning the turn event's t0s
xmax=225; deg_list=60;


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
    elseif strcmp(varargin{i},'input_cond')
        input_cond=string(varargin{i+1});
    end
end
input_cond
for ii=1:ii_num
    load(cond(ii));
    %% 1-1) Compute the turning frequency-->count of turning event/time
    % get the data ready for binning
    pre=vertcat(dat.pre_deg{:});
    post=vertcat(dat.post_deg{:});
    turn_x=vertcat(dat.turn_x{:});
    xbin2=xmax/2;

    pre_q=check_qua_turn(pre);
    %post_q=check_qua_turn(post);

    pre_cos_theta=cos(deg2rad(pre));
    post_cos_theta=cos(deg2rad(post));
    %% 1-2) find turning events in q1 and q3, plot the cos for pre and post as a scatter plot

    idx_q1=find(pre_q==1);
    idx_q3=find(pre_q==3);

    mlt_subplt(post_cos_theta(idx_q1),pre_cos_theta(idx_q1),100,row,col,ii,'Cos theta for pre and post hd in quadrant 1_scatter plot','cos(theta) for pre-turning orientation','cos(theta) for post-turning orientation',color,name,'scatter','alpha',0.8,'marker_size',2,'ylim',[-1 1],'xlim',[-1 1]);
    mlt_subplt(post_cos_theta(idx_q3),pre_cos_theta(idx_q3),101,row,col,ii,'Cos theta for pre and post hd in quadrant 3_scatter plot','cos(theta) for pre-turning orientation','cos(theta) for post-turning orientation',color,name,'scatter','alpha',0.8,'marker_size',2,'ylim',[-1 1],'xlim',[-1 1]);
    %     split the tunring event into prox and dist regions
    idx_prox_q1=find(pre_q==1&turn_x>=xbin2&turn_x<=xmax);
    idx_dist_q1=find(pre_q==1&turn_x<xbin2);

    idx_prox_q3=find(pre_q==3&turn_x>=xbin2&turn_x<=xmax);
    idx_dist_q3=find(pre_q==3&turn_x<xbin2);
    % find the turning events that fit the category based on the cos(tehta of the pre or post direction)

    cor_q1=find(pre_cos_theta(idx_q1)<post_cos_theta(idx_q1));
    cor_prox_q1=find(pre_cos_theta(idx_prox_q1)<post_cos_theta(idx_prox_q1));
    cor_dist_q1=find(pre_cos_theta(idx_dist_q1)<post_cos_theta(idx_dist_q1));

    cor_q3=find(pre_cos_theta(idx_q3)<post_cos_theta(idx_q3));
    cor_prox_q3=find(pre_cos_theta(idx_prox_q3)<post_cos_theta(idx_prox_q3));
    cor_dist_q3=find(pre_cos_theta(idx_dist_q3)<post_cos_theta(idx_dist_q3));


    % compute the acceptance rate based on the event number in the previous
    % step
    if ii==1
        %cos_diff is used to plot the violinplot about the difference between cos(post)-cos(pre)
        cos_diff_q1=cell(1,ii_num); cos_diff_prox_q1=cell(1,ii_num); cos_diff_dist_q1=cell(1,ii_num);
        cos_diff_q3=cell(1,ii_num); cos_diff_prox_q3=cell(1,ii_num); cos_diff_dist_q3=cell(1,ii_num);

%         cos_diff_q1_L=cell(1,ii_num); cos_diff_prox_q1_L=cell(1,ii_num); cos_diff_dist_q1_L=cell(1,ii_num);
%         cos_diff_q3_L=cell(1,ii_num); cos_diff_prox_q3_L=cell(1,ii_num); cos_diff_dist_q3_L=cell(1,ii_num);

        ac=zeros(ii_num,2);% ac_L=zeros(ii_num,2);
        ac_sem=zeros(ii_num,2);% ac_sem_L=zeros(ii_num,2);

        ac_prox=zeros(ii_num,2); % ac_prox_L=zeros(ii_num,2);
        ac_sem_prox=zeros(ii_num,2);% ac_sem_prox_L=zeros(ii_num,2);

        ac_dist=zeros(ii_num,2);% ac_dist_L=zeros(ii_num,2);
        ac_sem_dist=zeros(ii_num,2);% ac_sem_dist_L=zeros(ii_num,2);
    end


    q_legends=["Toward odor";"Away from odor"];
    ac(ii,1)=length(cor_q1)/length(idx_q1);
    ac(ii,2)=length(cor_q3)/length(idx_q3);
    ac_sem(ii,1)=1.96*sqrt(ac(ii,1)*(1-ac(ii,1))/length(idx_q1));
    ac_sem(ii,2)=1.96*sqrt(ac(ii,2)*(1-ac(ii,2))/length(idx_q3));

    ac_prox(ii,1)=length(cor_prox_q1)/length(idx_prox_q1);
    ac_prox(ii,2)=length(cor_prox_q3)/length(idx_prox_q3);
    ac_sem_prox(ii,1)=1.96*sqrt(ac_prox(ii,1)*(1-ac_prox(ii,1))/length(idx_prox_q1));
    ac_sem_prox(ii,2)=1.96*sqrt(ac_prox(ii,2)*(1-ac_prox(ii,2))/length(idx_prox_q3));

    ac_dist(ii,1)=length(cor_dist_q1)/length(idx_dist_q1);
    ac_dist(ii,2)=length(cor_dist_q3)/length(idx_dist_q3);
    ac_sem_dist(ii,1)=1.96*sqrt(ac_dist(ii,1)*(1-ac_dist(ii,1))/length(idx_dist_q1));
    ac_sem_dist(ii,2)=1.96*sqrt(ac_dist(ii,2)*(1-ac_dist(ii,2))/length(idx_dist_q3));

    %% save the difference between cos(post) and cos(pre) and used to plot the violinplot

    cos_diff_q1{1,ii}=post_cos_theta(idx_q1)-pre_cos_theta(idx_q1);
    cos_diff_prox_q1{1,ii}=post_cos_theta(idx_prox_q1)-pre_cos_theta(idx_prox_q1);
    cos_diff_dist_q1{1,ii}=post_cos_theta(idx_dist_q1)-pre_cos_theta(idx_dist_q1);

    cos_diff_q3{1,ii}=post_cos_theta(idx_q3)-pre_cos_theta(idx_q3);
    cos_diff_prox_q3{1,ii}=post_cos_theta(idx_prox_q3)-pre_cos_theta(idx_prox_q3);
    cos_diff_dist_q3{1,ii}=post_cos_theta(idx_dist_q3)-pre_cos_theta(idx_dist_q3);
    %% 1-3) find the turning event larger than rms
%     reori_deg=vertcat(dat.reorient_deg_abs{:});
%     if length(deg_list)==1
%         deg_l=deg_list;
%     else
%         deg_l=deg_list(ii);
%     end
%     disp(append("The cutoff for large turning event of the group ", name(ii)," is ",num2str(deg_l)," deg"));
%     idx=find(reori_deg>=deg_l);
% 
%     pre_cos_theta_L=pre_cos_theta(idx);
%     post_cos_theta_L=post_cos_theta(idx);
% 
%     pre_q_L=pre_q(idx);
%     %     post_q_L=post_q(idx);
%     turn_x_L=turn_x(idx);
% 
%     idx_q1_L=find(pre_q_L==1);
%     idx_q3_L=find(pre_q_L==3);
% 
%     idx_prox_q1_L=find(pre_q_L==1&turn_x_L>=xbin2&turn_x_L<=xmax);
%     idx_dist_q1_L=find(pre_q_L==1&turn_x_L<xbin2);
% 
%     idx_prox_q3_L=find(pre_q_L==3&turn_x_L>=xbin2&turn_x_L<=xmax);
%     idx_dist_q3_L=find(pre_q_L==3&turn_x_L<xbin2);
% 
%     mlt_subplt(post_cos_theta_L(idx_q1_L),pre_cos_theta_L(idx_q1_L),102,row,col,ii,'Turn larger than rms_Cos theta for pre and post hd in quadrant 1_scatter plot','cos(theta) for pre-turning orientation','cos(theta) for post-turning orientation',color,name,'scatter','alpha',0.8,'marker_size',2,'ylim',[-1 1],'xlim',[-1 1]);
%     mlt_subplt(post_cos_theta_L(idx_q3_L),pre_cos_theta_L(idx_q3_L),103,row,col,ii,'Turn larger than rms_Cos theta for pre and post hd in quadrant 3_scatter plot','cos(theta) for pre-turning orientation','cos(theta) for post-turning orientation',color,name,'scatter','alpha',0.8,'marker_size',2,'ylim',[-1 1],'xlim',[-1 1]);

    %% find the turning events that fit the category based on the cos(tehta of the pre or post direction)

%     cor_q1_L=find(pre_cos_theta_L(idx_q1_L)<post_cos_theta_L(idx_q1_L));
%     cor_prox_q1_L=find(pre_cos_theta_L(idx_prox_q1_L)<post_cos_theta_L(idx_prox_q1_L));
%     cor_dist_q1_L=find(pre_cos_theta_L(idx_dist_q1_L)<post_cos_theta_L(idx_dist_q1_L));
% 
%     cor_q3_L=find(pre_cos_theta_L(idx_q3_L)<post_cos_theta_L(idx_q3_L));
%     cor_prox_q3_L=find(pre_cos_theta_L(idx_prox_q3_L)<post_cos_theta_L(idx_prox_q3_L));
%     cor_dist_q3_L=find(pre_cos_theta_L(idx_dist_q3_L)<post_cos_theta_L(idx_dist_q3_L));


%     ac_L(ii,1)=length(cor_q1_L)/length(idx_q1_L);
%     ac_L(ii,2)=length(cor_q3_L)/length(idx_q3_L);
%     ac_sem_L(ii,1)=1.96*sqrt(ac_L(ii,1)*(1-ac_L(ii,1))/length(idx_q1_L));
%     ac_sem_L(ii,2)=1.96*sqrt(ac_L(ii,2)*(1-ac_L(ii,2))/length(idx_q3_L));
% 
%     ac_prox_L(ii,1)=length(cor_prox_q1_L)/length(idx_prox_q1_L);
%     ac_prox_L(ii,2)=length(cor_prox_q3_L)/length(idx_prox_q3_L);
%     ac_sem_prox_L(ii,1)=1.96*sqrt(ac_prox_L(ii,1)*(1-ac_prox_L(ii,1))/length(idx_prox_q1_L));
%     ac_sem_prox_L(ii,2)=1.96*sqrt(ac_prox_L(ii,2)*(1-ac_prox_L(ii,2))/length(idx_prox_q3_L));
% 
%     ac_dist_L(ii,1)=length(cor_dist_q1_L)/length(idx_dist_q1_L);
%     ac_dist_L(ii,2)=length(cor_dist_q3_L)/length(idx_dist_q3_L);
%     ac_sem_dist_L(ii,1)=1.96*sqrt(ac_dist_L(ii,1)*(1-ac_dist_L(ii,1))/length(idx_dist_q1_L));
%     ac_sem_dist_L(ii,2)=1.96*sqrt(ac_dist_L(ii,2)*(1-ac_dist_L(ii,2))/length(idx_dist_q3_L));

    %% save the difference between cos(post) and cos(pre) and used to plot the violinplot

%     cos_diff_q1_L{1,ii}=post_cos_theta_L(idx_q1_L)-pre_cos_theta_L(idx_q1_L);
%     cos_diff_prox_q1_L{1,ii}=post_cos_theta_L(idx_prox_q1_L)-pre_cos_theta_L(idx_prox_q1_L);
%     cos_diff_dist_q1_L{1,ii}=post_cos_theta_L(idx_dist_q1_L)-pre_cos_theta_L(idx_dist_q1_L);
% 
%     cos_diff_q3_L{1,ii}=post_cos_theta_L(idx_q3_L)-pre_cos_theta_L(idx_q3_L);
%     cos_diff_prox_q3_L{1,ii}=post_cos_theta_L(idx_prox_q3_L)-pre_cos_theta_L(idx_prox_q3_L);
%     cos_diff_dist_q3_L{1,ii}=post_cos_theta_L(idx_dist_q3_L)-pre_cos_theta_L(idx_dist_q3_L);
    %% 2-1) store data for chi-sqaured test
    if ii==1
        Y_q1=zeros(ii_num,1); % Y_q1_L=zeros(ii_num,1);
        Y_prox_q1=zeros(ii_num,1); % Y_prox_q1_L=zeros(ii_num,1);
        Y_dist_q1=zeros(ii_num,1); % Y_dist_q1_L=zeros(ii_num,1);

        Y_q3=zeros(ii_num,1); % Y_q3_L=zeros(ii_num,1);
        Y_prox_q3=zeros(ii_num,1);%  Y_prox_q3_L=zeros(ii_num,1);
        Y_dist_q3=zeros(ii_num,1); % Y_dist_q3_L=zeros(ii_num,1);


        N_q1=zeros(ii_num,1);%  N_q1_L=zeros(ii_num,1);
        N_prox_q1=zeros(ii_num,1); % N_prox_q1_L=zeros(ii_num,1);
        N_dist_q1=zeros(ii_num,1); % N_dist_q1_L=zeros(ii_num,1);

        N_q3=zeros(ii_num,1); % N_q3_L=zeros(ii_num,1);
        N_prox_q3=zeros(ii_num,1); % N_prox_q3_L=zeros(ii_num,1);
        N_dist_q3=zeros(ii_num,1); % N_dist_q3_L=zeros(ii_num,1);

        group_q1=cell(ii_num,1);  % group_q1_L=cell(ii_num,1);
        group_prox_q1=cell(ii_num,1); %  group_prox_q1_L=cell(ii_num,1);
        group_dist_q1=cell(ii_num,1);  % group_dist_q1_L=cell(ii_num,1);

        group_q3=cell(ii_num,1); % group_q3_L=cell(ii_num,1);
        group_prox_q3=cell(ii_num,1); % group_prox_q3_L=cell(ii_num,1);
        group_dist_q3=cell(ii_num,1); %  group_dist_q3_L=cell(ii_num,1);

    end
    % these variables save the length of acceptance or non acceptance
    % events
    Y_q1(ii,1)=length(cor_q1); N_q1(ii,1)=length(idx_q1)-length(cor_q1);
    Y_prox_q1(ii,1)=length(cor_prox_q1); N_prox_q1(ii,1)=length(idx_prox_q1)-length(cor_prox_q1);
    Y_dist_q1(ii,1)=length(cor_dist_q1); N_dist_q1(ii,1)=length(idx_dist_q1)-length(cor_dist_q1);

    Y_q3(ii,1)=length(cor_q3); N_q3(ii,1)=length(idx_q3)-length(cor_q3);
    Y_prox_q3(ii,1)=length(cor_prox_q3); N_prox_q3(ii,1)=length(idx_prox_q3)-length(cor_prox_q3);
    Y_dist_q3(ii,1)=length(cor_dist_q3); N_dist_q3(ii,1)=length(idx_dist_q3)-length(cor_dist_q3);

%     Y_q1_L(ii,1)=length(cor_q1_L); N_q1_L(ii,1)=length(idx_q1_L)-length(cor_q1_L);
%     Y_prox_q1_L(ii,1)=length(cor_prox_q1_L); N_prox_q1_L(ii,1)=length(idx_prox_q1_L)-length(cor_prox_q1_L);
%     Y_dist_q1_L(ii,1)=length(cor_dist_q1_L); N_dist_q1_L(ii,1)=length(idx_dist_q1_L)-length(cor_dist_q1_L);
% 
%     Y_q3_L(ii,1)=length(cor_q3_L); N_q3_L(ii,1)=length(idx_q3_L)-length(cor_q3_L);
%     Y_prox_q3_L(ii,1)=length(cor_prox_q3_L); N_prox_q3_L(ii,1)=length(idx_prox_q3_L)-length(cor_prox_q3_L);
%     Y_dist_q3_L(ii,1)=length(cor_dist_q3_L); N_dist_q3_L(ii,1)=length(idx_dist_q3_L)-length(cor_dist_q3_L);

    %the first column of the 'group' is the condition, the second column
    %is "Y" or "N"
    group_q1{ii,1}(:,1)=repmat(append(name(ii),"-Toward odor"),length(idx_q1),1);
    group_q1{ii,1}(:,2)=[repmat("Y",Y_q1(ii,1),1);repmat("N",N_q1(ii,1),1)];
    group_prox_q1{ii,1}(:,1)=repmat(append(name(ii),"-Toward odor-prox"),length(idx_prox_q1),1);
    group_prox_q1{ii,1}(:,2)=[repmat("Y",Y_prox_q1(ii,1),1);repmat("N",N_prox_q1(ii,1),1)];
    group_dist_q1{ii,1}(:,1)=repmat(append(name(ii),"-Toward odor-dist"),length(idx_dist_q1),1);
    group_dist_q1{ii,1}(:,2)=[repmat("Y",Y_dist_q1(ii,1),1);repmat("N",N_dist_q1(ii,1),1)];

    group_q3{ii,1}(:,1)=repmat(append(name(ii),"-Away from odor"),length(idx_q3),1);
    group_q3{ii,1}(:,2)=[repmat("Y",Y_q3(ii,1),1);repmat("N",N_q3(ii,1),1)];
    group_prox_q3{ii,1}(:,1)=repmat(append(name(ii),"-Away from  odor-prox"),length(idx_prox_q3),1);
    group_prox_q3{ii,1}(:,2)=[repmat("Y",Y_prox_q3(ii,1),1);repmat("N",N_prox_q3(ii,1),1)];
    group_dist_q3{ii,1}(:,1)=repmat(append(name(ii),"-Away from  odor-dist"),length(idx_dist_q3),1);
    group_dist_q3{ii,1}(:,2)=[repmat("Y",Y_dist_q3(ii,1),1);repmat("N",N_dist_q3(ii,1),1)];


%     group_q1_L{ii,1}(:,1)=repmat(append(name(ii),"-Toward odor"),length(idx_q1_L),1);
%     group_q1_L{ii,1}(:,2)=[repmat("Y",Y_q1_L(ii,1),1);repmat("N",N_q1_L(ii,1),1)];
%     group_prox_q1_L{ii,1}(:,1)=repmat(append(name(ii),"-Toward odor-prox"),length(idx_prox_q1_L),1);
%     group_prox_q1_L{ii,1}(:,2)=[repmat("Y",Y_prox_q1_L(ii,1),1);repmat("N",N_prox_q1_L(ii,1),1)];
%     group_dist_q1_L{ii,1}(:,1)=repmat(append(name(ii),"-Toward odor-dist"),length(idx_dist_q1_L),1);
%     group_dist_q1_L{ii,1}(:,2)=[repmat("Y",Y_dist_q1_L(ii,1),1);repmat("N",N_dist_q1_L(ii,1),1)];
% 
%     group_q3_L{ii,1}(:,1)=repmat(append(name(ii),"-Away from odor"),length(idx_q3_L),1);
%     group_q3_L{ii,1}(:,2)=[repmat("Y",Y_q3_L(ii,1),1);repmat("N",N_q3_L(ii,1),1)];
%     group_prox_q3_L{ii,1}(:,1)=repmat(append(name(ii),"-Away from  odor-prox"),length(idx_prox_q3_L),1);
%     group_prox_q3_L{ii,1}(:,2)=[repmat("Y",Y_prox_q3_L(ii,1),1);repmat("N",N_prox_q3_L(ii,1),1)];
%     group_dist_q3_L{ii,1}(:,1)=repmat(append(name(ii),"-Away from  odor-dist"),length(idx_dist_q3_L),1);
%     group_dist_q3_L{ii,1}(:,2)=[repmat("Y",Y_dist_q3_L(ii,1),1);repmat("N",N_dist_q3_L(ii,1),1)];
% 

    %% 3) clear variables
    clearvars -except cond row* col* ii* *bin* *max name color outdir tf* output_name input_cond *series* deg* t ac* cond_list group* Y* N* q_legends cos_diff*


end

%% 1-4) plot the acceptance rate
data={ac,ac_sem};
mlt_subplt(data,[],104,1,3,1,'',"Groups","Acceptance rate",color,name,"bar",'ii_num',ii_num,'ebar',1,'legends',q_legends,'ylim',[0 1],'bar_color',[1 0.5 0;0.5 0.25 1],'yline',[0.25;0.75]);

data_prox={ac_prox,ac_sem_prox};
mlt_subplt(data_prox,[],104,1,3,2,'',"Groups","Acceptance rate (prox)",color,name,"bar",'ii_num',ii_num,'ebar',1,'legends',q_legends,'ylim',[0 1],'bar_color',[1 0.5 0;0.5 0.25 1],'yline',[0.25;0.75]);

data_dist={ac_dist,ac_sem_dist};
mlt_subplt(data_dist,[],104,1,3,3,'Acceptance rate in q1 and 13',"Groups","Acceptance rate (dist)",color,name,"bar",'ii_num',ii_num,'ebar',1,'legends',q_legends,'ylim',[0 1],'bar_color',[1 0.5 0;0.5 0.25 1],'yline',[0.25;0.75]);
%% 1-5) violin plot of cos(post)-cos(pre)
for i=1:length(name)
    name1{1,i}=name(i,:);
end
figure(107)
[~,~,cos_diff_q1_mean,cos_diff_q1_med]=violin(cos_diff_q1,'xlabel',name1,'facecolor',color);
change_fig_prop(107,'title',"Toward odor",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre)",'name',"Violinplot_difference between cos post and pre_q1","ylim",[-2 2]);

figure(108)
[~,~,cos_diff_prox_q1_mean,cos_diff_prox_q1_med]=violin(cos_diff_prox_q1,'xlabel',name1,'facecolor',color);
change_fig_prop(108,'title',"Toward odor & prox",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre),prox",'name',"Violinplot_difference between cos post and pre_q1_prox","ylim",[-2 2]);

figure(109)
[~,~,cos_diff_dist_q1_mean,cos_diff_dist_q1_med]=violin(cos_diff_dist_q1,'xlabel',name1,'facecolor',color);
change_fig_prop(109,'title',"Toward odor & dist",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre),dist",'name',"Violinplot_difference between cos post and pre_q1_dist","ylim",[-2 2]);

figure(110)
[~,~,cos_diff_q3_mean,cos_diff_q3_med]=violin(cos_diff_q3,'xlabel',name1,'facecolor',color);
change_fig_prop(110,'title',"Away from odor",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre)",'name',"Violinplot_difference between cos post and pre_q3","ylim",[-2 2]);

figure(111)
[~,~,cos_diff_prox_q3_mean,cos_diff_prox_q3_med]=violin(cos_diff_prox_q3,'xlabel',name1,'facecolor',color);
change_fig_prop(111,'title',"Away from odor & prox",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre),prox",'name',"Violinplot_difference between cos post and pre_q3_prox","ylim",[-2 2]);

figure(112)
[~,~,cos_diff_dist_q3_mean,cos_diff_dist_q3_med]=violin(cos_diff_dist_q3,'xlabel',name1,'facecolor',color);
change_fig_prop(112,'title',"Away from odor & dist",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre),dist",'name',"Violinplot_difference between cos post and pre_q3_dist","ylim",[-2 2]);

% figure(113)
% [~,~,cos_diff_q1_L_mean,cos_diff_q1_L_med]=violin(cos_diff_q1_L,'xlabel',name1,'facecolor',color);
% change_fig_prop(113,'title',"Toward odor & larger than rms turn size",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre)",'name',"Violinplot_difference between cos post and pre_q1_larger than rms turn size","ylim",[-2 2]);
% 
% figure(114)
% [~,~,cos_diff_prox_q1_L_mean,cos_diff_prox_q1_L_med]=violin(cos_diff_prox_q1_L,'xlabel',name1,'facecolor',color);
% change_fig_prop(114,'title',"Toward odor & larger than rms turn size & prox ",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre),prox",'name',"Violinplot_difference between cos post and pre_q1_prox_larger than rms turn size","ylim",[-2 2]);
% 
% figure(115)
% [~,~,cos_diff_dist_q1_L_mean,cos_diff_dist_q1_L_med]=violin(cos_diff_dist_q1_L,'xlabel',name1,'facecolor',color);
% change_fig_prop(115,'title',"Toward odor & larger than rms turn size & dist",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre),dist",'name',"Violinplot_difference between cos post and pre_q1_dist_larger than rms turn size","ylim",[-2 2]);
% 
% figure(116)
% [~,~,cos_diff_q3_L_mean,cos_diff_q3_L_med]=violin(cos_diff_q3_L,'xlabel',name1,'facecolor',color);
% change_fig_prop(116,'title',"Away from odor & larger than rms turn size",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre)",'name',"Violinplot_difference between cos post and pre_q3_larger than rms turn size","ylim",[-2 2]);
% 
% figure(117)
% [~,~,cos_diff_prox_q3_L_mean,cos_diff_prox_q3_L_med]=violin(cos_diff_prox_q3_L,'xlabel',name1,'facecolor',color);
% change_fig_prop(117,'title',"Away from odor & larger than rms turn size & prox",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre),prox",'name',"Violinplot_difference between cos post and pre_q3_prox_larger than rms turn size","ylim",[-2 2]);
% 
% figure(118)
% [~,~,cos_diff_dist_q3_L_mean,cos_diff_dist_q3_L_med]=violin(cos_diff_dist_q3_L,'xlabel',name1,'facecolor',color);
% change_fig_prop(118,'title',"Away from odor & larger than rms turn size & dist",'xtitle',"Groups",'ytitle',"cos(post)-cos(pre),dist",'name',"Violinplot_difference between cos post and pre_q3_dist_larger than rms turn size","ylim",[-2 2]);
%% SAVE FIGURES
outdir1=fullfile(outdir,"acceptance rate");
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;
%% 2) chi-sqaure test between the ctr and exp group
if ~isempty(input_cond)
    idx1=find(contains(input_cond,"ctr"));
    [ctr_q1,exp_q1,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(group_q1,idx1,name,color,'input_cond',input_cond);
    [ctr_prox_q1,exp_prox_q1]=split_data_ctr_exp(group_prox_q1,idx1,name,color,'input_cond',input_cond);
    [ctr_dist_q1,exp_dist_q1]=split_data_ctr_exp(group_dist_q1,idx1,name,color,'input_cond',input_cond);

%     [ctr_q1_L,exp_q1_L]=split_data_ctr_exp(group_q1_L,idx1,name,color,'input_cond',input_cond);
%     [ctr_prox_q1_L,exp_prox_q1_L]=split_data_ctr_exp(group_prox_q1_L,idx1,name,color,'input_cond',input_cond);
%     [ctr_dist_q1_L,exp_dist_q1_L]=split_data_ctr_exp(group_dist_q1_L,idx1,name,color,'input_cond',input_cond);

    [ctr_q3,exp_q3]=split_data_ctr_exp(group_q3,idx1,name,color,'input_cond',input_cond);
    [ctr_prox_q3,exp_prox_q3]=split_data_ctr_exp(group_prox_q3,idx1,name,color,'input_cond',input_cond);
    [ctr_dist_q3,exp_dist_q3]=split_data_ctr_exp(group_dist_q3,idx1,name,color,'input_cond',input_cond);

%     [ctr_q3_L,exp_q3_L]=split_data_ctr_exp(group_q3_L,idx1,name,color,'input_cond',input_cond);
%     [ctr_prox_q3_L,exp_prox_q3_L]=split_data_ctr_exp(group_prox_q3_L,idx1,name,color,'input_cond',input_cond);
%     [ctr_dist_q3_L,exp_dist_q3_L]=split_data_ctr_exp(group_dist_q3_L,idx1,name,color,'input_cond',input_cond);

    comp=strings(length(exp_q1),1);

    p_q1=zeros(length(exp_q1),1); % p_q1_L=zeros(length(exp_q1_L),1);
    p_prox_q1=zeros(length(exp_prox_q1),1);  %p_prox_q1_L=zeros(length(exp_prox_q1_L),1);
    p_dist_q1=zeros(length(exp_dist_q1),1); % p_dist_q1_L=zeros(length(exp_dist_q1_L),1);

    p_q3=zeros(length(exp_q3),1); % p_q3_L=zeros(length(exp_q3_L),1);
    p_prox_q3=zeros(length(exp_prox_q3),1); % p_prox_q3_L=zeros(length(exp_prox_q3_L),1);
    p_dist_q3=zeros(length(exp_dist_q3),1);  %p_dist_q3_L=zeros(length(exp_dist_q3_L),1);

    for i =1:length(exp_q1)

        comp(i,1)=append(ctr_name(i,1),"-",exp_name(i,1));
        x=[ctr_q1{i,1}(:,1);exp_q1{i,1}(:,1)];
        y=[ctr_q1{i,1}(:,2);exp_q1{i,1}(:,2)];
        [~,~,p_q1(i,1)]=crosstab(x,y);
        clear x y

        x=[ctr_prox_q1{i,1}(:,1);exp_prox_q1{i,1}(:,1)];
        y=[ctr_prox_q1{i,1}(:,2);exp_prox_q1{i,1}(:,2)];
        [~,~,p_prox_q1(i,1)]=crosstab(x,y);
        clear x y

        x=[ctr_dist_q1{i,1}(:,1);exp_dist_q1{i,1}(:,1)];
        y=[ctr_dist_q1{i,1}(:,2);exp_dist_q1{i,1}(:,2)];
        [~,~,p_dist_q1(i,1)]=crosstab(x,y);
        clear x y

%         x=[ctr_q1_L{i,1}(:,1);exp_q1_L{i,1}(:,1)];
%         y=[ctr_q1_L{i,1}(:,2);exp_q1_L{i,1}(:,2)];
%         [~,~,p_q1_L(i,1)]=crosstab(x,y);
%         clear x y
% 
% 
%         x=[ctr_prox_q1_L{i,1}(:,1);exp_prox_q1_L{i,1}(:,1)];
%         y=[ctr_prox_q1_L{i,1}(:,2);exp_prox_q1_L{i,1}(:,2)];
%         [~,~,p_prox_q1_L(i,1)]=crosstab(x,y);
%         clear x y
% 
%         x=[ctr_dist_q1_L{i,1}(:,1);exp_dist_q1_L{i,1}(:,1)];
%         y=[ctr_dist_q1_L{i,1}(:,2);exp_dist_q1_L{i,1}(:,2)];
%         [~,~,p_dist_q1_L(i,1)]=crosstab(x,y);
%         clear x y

        x=[ctr_q3{i,1}(:,1);exp_q3{i,1}(:,1)];
        y=[ctr_q3{i,1}(:,2);exp_q3{i,1}(:,2)];
        [~,~,p_q3(i,1)]=crosstab(x,y);
        clear x y

        x=[ctr_prox_q3{i,1}(:,1);exp_prox_q3{i,1}(:,1)];
        y=[ctr_prox_q3{i,1}(:,2);exp_prox_q3{i,1}(:,2)];
        [~,~,p_prox_q3(i,1)]=crosstab(x,y);
        clear x y

        x=[ctr_dist_q3{i,1}(:,1);exp_dist_q3{i,1}(:,1)];
        y=[ctr_dist_q3{i,1}(:,2);exp_dist_q3{i,1}(:,2)];
        [~,~,p_dist_q3(i,1)]=crosstab(x,y);
        clear x y

%         x=[ctr_q3_L{i,1}(:,1);exp_q3_L{i,1}(:,1)];
%         y=[ctr_q3_L{i,1}(:,2);exp_q3_L{i,1}(:,2)];
%         [~,~,p_q3_L(i,1)]=crosstab(x,y);
%         clear x y
% 
%         x=[ctr_prox_q3_L{i,1}(:,1);exp_prox_q3_L{i,1}(:,1)];
%         y=[ctr_prox_q3_L{i,1}(:,2);exp_prox_q3_L{i,1}(:,2)];
%         [~,~,p_prox_q3_L(i,1)]=crosstab(x,y);
%         clear x y
% 
%         x=[ctr_dist_q3_L{i,1}(:,1);exp_dist_q3_L{i,1}(:,1)];
%         y=[ctr_dist_q3_L{i,1}(:,2);exp_dist_q3_L{i,1}(:,2)];
%         [~,~,p_dist_q3_L(i,1)]=crosstab(x,y);
%         clear x y

    end
end
%% 2-3) save data
Data=table(name,Y_q1,N_q1,Y_prox_q1,N_prox_q1,Y_dist_q1,N_dist_q1,Y_q3,N_q3,Y_dist_q3,N_prox_q3,Y_dist_q3,N_dist_q3,ac,ac_sem,ac_prox,ac_sem_prox,ac_dist,ac_sem_dist);
filename=fullfile(outdir1,'acceptance_rate_q1_q3.xlsx');
writetable(Data,filename);

if ~isempty(comp)
    p_val=table(comp,p_q1,p_prox_q1,p_dist_q1,p_q3,p_prox_q3,p_dist_q3);
    filename=fullfile(outdir1,'chi_squared_test_p_q1_q3.xlsx');
    writetable(p_val,filename);
end
clear
end