function []=acceptance_rate(tracker_num,genotype,condition,filename,output_name,varargin)
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
hbin=10;% for binning data-->frequency of reorientation angle
xmax=225;deg_list=60;
input_cond=[];%it will have a same length as name, telling us whether the input is 'ctr' or 'exp'

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
    elseif strcmp(varargin{i},'xmax')
        xmax=varargin{i+1};
    elseif strcmp(varargin{i},'deg_l')
        deg_list=varargin{i+1};
    elseif strcmp(varargin{i},'input_cond')
       input_cond=string(varargin{i+1});
    end
end
input_cond
for ii=1:ii_num
    %% load data
    load(cond(ii));

    %% get the data ready for binning
    pre=vertcat(dat.pre_deg{:});
    post=vertcat(dat.post_deg{:});
    turn_x=vertcat(dat.turn_x{:});
    %     turn_y=vertcat(dat.turn_y{:});
    reori_deg=vertcat(dat.reorient_deg_abs{:});
    %     t0s1=vertcat(dat.t0s{:});t0s=vertcat(t0s1{:});
    xbin2=xmax/2;
     %% get the cutoff for large turning events
    if length(deg_list)==1
        deg_l=deg_list;%later may be use idx for the large turning event
    else
        deg_l=deg_list(ii);
    end
    disp(append("The cutoff for large turning event of the group ", name(ii)," is ",num2str(deg_l)," deg"));
    %% 1-1) Check the quadrants of the pre and post turning angle
    % Quadrant1: -45<=deg<45; qudrant 2: 45<=deg<135;quadrant 4:
    % -135<=deg<-45; quadrant 3: -180<=deg<135 && 135<=deg<=180
    %q_post=check_qua_turn(post);
    q_pre=check_qua_turn(pre);
    %% 1-2) Compute the acceptance rate for turning events in quadrants 2 and 4
    if ii==1
        ac=zeros(ii_num,2);ac_sem=zeros(ii_num,2);
        ac_prox=zeros(ii_num,2);ac_sem_prox=zeros(ii_num,2);
        ac_dist=zeros(ii_num,2);ac_sem_dist=zeros(ii_num,2);
    end
    q_ver=find(q_pre==2|q_pre==4);
    idx=find((q_pre==2|q_pre==4)&abs(pre)>abs(post));
    ac(ii,1)=length(idx)/length(q_ver);
    ac_sem(ii,1)=1.96*sqrt(ac(ii,1)*(1-ac(ii,1))/length(q_ver));
    
    % for turning events larger than 60 deg
    q_ver1=find((q_pre==2|q_pre==4)&reori_deg>deg_l);
    idx1=find((q_pre==2|q_pre==4)&abs(pre)>abs(post)&reori_deg>deg_l);
    ac(ii,2)=length(idx1)/length(q_ver1);
    ac_sem(ii,2)=1.96*sqrt(ac(ii,2)*(1-ac(ii,2))/length(q_ver1));

    %% 1-3) Acceptance rate for proximate and distant regions
    %proximate region
    q_ver_prox=find((q_pre==2|q_pre==4)&turn_x<=xmax&turn_x>=xbin2);
    idx_prox=find((q_pre==2|q_pre==4)&turn_x<=xmax&turn_x>=xbin2&abs(pre)>abs(post));
    ac_prox(ii,1)=length(idx_prox)/length(q_ver_prox);
    ac_sem_prox(ii,1)=1.96*sqrt(ac_prox(ii,1)*(1-ac_prox(ii,1))/length(q_ver_prox));

    q_ver_prox1=find((q_pre==2|q_pre==4)&reori_deg>deg_l&turn_x<=xmax&turn_x>=xbin2);
    idx_prox1=find((q_pre==2|q_pre==4)&turn_x<=xmax&turn_x>=xbin2&abs(pre)>abs(post)&reori_deg>deg_l);
    ac_prox(ii,2)=length(idx_prox1)/length(q_ver_prox1);
    ac_sem_prox(ii,2)=1.96*sqrt(ac_prox(ii,2)*(1-ac_prox(ii,2))/length(q_ver_prox1));

    % Distant region
    q_ver_dist=find((q_pre==2|q_pre==4)&turn_x>=0&turn_x<xbin2);
    idx_dist=find((q_pre==2|q_pre==4)&turn_x>=0&turn_x<xbin2&abs(pre)>abs(post));
    ac_dist(ii,1)=length(idx_dist)/length(q_ver_dist);
    ac_sem_dist(ii,1)=1.96*sqrt(ac_dist(ii,1)*(1-ac_dist(ii,1))/length(q_ver_dist));

    q_ver_dist1=find((q_pre==2|q_pre==4)&reori_deg>deg_l&turn_x>=0&turn_x<xbin2);
    idx_dist1=find((q_pre==2|q_pre==4)&turn_x>=0&turn_x<xbin2&abs(pre)>abs(post)&reori_deg>deg_l);
    ac_dist(ii,2)=length(idx_dist1)/length(q_ver_dist1);
    ac_sem_dist(ii,2)=1.96*sqrt(ac_dist(ii,2)*(1-ac_dist(ii,2))/length(q_ver_dist1));

    %% 2-1) save data for chi-squared test
    if ii==1
        Y_all=zeros(ii_num,1); Y_large_turn=zeros(ii_num,1);
        Y_all_prox=zeros(ii_num,1); Y_large_turn_prox=zeros(ii_num,1);
        Y_all_dist=zeros(ii_num,1); Y_large_turn_dist=zeros(ii_num,1);

        N_all=zeros(ii_num,1); N_large_turn=zeros(ii_num,1);
        N_all_prox=zeros(ii_num,1); N_large_turn_prox=zeros(ii_num,1);
        N_all_dist=zeros(ii_num,1); N_large_turn_dist=zeros(ii_num,1);

    end

    Y_all(ii,1)=length(idx); N_all(ii,1)=length(q_ver)-length(idx);
    Y_large_turn(ii,1)=length(idx1); N_large_turn(ii,1)=length(q_ver1)-length(idx1);

    Y_all_prox(ii,1)=length(idx_prox); N_all_prox(ii,1)=length(q_ver_prox)-length(idx_prox);
    Y_large_turn_prox(ii,1)=length(idx_prox1); N_large_turn_prox(ii,1)=length(q_ver_prox1)-length(idx_prox1);

    Y_all_dist(ii,1)=length(idx_dist); N_all_dist(ii,1)=length(q_ver_dist)-length(idx_dist);
    Y_large_turn_dist(ii,1)=length(idx_dist1); N_large_turn_dist(ii,1)=length(q_ver_dist1)-length(idx_dist1);
    %% 2-2) Store the data for chi-squared test
    if ii==1
        %the first column of the group is the condition, the second column
        %is "Y" or "N"
        group=cell(ii_num,1);
        group_prox=cell(ii_num,1);
        group_dist=cell(ii_num,1);

        group1=cell(ii_num,1); %FOR TURNING LARGER THAN 60Deg
        group_prox1=cell(ii_num,1);
        group_dist1=cell(ii_num,1);
    end 
    group{ii,1}(:,1)=repmat(name(ii),length(q_ver),1);
    group{ii,1}(:,2)=[repmat("Y",Y_all(ii,1),1);repmat("N",N_all(ii,1),1)];
    group1{ii,1}(:,1)=repmat(name(ii),length(q_ver1),1);
    group1{ii,1}(:,2)=[repmat("Y",Y_large_turn(ii,1),1);repmat("N",N_large_turn(ii,1),1)];

    group_prox{ii,1}(:,1)=repmat(name(ii),length(q_ver_prox),1);
    group_prox{ii,1}(:,2)=[repmat("Y",Y_all_prox(ii,1),1);repmat("N",N_all_prox(ii,1),1)];
    group_prox1{ii,1}(:,1)=repmat(name(ii),length(q_ver_prox1),1);
    group_prox1{ii,1}(:,2)=[repmat("Y",Y_large_turn_prox(ii,1),1);repmat("N",N_large_turn_prox(ii,1),1)];

    group_dist{ii,1}(:,1)=repmat(name(ii),length(q_ver_dist),1);
    group_dist{ii,1}(:,2)=[repmat("Y",Y_all_dist(ii,1),1);repmat("N",N_all_dist(ii,1),1)];
    group_dist1{ii,1}(:,1)=repmat(name(ii),length(q_ver_dist1),1);
    group_dist1{ii,1}(:,2)=[repmat("Y",Y_large_turn_dist(ii,1),1);repmat("N",N_large_turn_dist(ii,1),1)];
    %% 3) clear variables
    clearvars -except cond row* col* ii* *bin* *max name color outdir output_name input_cond ac* deg* Y* N* group* input_cond

end

%% 1-4) Plot the acceptance rate
data={ac,ac_sem};
mlt_subplt(data,[],100,1,1,1,'Acceptance Rate',"Groups","Acceptance rate",color,name,"bar",'ii_num',ii_num,'ebar',1,'legends',["All turning events";"Turning event > rms turn size"],'ylim',[0 1],'bar_color',[0.8 0.8 0.8;0.47 0.67 0.19],'yline',0.5);

data_prox={ac_prox,ac_sem_prox};
mlt_subplt(data_prox,[],101,1,1,1,'Acceptance Rate (Prox)',"Groups","Acceptance rate (prox)",color,name,"bar",'ii_num',ii_num,'ebar',1,'legends',["All turning events";"Turning event > rms turn size"],'ylim',[0 1],'bar_color',[0.8 0.8 0.8;0.47 0.67 0.19],'yline',0.5);

data_dist={ac_dist,ac_sem_dist};
mlt_subplt(data_dist,[],102,1,1,1,'Acceptance Rate (Dist)',"Groups","Acceptance rate (dist)",color,name,"bar",'ii_num',ii_num,'ebar',1,'legends',["All turning events";"Turning event > rms turn size"],'ylim',[0 1],'bar_color',[0.8 0.8 0.8;0.47 0.67 0.19],'yline',0.5);

%% SAVE FIGURES
outdir1=fullfile(outdir,"acceptance rate");
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;
%% Chi-squared test
if ~isempty(input_cond)
    idx1=find(contains(input_cond,"ctr"));
    [ctr,exp,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(group,idx1,name,color,'input_cond',input_cond);
    [ctr1,exp1]=split_data_ctr_exp(group1,idx1,name,color,'input_cond',input_cond);
    [ctr_prox,exp_prox]=split_data_ctr_exp(group_prox,idx1,name,color,'input_cond',input_cond);
    [ctr_prox1,exp_prox1]=split_data_ctr_exp(group_prox1,idx1,name,color,'input_cond',input_cond);
    [ctr_dist,exp_dist]=split_data_ctr_exp(group_dist,idx1,name,color,'input_cond',input_cond);
    [ctr_dist1,exp_dist1]=split_data_ctr_exp(group_dist1,idx1,name,color,'input_cond',input_cond);

    comp=strings(length(exp),1);

    p_all=zeros(length(exp),1);
    p_prox=zeros(length(exp_prox),1);
    p_dist=zeros(length(exp_dist),1);

    p_large_turn=zeros(length(exp1),1);
    p_large_turn_prox=zeros(length(exp_prox1),1);
    p_large_turn_dist=zeros(length(exp_dist1),1);

    for i =1:length(exp)
        comp(i,1)=append(ctr_name(i,1),"-",exp_name(i,1));
        x=[ctr{i,1}(:,1);exp{i,1}(:,1)];
        y=[ctr{i,1}(:,2);exp{i,1}(:,2)];
        [~,~,p_all(i,1)]=crosstab(x,y);
        clear x y

        x=[ctr1{i,1}(:,1);exp1{i,1}(:,1)];
        y=[ctr1{i,1}(:,2);exp1{i,1}(:,2)];
        [~,~,p_large_turn(i,1)]=crosstab(x,y);
        clear x y

        x=[ctr_prox{i,1}(:,1);exp_prox{i,1}(:,1)];
        y=[ctr_prox{i,1}(:,2);exp_prox{i,1}(:,2)];
        [~,~,p_prox(i,1)]=crosstab(x,y);
        clear x y

        x=[ctr_prox1{i,1}(:,1);exp_prox1{i,1}(:,1)];
        y=[ctr_prox1{i,1}(:,2);exp_prox1{i,1}(:,2)];
        [~,~,p_large_turn_prox(i,1)]=crosstab(x,y);
        clear x y

        x=[ctr_dist{i,1}(:,1);exp_dist{i,1}(:,1)];
        y=[ctr_dist{i,1}(:,2);exp_dist{i,1}(:,2)];
        [~,~,p_dist(i,1)]=crosstab(x,y);
        clear x y

        x=[ctr_dist1{i,1}(:,1);exp_dist1{i,1}(:,1)];
        y=[ctr_dist1{i,1}(:,2);exp_dist1{i,1}(:,2)];
        [~,~,p_large_turn_dist(i,1)]=crosstab(x,y);
        clear x y

        
    end 
end 
%% 2-2) save data
Data=table(name,Y_all,N_all,Y_large_turn,N_large_turn,Y_all_prox,N_all_prox,Y_large_turn_prox,N_large_turn_prox, ...
    Y_all_dist,N_all_dist,Y_large_turn_dist,N_large_turn_dist,ac,ac_sem,ac_prox,ac_sem_prox,ac_dist,ac_sem_dist);

filename=fullfile(outdir1,'acceptance_rate.xlsx');
writetable(Data,filename);
if ~isempty(comp)
    p_val=table(comp,p_all,p_large_turn,p_prox,p_large_turn_prox,p_dist,p_large_turn_dist);
    filename=fullfile(outdir1,'chi_squared_test_p.xlsx');
    writetable(p_val,filename);
end
clear
end