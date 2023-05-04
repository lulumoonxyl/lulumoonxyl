function []=turn_freq_ratio_sL(tracker_num,genotype,condition,filename,output_name,varargin)
%% this function will plot a tacked bar plot of the turning freq, showing how many is small turning events and how many is turning events larger than rms turn size
%% 1) all turning freq
% 2) tf based on position 3) tf across time 4) across both time and
% position
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
color=[0.47 0.67 0.19;0.85 0.33 0.1];
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
    % get the turning freq based on the reorientation angle: 1) small vs 2)
    % large
    %% 1-1) Compute the turning frequency-->count of turning event/time--> small and large
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
    et=dat.et;

    if length(deg_list)==1
        deg_l=deg_list;%later may be use idx for the large turning event
    else
        deg_l=deg_list(ii);
    end
    disp(append("The cutoff for large turning event of the group ", name(ii)," is ",num2str(deg_l)," deg"));
    idx=find(reori_deg>=deg_l);%later may be use idx for the large turning event

    xbin2=xmax/2;
    hdseries=[-180:hbin:180]';

    t=0;
    for i=1:length(et)
      
            t=et{i,1}(end)-et{i,1}(1)+t;
    end

    %% compute the turning freq for small and large turning events
    legend=["Turning event (< rms turn size)";"Turning event (>= rms turn size)"];

    %initialize the viariables
    if ii==1
        tf_sl=zeros(ii_num,2); p_sl=cell(1,ii_num*2);
        tf_sl_q3=zeros(ii_num,2,l); p_sl_q3=cell(l,ii_num*2);
        tf_sl_q3_prox=zeros(ii_num,2,l); p_sl_q3_prox=cell(l,ii_num*2);
        tf_sl_q3_dist=zeros(ii_num,2,l); p_sl_q3_dist=cell(l,ii_num*2);
    end
    tf_sl(ii,1)=((length(reori_deg)-length(idx))/t)*60;
    tf_sl(ii,2)=(length(idx)/t)*60;
    % calculate the ratio between large and small turning events
    p_sl{1,2*ii-1}=append(num2str(round((length(reori_deg)-length(idx))/length(reori_deg)*100,2)),"%");
    p_sl{1,2*ii}=append(num2str(round(length(idx)/length(reori_deg)*100,2)),"%");
    clear t 
    %% 1-2) Compute the turning frequency based on heading direction -> small and large
    t_x1=bin_data_sum_t(dat.orientation,dat.et,hdseries);
    t_x=sum(t_x1,1,'omitnan')';

     %simply combine two 45 deg quadrants into larger one
    if l==4
        t_q3=assign_num_4quadrants(t_x);
    elseif l==3
        t_q3=assign_num_3quadrants(t_x);
    end
    %count the number of turning events based on the pre-turning heading
    %direction
     c_x=bin_data_count(pre,hdseries);

    if l==4
        c_q3=assign_num_4quadrants(c_x);
    elseif l==3
        c_q3=assign_num_3quadrants(c_x);
    end

    % finr the number of large turning events
    c_x_L=bin_data_count(pre(idx),hdseries);

    if l==4
        c_q3_L=assign_num_4quadrants(c_x_L);
    elseif l==3
        c_q3_L=assign_num_3quadrants(c_x_L);
    end

    %assign the computed tf into the existing variable

    tf_sl_q3(ii,1,:)=((c_q3-c_q3_L)./t_q3)*60;
    tf_sl_q3(ii,2,:)=(c_q3_L./t_q3)*60;
    % calculate the ratio between large and small turning events

    p_s=round((c_q3-c_q3_L)./c_q3.*100,2);
    p_L=round(c_q3_L./c_q3.*100,2);

    for i=1:l
        p_sl_q3{i,2*ii-1}=append(num2str(p_s(i)),"%");
        p_sl_q3{i,ii*2}=append(num2str(p_L(i)),"%");
    end

    clear c_x_L t_q3 t_x t_x1 c_q3 p_s p_L
    %% 1-3) Compute the turning frequency based on different pos and heading direction --> dist vs prox
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

    c_x2_L=bin_data_count2(pre(idx),hdseries,turn_x(idx),[0:xbin2:xmax]');

    if l==4
        c_q3_L=assign_num_4quadrants(c_x2_L);
    elseif l==3
        c_q3_L=assign_num_3quadrants(c_x2_L);
    end
    
    % assign the tf in the variables
      tf_L=(c_q3_L./t_q3)*60;
      tf_s=((c_q3-c_q3_L)./t_q3)*60;
 
      tf_sl_q3_dist(ii,1,:)=tf_s(:,1);
      tf_sl_q3_dist(ii,2,:)=tf_L(:,1);

      tf_sl_q3_prox(ii,1,:)=tf_s(:,2);
      tf_sl_q3_prox(ii,2,:)=tf_L(:,2);

      % calculate the ratio between large and small turning events
      p_s=round((c_q3-c_q3_L)./c_q3.*100,2);
      p_L=round(c_q3_L./c_q3.*100,2);

      for i=1:l
          p_sl_q3_prox{i,2*ii-1}=append(num2str(p_s(i,2)),"%");
          p_sl_q3_prox{i,ii*2}=append(num2str(p_L(i,2)),"%");

          p_sl_q3_dist{i,2*ii-1}=append(num2str(p_s(i,1)),"%");
          p_sl_q3_dist{i,ii*2}=append(num2str(p_L(i,1)),"%");
      end
      clear tf_L tf_s c_x_L t_q3 t_x t_x1 c_q3 p_s p_L
    %% clear variables
      clearvars -except cond row* col* ii* *bin* *max name color outdir tf* output_name input_cond *series* deg* t l q_legends tf* legend p*

end
% plot the turning freq as a stacked bar

mlt_subplt(tf_sl,[1:ii_num]',100,1,1,1,'Total tf with ratio','Groups','Turning frequency (#/min)',color,name,'stacked bar','legends',legend,'ii_num',ii_num,'ylim',[0 7],'label',p_sl);

if l==4
    for i=1:l
        mlt_subplt(tf_sl_q3(:,:,i),[1:ii_num]',101,2,2,i,'q4_tf with ratio','Groups','Turning frequency (#/min)',color,name,'stacked bar','legends',legend,'ii_num',ii_num,'ylim',[0 7],'title',q_legends(i,1),'label',p_sl_q3(i,:));
        mlt_subplt(tf_sl_q3_dist(:,:,i),[1:ii_num]',102,2,2,i,'q4_distal_tf with ratio','Groups','Turning frequency (#/min,dist)',color,name,'stacked bar','legends',legend,'ii_num',ii_num,'ylim',[0 7],'title',q_legends(i,1),'label',p_sl_q3_dist(i,:));
        mlt_subplt(tf_sl_q3_prox(:,:,i),[1:ii_num]',103,2,2,i,'q4_proximal_tf with ratio','Groups','Turning frequency (#/min,prox)',color,name,'stacked bar','legends',legend,'ii_num',ii_num,'ylim',[0 7],'title',q_legends(i,1),'label',p_sl_q3_prox(i,:));
        
    end
elseif l==3
    for i=1:l

        mlt_subplt(tf_sl_q3(:,:,i),[1:ii_num]',101,1,3,i,'q3_tf with ratio','Groups','Turning frequency (#/min)',color,name,'stacked bar','legends',legend,'ii_num',ii_num,'ylim',[0 7],'title',q_legends(i,1),'label',p_sl_q3(i,:));
        mlt_subplt(tf_sl_q3_dist(:,:,i),[1:ii_num]',102,1,3,i,'q3_distal tf with ratio','Groups','Turning frequency (#/min,dist)',color,name,'stacked bar','legends',legend,'ii_num',ii_num,'ylim',[0 7],'title',q_legends(i,1),'label',p_sl_q3_dist(i,:));
        mlt_subplt(tf_sl_q3_prox(:,:,i),[1:ii_num]',103,1,3,i,'q3_proximal tf with ratio','Groups','Turning frequency (#/min,prox)',color,name,'stacked bar','legends',legend,'ii_num',ii_num,'ylim',[0 7],'title',q_legends(i,1),'label',p_sl_q3_prox(i,:));
      
    end
end

%% SAVE FIGURES
outdir1=fullfile(outdir,"turn freq","rms turn size",append(num2str(l)," quadrants"));
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;
end
