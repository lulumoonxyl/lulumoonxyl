function []=cal_PI(tracker_num,genotype,condition,filename,output_name,varargin)


% load("C:\Users\feihu\OneDrive - McGill University\matlab\Analysis based on larvae\olfactory_output\data.mat")
w=1;
%    name=["10E-1GA","10E-2GA","10E-3GA","10^-^5EA","H2O"];
for i=1:length(genotype)
    for j=1:length(condition)
        cond(w,1)=fullfile("/project/6010970/screen/olfactory_output/choreography",tracker_num,genotype{i,1},condition{j,1},filename);
        name(w,1)=append(genotype{i,1},"@",condition{j,1});
        w=w+1;
    end
end

outdir=fullfile("/project/6010970/screen/olfactory_output/Preferential index",output_name)

%set the properties for plots
row=3;col=3;
color=[0,0,0;0.627,0.82,1;0.439,0.733,1;0.255,0.643,1;0.871,0.416,0.451];
input_cond=[];%it will have a same length as name, telling us whether the input is 'ctr' or 'exp'
%if there is more than one ctr, needs to input like, spelling needs also be
%correct
%['exp1';'ctr1';'exp2';'ctr2'];
vector=[];%vector indicating the position of the strongest stimulus-->vector=[1 0] if the odor is in the right side after group_choreo_data-->better to check the traj if you are not sure
xmax=230;ymax=250;tmax=900;bin=5;tbin=30;traj=0;hp=0;odor_pos=215;t=300; sti="odor";
xbin=3;ybin=3;%this bin is used for the PI: v_x/sp vs x position or vy/sp in y
grad_min=0; grad_d_min=0;%the min value of gradient (opto or odor)
grad_max=3; grad_d_max=0.05;
for i=1:2:length(varargin)
    if strcmp ('row',varargin{i})
        row=varargin{i+1};
    elseif strcmp ('col',varargin{i})
        col=varargin{i+1};
    elseif strcmp ('color',varargin{i})
        color=varargin{i+1};
    elseif strcmp(varargin{i},'xmax')
        xmax=varargin{i+1};
    elseif strcmp(varargin{i},'ymax')
        ymax=varargin{i+1};
    elseif strcmp(varargin{i},'tmax')
        ymax=varargin{i+1};
    elseif strcmp(varargin{i},'tbin')% tbin is used for the plot of PI vs tsum
        tbin=varargin{i+1};
    elseif strcmp(varargin{i},'bin')
        bin=varargin{i+1};
    elseif strcmp(varargin{i},'xbin')
        xbin=varargin{i+1};
    elseif strcmp(varargin{i},'ybin')
        ybin=varargin{i+1};
    elseif strcmp(varargin{i},'plot traj')
        traj=varargin{i+1};
    elseif strcmp(varargin{i},'plot heatmap')
        hp=varargin{i+1};
    elseif strcmp(varargin{i},'odor pos')
        odor_pos=varargin{i+1};
    elseif strcmp(varargin{i},'PI t')
        t=varargin{i+1};
    elseif strcmp(varargin{i},'name')
        clear name
        name=string(varargin{i+1})
    elseif strcmp(varargin{i},'stimulus')
        %either opto or odor
        sti=string(varargin{i+1});
    elseif strcmp(varargin{i},'grad_min')
        grad_min=varargin{i+1};
    elseif strcmp(varargin{i},'grad_max')
        grad_max=varargin{i+1};
    elseif strcmp(varargin{i},'grad_d_min')
        grad_d_min=varargin{i+1};
    elseif strcmp(varargin{i},'grad_d_max')
        grad_d_max=varargin{i+1};
    elseif strcmp(varargin{i},'input_cond')
        input_cond=string(varargin{i+1});
    elseif strcmp(varargin{i},'vector')
        vector=varargin{i+1};
    end
end

xseries=[0:bin:xmax]'; x_series=[0+bin/2:bin:xmax-bin/2]';
yseries=[0:bin:ymax]'; y_series=[0+bin/2:bin:ymax-bin/2]';
tseries=[0:tbin:tmax]';  t_series=[0+tbin/2:tbin:tmax-tbin/2]';


ii_num=length(cond)
for ii=1:ii_num
    load(cond(ii));
    %Note

    % 1) if you want to calculate the preferential index (v_x/speed) based
    % on different interval; change the t value and enter it as an input to the
    % function 'bin_PI_cell'
    % 2) change ii_num,name,color, row,col based on how many conditions you
    % have
    %% check whether all the AN is an integer
    idx=find(round(dat_grouped.AN(:))~=dat_grouped.AN(:));
    if ~isempty(idx)
        warning(append('For data in ',cond(ii),' in row=',num2str(idx),', there are abnormal animal numbers, solutions: run choreography/JB/JAABA again'));
        quit
    end
    %% check the number of animals tracked across time
    t_series2=[0:1:tmax]';
    c_tracked=bin_data_count_tracked(dat_grouped.et,t_series2);
    c_tra(:,ii)=c_tracked;

    %% 1) Plot the tracjectory
    x=vertcat(dat_grouped.x{:});
    y=vertcat(dat_grouped.y{:});

    if traj==1
        mlt_subplt(y,x,100,row,col,ii,'Trajectory','X(mm)','Y(mm)',color,name,'scatter','alpha',0.3);
        %plot the centered traj, all traj starts from [0,0]
        xcen=vertcat(dat_grouped.xcentered{:});
        ycen=vertcat(dat_grouped.ycentered{:});
        mlt_subplt(ycen,xcen,102,row,col,ii,'Centered Trajectory','X(mm)','Y(mm)',color,name,'scatter','alpha',0.3);
    end



    %% 2)plot heatmap (another representation of trajectory)
    xlen=length(xseries)-1;
    ylen=length(yseries)-1;
    tlen=length(tseries)-1;
    if hp==1
        if ii==1
            lim=zeros(ii_num,2);
        end

        num_hp=bin_data_count2(vertcat(dat_grouped.y{:}),yseries,vertcat(dat_grouped.x{:}),xseries);

        lim(ii,:)= mlt_subplt(flipud(num_hp),x_series,103,row,col,ii,'Heatmap of traj','X(mm)','Y(mm)',color,name(ii),'heatmap','series2',y_series);
        if ii==ii_num
            axis_max=max(lim(:,2));
            axis_min=min(lim(:,1));
            for i=1:ii_num
                figure(103)
                subplot(row,col,i)
                caxis([axis_min axis_max]);
            end
        end
    end
    clear x y xcen ycen
    %% 3) calculate the distance from the center of mass to the odor
    for i=1:length(dat_grouped.x)
        dis_to{i,1}=odor_pos-dat_grouped.x{i,1};
    end

    if ii==1
        dis_t=cell(ii_num,1);
    end
    tseries1=0:1:tmax;
    dis_t{ii,1}=bin_data_mean(vertcat(dat_grouped.et{:}),tseries1,vertcat(dis_to{:}));

    %% 4) Calculate the navigational index v_x/sp

    if ii==1
        PI_mean=zeros(ii_num,2);PI_mean_y=zeros(ii_num,2);
        vx=cell(ii_num,1);sp=cell(ii_num,1);PI_tpl=cell(ii_num,1);
        str1=[];PI_larvae=[]; PI_size=[];
        str1_y=[];PI_larvae_y=[]; PI_size_y_all=[];
    end
    %% 4-1)calculate the navigational index:v_x/speed for each timepoint, input

    if strcmp(sti,'opto')
        [dat_grouped.PI_larva,v_x,v_t,spd,PI_tp,PI_tp_x,PI_tp_t,PI_tp_grad,PI_tp_grad_diff]=cal_PI_larva(dat_grouped.x,dat_grouped.y,dat_grouped.et,'grad',dat_grouped.grad,'grad_diff',dat_grouped.grad_diff);
    else
        [dat_grouped.PI_larva,v_x,v_t,spd,PI_tp,PI_tp_x,PI_tp_t]=cal_PI_larva(dat_grouped.x,dat_grouped.y,dat_grouped.et);
        %compute the PI using v_y/v
        [dat_grouped.PI_larva_y,~,~,~,PI_tp_y1,PI_tp_y,PI_tp_t_y]=cal_PI_larva(dat_grouped.y,dat_grouped.x,dat_grouped.et);
    end

    PI_mean(ii,1)=mean(dat_grouped.PI_larva);
    PI_mean(ii,2)=std(dat_grouped.PI_larva)/sqrt(length(dat_grouped.PI_larva));
    PI_mean_y(ii,1)=mean(dat_grouped.PI_larva_y);
    PI_mean_y(ii,2)=std(dat_grouped.PI_larva_y)/sqrt(length(dat_grouped.PI_larva_y));
    %% bin the v_x and spd based on time bin
    v_x=vertcat(v_x{:});
    v_t=vertcat(v_t{:});
    spd=vertcat(spd{:});
    PI_tp1=vertcat(PI_tp{:});

    vx{ii,1}=bin_data_mean(v_t,tseries1,v_x);
    sp{ii,1}=bin_data_mean(v_t,tseries1,spd);
    PI_tpl{ii,1}=bin_data_mean(v_t,[0:1:tmax],PI_tp1);
    %% 4-2) compute the PI based on different duration bin (t=300s is the default), and region (xbin=3 is the default)
    if ~isempty(t)&& ~isempty(xbin)
        time=0:t:tmax;lt=length(time)-1;
        xpos=0:xmax/xbin:xmax;lx=length(xpos)-1;
        ypos=0:ymax/ybin:ymax;ly=length(ypos)-1;
        if ii==1
            PI_t_mean=zeros(ii_num,lt);PI_t_sem=zeros(ii_num,lt);
            PI_t_mean_y=zeros(ii_num,lt);PI_t_sem_y=zeros(ii_num,lt);
            PI_x_mean=zeros(ii_num,lx);PI_x_sem=zeros(ii_num,lx);
            PI_y_mean=zeros(ii_num,ly);PI_y_sem=zeros(ii_num,ly);
            PI_x=cell(ii_num,1); PI_y=cell(ii_num,1);
        end


       [PI_t_mean(ii,:),PI_t_sem(ii,:),PI_larva_t,PI_size_t,group_t]=bin_PI_cell(PI_tp_t,time,PI_tp,name(ii));
       [PI_t_mean_y(ii,:),PI_t_sem_y(ii,:),PI_larva_t_y,PI_size_t_y,group_t_y]=bin_PI_cell(PI_tp_t_y,time,PI_tp_y1,name(ii));
        % if bin is too huge, change another way to plot-->plot line
        % instead of bar

        if xbin<=6
            [PI_x_mean(ii,:),PI_x_sem(ii,:),PI_larva_x,PI_size_x,group_x]=bin_PI_cell(PI_tp_x,xpos,PI_tp,name(ii));
            [PI_y_mean(ii,:),PI_y_sem(ii,:),PI_larva_y,PI_size_y,group_y]=bin_PI_cell(PI_tp_y,ypos,PI_tp_y1,name(ii));
        else
            [a,b,PI_larva_x,PI_size_x,group_x]=bin_PI_cell(PI_tp_x,xpos,PI_tp,name(ii));
            PI_x{ii,1}(:,1)=a';PI_x{ii,1}(:,2)=b';
            clear a b
            [a,b,PI_larva_y,PI_size_y,group_y]=bin_PI_cell(PI_tp_y,ypos,PI_tp_y1,name(ii));
            PI_y{ii,1}(:,1)=a';PI_y{ii,1}(:,2)=b';
            clear a b
        end
        %% 4-3) bin the PI based on both duration and postion
        [PI_mean_xt,PI_sem_xt]=bin_PI_cell2(PI_tp_x,xpos,PI_tp_t,[0:5:tmax]',PI_tp);
        [PI_mean_xy,PI_sem_xy]=bin_PI_cell2(PI_tp_y,ypos,PI_tp_t,[0:5:tmax]',PI_tp_y1);
    end
    %     if strcmp(sti,'opto')
    %% 4-4) compute the preferential index for different gradient bin
    %         if ii==1
    %             PI_grad=cell(ii_num,1); PI_grad_diff=cell(ii_num,1);
    %         end
    %         grad_series=[grad_min:0.1:grad_max]';
    %         [c,d,PI_larva_grad,PI_size_grad,group_grad]=bin_PI_cell(PI_tp_grad,grad_series,PI_tp,name(ii));
    %         PI_grad{ii,1}(:,1)=c';PI_grad{ii,1}(:,2)=d';
    %
    %         % save data for stat test
    %         for i=1:length(PI_size_grad)
    %             PI_larva_grad_all{i,1}{ii,1}=PI_larva_grad{i,1};
    %             group_grad_all{i,1}{ii,1}=group_grad{i,1};
    %             size_grad{i,1}{ii,1}=PI_size_grad(i,1);
    %         end
    %         clear c d grad grad_series idx
    %% 4-5) compute the PI for different gradient diff bin
    %         grad_d_series=[grad_d_min:0.001:grad_d_max]';
    %         [c,d,PI_larva_grad_d,PI_size_grad_d,group_grad_d]=bin_PI_cell(PI_tp_grad_diff,grad_d_series,PI_tp,name(ii));
    %         PI_grad_diff{ii,1}(:,1)=c';PI_grad_diff{ii,1}(:,2)=d';
    %         % save data for stat test
    %         for i=1:length(PI_size_grad_d)
    %             PI_larva_grad_d_all{i,1}{ii,1}=PI_larva_grad_d{i,1};
    %             group_grad_d_all{i,1}{ii,1}=group_grad_d{i,1};
    %             size_grad_d{i,1}{ii,1}=PI_size_grad_d(i,1);
    %         end
    %         clear c d grad grad_series idx grad_d_series
    %     end
    %% 4-5) save the data for stat test
    for i=1:lt
        PI_larva_t_all{i,1}{ii,1}=PI_larva_t{i,1};
        group_t_all{i,1}{ii,1}=group_t{i,1};
        size_t{i,1}{ii,1}=PI_size_t(i,1);
    end

    for i=1:lx
        PI_larva_x_all{i,1}{ii,1}=PI_larva_x{i,1};
        group_x_all{i,1}{ii,1}=group_x{i,1};
        size_x{i,1}{ii,1}=PI_size_x(i,1);
    end

    for i=1:ly
        PI_larva_y_all{i,1}{ii,1}=PI_larva_y{i,1};
        group_y_all{i,1}{ii,1}=group_y{i,1};
        size_y{i,1}{ii,1}=PI_size_y(i,1);
    end
    %% save the PI for each region across time for each condition
    %     for i=1:xbin
    %         PI_xt{i,1}{ii,1}(:,1)=PI_mean_xt(i,:)';
    %         PI_xt{i,1}{ii,1}(:,2)=PI_sem_xt(i,:)';
    %     end
    %% save data for stat test for PI
    PI_larvae_y=vertcat(PI_larvae_y,dat_grouped.PI_larva_y);
    PI_size_y_all(ii,1)=length(dat_grouped.PI_larva_y);
    str1_y=vertcat(str1_y,repmat({name(ii)},PI_size_y_all(ii,1),1));

    PI_larvae=vertcat(PI_larvae,dat_grouped.PI_larva);
    PI_size(ii,1)=length(dat_grouped.PI_larva);
    str1=vertcat( str1,repmat({name(ii)},PI_size(ii,1),1));
    %% 5) Caculate the timeseries for percentage distribution ((number of larvae close to the odor-number of larvae away from the odor)/total number of larvare)
    if ii==1
        p_timeseries=cell(ii_num,1);PI_timeseries=cell(ii_num,1);
        p_timeseries_y=cell(ii_num,1);PI_timeseries_y=cell(ii_num,1);
    end

    timebins=5;time=0:timebins:tmax;
    time1=[0+timebins/2:timebins:tmax-timebins/2]';

    edge=[0 xmax/2-xmax*0.1 xmax/2+xmax*0.1 xmax];%separate x pos into three large bins and count the animal in that bins every 5s (timebins)
    edge_y=[0 ymax/2-ymax*0.1 ymax/2+ymax*0.1 ymax];
    [p_timeseries{ii,1},PI_timeseries{ii,1}]=cal_navindex_vs_t(dat_grouped.x,dat_grouped.et,time,edge);
    [p_timeseries_y{ii,1},PI_timeseries_y{ii,1}]=cal_navindex_vs_t(dat_grouped.y,dat_grouped.et,time,edge_y);

    color3=[0.5 1 0;0.5 0 1;1 0 0];
    mlt_subplt(p_timeseries{ii,1},time1,105,row,col,ii,'Frequency of animal distribution across time','Time(s)','Frequency',color3,name(ii),'multiple lines','legends',["Distant";"Central";"Proximate"],'ii_num',ii_num,'lw',2,'ylim',[0 1]);
    mlt_subplt(p_timeseries_y{ii,1},time1,101,row,col,ii,'Frequency of animal distribution across time in y region','Time(s)','Frequency',color3,name(ii),'multiple lines','legends',["Lower";"Central";"Upper"],'ii_num',ii_num,'lw',2,'ylim',[0 1]);

    %% plot PI vs its own tracking time (see whether tsum affect the PI)
    %     tt=0:tbin:tmax;
    %     tt1=tt(1)+tbin/2:tbin:tt(end)-tbin/2;
    %     l=length(tt)-1;
    %     for i=1:l
    %         if i==l
    %             idx=find(dat_grouped.tsum>=tt(i)&dat_grouped.tsum<=tt(i+1));
    %         else
    %             idx=find(dat_grouped.tsum>=tt(i)&dat_grouped.tsum<=tt(i+1));
    %         end
    %         tt_f(i,1)=length(idx)/length(dat_grouped.tsum);
    %
    %     end
    %     mlt_subplt(dat_grouped.tsum,dat_grouped.PI_larva,106,row,col,ii,"PI vs tracking time","Tracking time(s)","PI: v_x/s",color,name(ii),"fitted line");
    %     mlt_subplt(tt_f,tt1,107,row,col,ii,"Frequency of Tracking Interval","Tracking duration(s)","Frequency",color,name(ii),"bar",'ii_num',ii_num);
    %% 6. compute the PI using thetha instead of vx/sp
    if ~isempty(vector)
        for i=1:length(dat_grouped.x)
            x1=dat_grouped.x{i,1}(2:end);
            x2=dat_grouped.x{i,1}(1:end-1);

            y1=dat_grouped.y{i,1}(2:end);
            y2=dat_grouped.y{i,1}(1:end-1);
            deg=cal_ori_deg(x1,x2,y1,y2,vector);
            theta=deg2rad(deg);
            PI_cos{i,1}=cos(theta);
            PI_theta_cos(i,1)=mean(PI_cos{i,1});

            PI_sin{i,1}=sin(theta);
            PI_theta_sin(i,1)=mean(PI_sin{i,1});
            clear x1 x2 y1 y2
        end

        if ii==1
            PI_mean_theta_cos=zeros(ii_num,2);
            PI_mean_theta_sin=zeros(ii_num,2);
        end

        PI_mean_theta_cos(ii,1)=mean(PI_theta_cos);
        PI_mean_theta_cos(ii,2)=std(PI_theta_cos)/length(PI_theta_cos);

        PI_mean_theta_sin(ii,1)=mean(PI_theta_sin);
        PI_mean_theta_sin(ii,2)=std(PI_theta_sin)/length(PI_theta_sin);
        clear PI PI_theta_cos PI_theta_sin PI_sin PI_cos
    end
    %% keep certain variables
    clearvars -except PI* ii* lim*  str1* *series outdir str color* name t time* *bin data dis_* cond *max hp traj t odor_pos row col vx sp ceof grad* c_tra sti size* group* input_cond vector
end
%% plot the number of larvae tracked across time
mlt_subplt(c_tra,[1/2:1:tmax-1/2]',108,1,1,1,'Number of objects tracked across time','Time(s)','Number of objects tracked',color,'','multiple lines','legends',name);

%% plot the navigational index
%set the legends for duration and position
for i=1:xbin
    xleg(1,i)=append("Region ", num2str(i));
end
for i=1:ybin
    yleg(1,i)=append("Region ", num2str(i));
end
dur=[0:t:tmax]';
for i=1:length(dur)-1
    legends(i,1)=append(num2str(dur(i)),"-",num2str(dur(i+1)),"s");
end

%% total preferential index for different conditions
mlt_subplt(PI_mean,[],109,1,1,1,"Preferential Index","Groups","PI:v_x/sp",color,name,"bar",'ii_num',ii_num,'ebar',1,'ylim',[-0.25 0.25]);
mlt_subplt(PI_mean_y,[],180,1,1,1,"Preferential Indexof vy","Groups","PI:v_y/sp",color,name,"bar",'ii_num',ii_num,'ebar',1,'ylim',[-0.25 0.25]);
if ~isempty(vector)
    mlt_subplt(PI_mean_theta_cos,[],170,1,1,1,"Preferential Index,using cos theta","Groups","PI:cos(theta)",color,name,"bar",'ii_num',ii_num,'ebar',1,'ylim',[-0.25 0.25]);
    mlt_subplt(PI_mean_theta_sin,[],171,1,1,1,"Preferential Index,using sine theta","Groups","PI:sin(theta)",color,name,"bar",'ii_num',ii_num,'ebar',1,'ylim',[-0.25 0.25]);

end
data={PI_t_mean,PI_t_sem};
data_t_y={PI_t_mean_y,PI_t_sem_y};
%% preferential index based on time
mlt_subplt(data,[],110,1,1,1,append("Preferential Index vx  for different duration_",num2str(t),"s"),"Groups","PI:v_x/sp",color,name,"bar",'ii_num',ii_num,'title_lab','Different duration','ebar',1,'legends',legends,'ylim',[-0.2 0.3]);
mlt_subplt(data_t_y,[],201,1,1,1,append("Preferential Index vy  for different duration_",num2str(t),"s"),"Groups","PI:v_y/sp",color,name,"bar",'ii_num',ii_num,'title_lab','Different duration','ebar',1,'legends',legends,'ylim',[-0.2 0.3]);

data_x={PI_x_mean,PI_x_sem};
data_y={PI_y_mean,PI_y_sem};

if xbin<=6 %plot the preferential index based on regions
    mlt_subplt(data_x,[],111,1,1,1,append("Preferential Index for different x regions_",num2str(xbin)," regions"),"Groups","PI:v_x/sp",color,name,"bar",'ii_num',ii_num,'title_lab','Different x positions','ebar',1,'legends',xleg,'ylim',[-0.2 0.3]);
    mlt_subplt(data_y,[],200,1,1,1,append("Preferential Index for different y regions_",num2str(xbin)," regions"),"Groups","PI:v_y/sp",color,name,"bar",'ii_num',ii_num,'title_lab','Different y positions','ebar',1,'legends',yleg,'ylim',[-0.2 0.3]);

else
    mlt_subplt(PI_x,[1:xbin]',111,1,1,1,append("Preferential Index for different x regions_",num2str(xbin)," regions"),"Regions","PI:v_x/sp",color,name,"line",'ii_num',ii_num,'ebar',1,'ylim',[-0.2 0.3]);
    mlt_subplt(PI_y,[1:ybin]',200,1,1,1,append("Preferential Index for different y regions_",num2str(xbin)," regions"),"Regions","PI:v_y/sp",color,name,"line",'ii_num',ii_num,'ebar',1,'ylim',[-0.2 0.3]);

end

%% compute the PI (substract mean of ctr from mean of exp)
% if ~isempty(input_cond)
%     idx=find(contains(input_cond,"ctr"));
%     if xbin<=6
%         [data_x_sub,name1,color1]=sub_ctr_data(PI_x_mean,idx,name,color,'sem',PI_x_sem,'input_cond',input_cond);
%         mlt_subplt(data_x_sub,[],160,1,1,1,append("Preferential Index for different regions_",num2str(xbin)," regions,substract mean of the control group"),"Groups","PI:v_x/sp (exp-ctr)",color1,name1,"bar",'ii_num',length(name1),'title_lab','Different positions','ebar',1,'legends',xleg);
%     else
%         [data_x_sub,name1,color1]=sub_ctr_data(PI_x,idx,name,color,'input_cond',input_cond);
%         mlt_subplt(data_x_sub,[1:xbin]',160,1,1,1,append("Preferential Index for different regions_",num2str(xbin)," regions,substract mean of the control group"),"Regions","PI:v_x/sp (exp-ctr)",color1,name1,"line",'ii_num',length(name1),'ebar',1);
%
%     end
%     clear idx
% end
%% plot the preferential index vs different gradient either opto or odor
% if sti=="opto"
%     grad_series=[grad_min+0.05:0.1:grad_max-0.05]';
%     grad_series1=[grad_d_min+0.0005:0.001:grad_d_max-0.0005]';
%     mlt_subplt(PI_grad,grad_series,153,1,1,1,"Preferential Index across light intensity","Light Intensity (uW/mm^2)","PI:v_x/sp",color,name,"line",'ii_num',ii_num,'ebar',1);
%     mlt_subplt(PI_grad_diff,grad_series1,154,1,1,1,"Preferential Index across light gradient","Light Gradient (uW/mm^2)","PI:v_x/sp",color,name,"line",'ii_num',ii_num,'ebar',1);
%
%     if ~isempty(input_cond)
%         idx=find(contains(input_cond,"ctr"));
%
%         [data_grad_sub,name1,color1]=sub_ctr_data(PI_grad,idx,name,color,'input_cond',input_cond);
%         mlt_subplt(data_grad_sub,grad_series,161,1,1,1,append("Preferential Index for different ",sti," intensity "," ,substract mean of the control group"),"Light Intensity (uW/mm^2)","PI:v_x/sp (exp-ctr)",color1,name1,"line",'ii_num',length(name1),'ebar',1);
%
%         [data_grad_diff_sub,name1,color1]=sub_ctr_data(PI_grad_diff,idx,name,color,'input_cond',input_cond);
%         mlt_subplt(data_grad_diff_sub,grad_series1,162,1,1,1,append("Preferential Index for different ",sti," gradient "," ,substract mean of the control group"),"Light Gradient (uW/mm^2)","PI:v_x/sp (exp-ctr)",color1,name1,"line",'ii_num',length(name1),'ebar',1);
%         clear idx
%     end
% end
%% plot the preferential index for each region across time
% for i=1:xbin
%     %remember to change the series if we change the one for calculate the
%     %PI_xt
%     mlt_subplt(PI_xt{i,1},[2.5:5:tmax-2.5]',111+i,1,1,1,append("Preferential Index for region_",num2str(i),"_across time"),"Time","PI:v_x/sp",color,name,"line",'ii_num',ii_num,'ebar',1);
% end
%% plot the preferential index across time (using the position)
mlt_subplt(PI_timeseries,time1,150,1,1,1,"Preferential Index based on position","Time(s)","Preferential index: (N_p_r_o_x-N_d_i_s_t)/(N_p_r_o_x+N_d_i_s_t)",color,name,"line",'ii_num',ii_num,'ebar',1,'ylim',[-1 1]);

%% plot dist to odor vs time
    t_series1=[0.5:1:tmax-0.5]';
    for i=1:length(dis_t)
        dis_t1{i,1}(:,1)=dis_t{i,1}(:,1)-dis_t{i,1}(1,1);
        dis_t1{i,1}(:,2)=dis_t{i,1}(:,2);
    end
    mlt_subplt(dis_t1,t_series1,151,1,1,1,'Distant to Odor','Time(s)','Distance to odor (mm)',color,name,'line','ebar',1,'ii_num',ii_num,'ylim',[-70 30],'ydir',1);

%% plot th vx, speed and PI vs time
mlt_subplt(vx,t_series1,152,2,2,1,'','Time(s)','v_x(mm/s)',color,name,'line','ebar',1,'ii_num',ii_num);
mlt_subplt(sp,t_series1,152,2,2,2,'','Time(s)','speed(mm/s)',color,name,'line','ebar',1,'ii_num',ii_num);
mlt_subplt(PI_tpl,t_series1,152,2,2,3,'vx,speed, and preferential index vs time','Time(s)','Preferential index (v_x/sp)',color,name,'line','ebar',1,'ii_num',ii_num);

%save figures
save_all_figures(outdir);
close all

%% stat test
%% 1) first, test the general PI for each group
str1=vertcat(str1{:});
[p_mwu,group_comp,p_k]=non_parametric_test(PI_larvae,str1,PI_size);
figgcf;
n=fig.Number;
outdir1=fullfile(outdir,"Stat_result");
change_fig_prop(n,'title',"Kruskal Wallis test of the preferential index of different conditions",'xtitle',"Groups",'ytitle',"PI(v_x/sp)", ...
    'name',"Kruskal_Wallis_test_groups_boxplot");
change_fig_prop(n-1,'name',"Kruskal_Wallis_test_groups_table");
%save figures
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all
% same thing but for vy/sp
str1_y=vertcat(str1_y{:});
[p_mwu_y,group_comp_y,p_k_y]=non_parametric_test(PI_larvae_y,str1_y,PI_size_y_all);
fig=gcf;
n=fig.Number;
change_fig_prop(n,'title',"Kruskal Wallis test of the preferential index of different conditions for vy",'xtitle',"Groups",'ytitle',"PI(v_y/sp)", ...
    'name',"Kruskal_Wallis_test_groups_boxplot");
change_fig_prop(n-1,'name',"Kruskal_Wallis_test_groups_table_for_vy");
%save figures

save_all_figures(outdir1);
close all

%% 2) check the PI for each timebin
outdir2=fullfile(outdir1,append("Duration_",num2str(t),"s"));
for i=1:length(size_t)
    size=vertcat(size_t{i,1}{:});
    PI_larva=vertcat(PI_larva_t_all{i,1}{:});
    group=vertcat(group_t_all{i,1}{:});
    group=vertcat(group{:});

    [p_mwut{i,1},group_compt{i,1},p_kt(i,1)]=non_parametric_test(PI_larva,group,size);
    fig=gcf;
    n=fig.Number;
    change_fig_prop(n,'title',append("Kruskal Wallis test of the preferential index of different conditions from ",num2str(dur(i)),"-",num2str(dur(i+1)),"s"),...
        'xtitle',"Groups",'ytitle',"PI(v_x/sp)",'name',append("Kruskal_Wallis_test_groups_duration_",num2str(dur(i)),"-",num2str(dur(i+1)),"s_boxplot"));
    change_fig_prop(n-1,'name',append("Kruskal_Wallis_test_groups_duration_",num2str(dur(i)),"-",num2str(dur(i+1)),"s_table"));
end
%save figures
if ~isfolder(outdir2)
    mkdir(outdir2);
end
save_all_figures(outdir2);
close all
%% 3) check the PI for each xbin
outdir3=fullfile(outdir1,append(num2str(xbin),"_regions"));

for i=1:length(size_x)
    size=vertcat(size_x{i,1}{:});
    PI_larva=vertcat(PI_larva_x_all{i,1}{:});
    group=vertcat(group_x_all{i,1}{:});
    group=vertcat(group{:});

    [p_mwux{i,1},group_compx{i,1},p_kx(i,1)]=non_parametric_test(PI_larva,group,size);

    fig=gcf;
    n=fig.Number;
    change_fig_prop(n,'title',append("Kruskal Wallis test of the preferential index of different region ",num2str(i)),...
        'xtitle',"Groups",'ytitle',"PI(v_x/sp)",'name',append("Kruskal_Wallis_test_groups_region_",num2str(i),"_boxplot"));
    change_fig_prop(n-1,'name',append("Kruskal_Wallis_test_groups_region_",num2str(i),"_table"));

end
%save figures
if ~isfolder(outdir3)
    mkdir(outdir3);
end
save_all_figures(outdir3);
close all
%% 3) check the PI for each intensity
%  if strcmp(sti,"opto")
%     outdir4=fullfile(outdir1,append(num2str(length(size_grad)),"_intensity"));
%     grad_series=[grad_min:0.1:grad_max]';
%     for i=1:length(size_grad)
%         size=vertcat(size_grad{i,1}{:});
%         PI_larva=vertcat(PI_larva_grad_all{i,1}{:});
%         group=vertcat(group_grad_all{i,1}{:});
%         group=vertcat(group{:});
%
%         [p_mwu_int{i,1},group_comp_int{i,1},p_k_int(i,1)]=non_parametric_test(PI_larva,group,size);
%
%         fig=gcf;
%         n=fig.Number;
%         change_fig_prop(n,'title',append("Kruskal Wallis test of the preferential index of",sti," intensity,",num2str(grad_series(i)),"-",num2str(grad_series(i+1)),"uW_per_mm2"),...
%             'xtitle',"Groups",'ytitle',"PI(v_x/sp)",'name',append("Kruskal_Wallis_test_groups_",sti,"_intensity_",num2str(grad_series(i)),"-",num2str(grad_series(i+1)),"uW_per_mm2_boxplot"));
%         change_fig_prop(n-1,'name',append("Kruskal_Wallis_test_groups_intensity_",num2str(grad_series(i)),"-",num2str(grad_series(i+1)),"uW_permm2_table"));
%
%     end
%     %save figures
%     if ~isfolder(outdir4)
%         mkdir(outdir4);
%     end
%     save_all_figures(outdir4);
%     close all
%     %% 4) check the PI for each light gradient
%     outdir5=fullfile(outdir1,append(num2str(length(size_grad_d)),"_gradient"));
%     grad_d_series=[grad_d_min:0.001:grad_d_max]';
%     for i=1:length(size_grad_d)
%         size=vertcat(size_grad_d{i,1}{:});
%         PI_larva=vertcat(PI_larva_grad_d_all{i,1}{:});
%         group=vertcat(group_grad_d_all{i,1}{:});
%         group=vertcat(group{:});
%
%         [p_mwu_grad{i,1},group_comp_grad{i,1},p_k_grad(i,1)]=non_parametric_test(PI_larva,group,size);
%
%         fig=gcf;
%         n=fig.Number;
%         change_fig_prop(n,'title',append("Kruskal Wallis test of the preferential index of",sti," gradient,",num2str(grad_d_series(i)),"-",num2str(grad_d_series(i+1)),"uW_per_mm2"),...
%             'xtitle',"Groups",'ytitle',"PI(v_x/sp)",'name',append("Kruskal_Wallis_test_groups_",sti,"_gradient_",num2str(grad_d_series(i)),"-",num2str(grad_d_series(i+1)),"uW_per_mm2_boxplot"));
%         change_fig_prop(n-1,'name',append("Kruskal_Wallis_test_groups_gradient_",num2str(grad_d_series(i)),"-",num2str(grad_d_series(i+1)),"uW_per_mm2_table"));
%
%     end
%     %save figures
%     if ~isfolder(outdir5)
%         mkdir(outdir5);
%     end
%     save_all_figures(outdir5);
%     close all
% end
%% save data

filename=fullfile(outdir1,'p_value.mat');

if isfile(filename)
    delete(filename);
end
save(filename,'p_k*','*timeseries','p_mwu*','group_comp*');
clear filename

close all
clear
end
