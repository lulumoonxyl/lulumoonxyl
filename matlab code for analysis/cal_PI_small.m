function []=cal_PI_small(tracker_num,genotype,condition,filename,output_name,varargin)
%cal PI for exp using small agar plate, location of odor is a dot instead
%of a line
w=1;
for i=1:length(genotype)
    for j=1:length(condition)
        cond(w,1)=fullfile("/project/6010970/screen/olfactory_output/choreography",tracker_num,genotype{i,1},condition{j,1},filename);
        name(w,1)=append(genotype{i,1},"@",condition{j,1});
        w=w+1;
    end
end

outdir=fullfile("/project/6010970/screen/olfactory_output/Preferential index",output_name)
%set the properties for plots
row=2;col=3;
color=[0.2,0.2,0.1;0.3,0.75,0.93;0.47,0.67,0.19;0.53,0.53,0.02;1,0.41,0.16]
% check the traj plot to input xmax, ymax
xmax=85;ymax=85;tmax=300;bin=5;tbin=10;traj=1;hp=1;
odor_pos=[83 42];

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
    elseif strcmp(varargin{i},'plot traj')
        traj=varargin{i+1};
    elseif strcmp(varargin{i},'plot heatmap')
        hp=varargin{i+1};
    elseif strcmp(varargin{i},'odor pos')
        odor_pos=varargin{i+1};
    elseif strcmp(varargin{i},'name')
        clear name
        name=string(varargin{i+1})
        
    end
end

xseries=[0:bin:xmax]'; x_series=[0+bin/2:bin:xmax-bin/2]';
yseries=[0:bin:ymax]'; y_series=[0+bin/2:bin:ymax-bin/2]';
tseries=[0:tbin:tmax]';  t_series=[0+tbin/2:tbin:tmax-tbin/2]';
ii_num=length(cond)
for ii=1:ii_num
    load(cond(ii));
    
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
        dis_to{i,1}=sqrt((odor_pos(1)-dat_grouped.x{i,1}).^2+(odor_pos(2)-dat_grouped.y{i,1}).^2);
    end
    
    if ii==1
        dis_t=cell(ii_num,1);
    end
    tseries1=0:1:tmax;
    dis_t{ii,1}=bin_data_mean(vertcat(dat_grouped.et{:}),tseries1,vertcat(dis_to{:}));
    
    %% 4) Calculate the navigational index using cos or sin(theta)
    % 4-1)calculate the navigational index, using the odor_pos and the sin theta
    if ~isempty(odor_pos)
        for i=1:length(dat_grouped.x)
            x1=dat_grouped.x{i,1}(2:end);
            x2=dat_grouped.x{i,1}(1:end-1);
            
            y1=dat_grouped.y{i,1}(2:end);
            y2=dat_grouped.y{i,1}(1:end-1);
            
            vector(:,1)=odor_pos(1)-x2;
            vector(:,2)=odor_pos(2)-y2;
            deg=cal_ori_deg(x1,x2,y1,y2,vector);
            theta=deg2rad(deg);
            
            PI_cos{i,1}=cos(theta);
            % Preferential index toward or away odorant
            PI_theta_cos(i,1)=mean(PI_cos{i,1});
            
            PI_sin{i,1}=sin(theta);
            % Preferential index perpendicular to the odor pos
            PI_theta_sin(i,1)=mean(PI_sin{i,1});
            clear x1 x2 y1 y2 theta deg vector
        end
    end
    
    if ii==1
        PI_mean_theta_cos=zeros(ii_num,2);
        PI_mean_theta_sin=zeros(ii_num,2);
        
        str_cos=[];str_sin=[];
        PI_larvae_cos=[]; PI_larvae_sin=[];
        PI_size_cos=[];PI_size_sin=[];
    end
    
    PI_mean_theta_cos(ii,1)=mean(PI_theta_cos);
    PI_mean_theta_cos(ii,2)=std(PI_theta_cos)/sqrt(length(PI_theta_cos));
    
    PI_mean_theta_sin(ii,1)=mean(PI_theta_sin);
    PI_mean_theta_sin(ii,2)=std(PI_theta_sin)/sqrt(length(PI_theta_sin));
    
    %% 4-2) save data for stat test of PI
    PI_larvae_sin=vertcat(PI_larvae_sin,PI_theta_sin);
    PI_size_sin(ii,1)=length(PI_theta_sin);
    str_sin=vertcat(str_sin,repmat({name(ii)},PI_size_sin(ii,1),1));
    
    PI_larvae_cos=vertcat(PI_larvae_cos,PI_theta_cos);
    PI_size_cos(ii,1)=length(PI_theta_cos);
    str_cos=vertcat(str_cos,repmat({name(ii)},PI_size_cos(ii,1),1));
    
    clear PI PI_theta_cos PI_theta_sin PI_sin PI_cos
    %% 5) Compute the PI using position
    if ii==1
        p_timeseries=cell(ii_num,1);PI_timeseries=cell(ii_num,1);
        p_timeseries_y=cell(ii_num,1);PI_timeseries_y=cell(ii_num,1);
    end
    
    timebins=2;time=0:timebins:tmax;
    time1=[0+timebins/2:timebins:tmax-timebins/2]';
    edge=[0 xmax/2-xmax*0.05 xmax/2+xmax*0.05 xmax];
    %separate x pos into three large bins and count the animal in that bins every 5s (timebins)
    edge_y=[0 ymax/2-ymax*0.05 ymax/2+ymax*0.05 ymax];
    [p_timeseries{ii,1},PI_timeseries{ii,1}]=cal_navindex_vs_t(dat_grouped.x,dat_grouped.et,time,edge);
    [p_timeseries_y{ii,1},PI_timeseries_y{ii,1}]=cal_navindex_vs_t(dat_grouped.y,dat_grouped.et,time,edge_y);
    
    color3=[0.5 1 0;0.5 0 1;1 0 0];
    mlt_subplt(p_timeseries{ii,1},time1,105,row,col,ii,'Frequency of animal distribution across time','Time(s)','Frequency',color3,name(ii),'multiple lines','legends',["Distant";"Central";"Proximate"],'ii_num',ii_num,'lw',2,'ylim',[0 1]);
    mlt_subplt(p_timeseries_y{ii,1},time1,101,row,col,ii,'Frequency of animal distribution across time in y region','Time(s)','Frequency',color3,name(ii),'multiple lines','legends',["Lower";"Central";"Upper"],'ii_num',ii_num,'lw',2,'ylim',[0 1]);
    %% keep certain variables
    clearvars -except PI* ii* lim* str* *series outdir str color* name t time* *bin dis_* cond *max hp traj t odor_pos row col c_tra size*
end


%% plot the number of larvae tracked across time
mlt_subplt(c_tra,[1/2:1:tmax-1/2]',108,1,1,1,'Number of objects tracked across time','Time(s)','Number of objects tracked',color,'','multiple lines','legends',name);
%% plot the preferential index across time (using the position)
mlt_subplt(PI_timeseries,time1,109,1,1,1,"Preferential Index based on position","Time(s)","Preferential index: (N_p_r_o_x-N_d_i_s_t)/(N_p_r_o_x+N_d_i_s_t)",color,name,"line",'ii_num',ii_num,'ebar',1,'ylim',[-1 1]);

%% plot dist to odor vs time
t_series1=[0.5:1:tmax-0.5]';
for i=1:length(dis_t)
    dis_t1{i,1}(:,1)=dis_t{i,1}(:,1)-dis_t{i,1}(1,1);
    dis_t1{i,1}(:,2)=dis_t{i,1}(:,2);
end
mlt_subplt(dis_t1,t_series1,110,1,1,1,'Distant to Odor','Time(s)','Distance to odor (mm)',color,name,'line','ebar',1,'ii_num',ii_num,'ylim',[-30 30],'ydir',1);

%% plot the preferenrial index
mlt_subplt(PI_mean_theta_cos,[],111,1,1,1,"Preferential Index (cos theta)","Groups","Preferential Index:cos(theta)",color,name,"bar",'ii_num',ii_num,'ebar',1,'ylim',[-0.25 0.25]);
mlt_subplt(PI_mean_theta_sin,[],112,1,1,1,"Preferential Index (sine theta)","Groups","Preferential Index:sin(theta)",color,name,"bar",'ii_num',ii_num,'ebar',1,'ylim',[-0.25 0.25]);

%save figures
if ~isfolder(outdir)
    mkdir(outdir);
end
save_all_figures(outdir);
close all
%% stat test

str_cos=vertcat(str_cos{:});
[p_mwu_cos,group_comp_cos,p_k_cos]=non_parametric_test(PI_larvae_cos,str_cos,PI_size_cos,0);

str_sin=vertcat(str_sin{:});
[p_mwu_sin,group_comp_sin,p_k_sin]=non_parametric_test(PI_larvae_sin,str_sin,PI_size_sin,0);

outdir1=fullfile(outdir,"Stat_result");
filename=fullfile(outdir1,'p_value.mat');
if ~isfolder(outdir1)
    mkdir(outdir1);
end
if isfile(filename)
    delete(filename);
end
save(filename,'p_k*','*timeseries','p_mwu*','group_comp*');
clear filename

close all
clear
end