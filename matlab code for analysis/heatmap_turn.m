function []=heatmap_turn(tracker_num,genotype,condition,filename,output_name,varargin)
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
outdir=fullfile("/project/6010970/screen/olfactory_output/heatmap about turning events",output_name);
ii_num=length(cond);
%% set the properties for plots
row=3;col=3;
color=[0.75,0.75,0.75;0.627,0.82,1;0.659,0.659,0.659;0.871,0.416,0.451];
xmax=225;ymax=225;tmax=900;
bin=10;tbin=30;hbin=10;% for binning data
deg_list=60;

for i=1:2:length(varargin)
    if strcmp ('row',varargin{i})
        row=varargin{i+1};
    elseif strcmp ('col',varargin{i})
        col=varargin{i+1};
    elseif strcmp ('tbin',varargin{i})
        tbin=varargin{i+1};
    elseif strcmp ('hbin',varargin{i})
        hbin=varargin{i+1};
    elseif strcmp ('tmax',varargin{i})
        tmax=varargin{i+1};
    elseif strcmp ('color',varargin{i})
        color=varargin{i+1};
    elseif strcmp(varargin{i},'name')
        clear name;
        name=string(varargin{i+1})
    elseif strcmp ('xmax',varargin{i})
        xmax=varargin{i+1};
    elseif strcmp ('ymax',varargin{i})
        ymax=varargin{i+1};
    elseif strcmp ('bin',varargin{i})
        bin=varargin{i+1};
    elseif strcmp ('deg_l',varargin{i})
        deg_list=varargin{i+1};
    end
end
xseries=[0:bin:xmax]'; x_series=[0+bin/2:bin:xmax-bin/2]';
yseries=[0:bin:ymax]'; y_series=[0+bin/2:bin:ymax-bin/2]';
tseries=[0:tbin:tmax]';  t_series=[0+tbin/2:tbin:tmax-tbin/2]';
hdseries=[-180:hbin:180]'; hd_series=[-180+hbin/2:hbin:180-hbin/2]';

for ii=1:ii_num
    load(cond(ii));
    %% 1-1) Compute probability of orientation spend in x/y/t bins
    [t_x,p_x,p_bin_x]=bin_data_sum_t2(dat.orientation,dat.x,hdseries,xseries,dat.et);
    [t_y,p_y,p_bin_y,]=bin_data_sum_t2(dat.orientation,dat.y,hdseries, yseries,dat.et);
    [t_t,p_t,p_bin_t]=bin_data_sum_t2(dat.orientation,dat.et,hdseries,tseries,dat.et);
    %save all the data for plotting later
    %fig_num=[1:13]';
    %% 1-2) Plot-->lim is saved for changing the caxis of heatmap later
    lim{1,1}(ii,:)=mlt_subplt(flipud(p_x),x_series',1,row,col,ii,append('Heatmap@relative p of orientation vs X pos@bin=',num2str(bin)), ...
        'X(mm)','Heading Direction (deg)',color,append('p_o_r_i_e_n_t_a_t_i_o_n-',name(ii)),'heatmap','series2',hd_series');
    
    lim{2,1}(ii,:)=mlt_subplt(flipud(p_y),y_series',2,row,col,ii,append('Heatmap@relative p of orientation vs Y pos@bin=',num2str(bin)), ...
        'Y(mm)','Heading Direction (deg)',color,append('p_o_r_i_e_n_t_a_t_i_o_n-',name(ii)),'heatmap','series2',hd_series');
    
    lim{3,1}(ii,:)=mlt_subplt(flipud(p_t),t_series',3,row,col,ii,append('Heatmap@relative p of orientation vs TIME@bin=',num2str(bin)), ...
        'Time(s)','Heading Direction (deg)',color,append('p_o_r_i_e_n_t_a_t_i_o_n-',name(ii)),'heatmap','series2',hd_series');
    
    lim{4,1}(ii,:)=mlt_subplt(flipud(p_bin_x),x_series',4,row,col,ii,append('Heatmap@relative p of orientation vs X pos in each bin@bin=',num2str(bin)), ...
        'X(mm)','Heading Direction (deg)',color,append('p_o_r_i_e_n_t_a_t_i_o_n-X bin',name(ii)),'heatmap','series2',hd_series');
    
    lim{5,1}(ii,:)=mlt_subplt(flipud(p_bin_y),y_series',5,row,col,ii,append('Heatmap@relative p of orientation vs Y pos in each bin@bin=',num2str(bin)), ...
        'Y(mm)','Heading Direction (deg)',color,append('p_o_r_i_e_n_t_a_t_i_o_n-Y bin',name(ii)),'heatmap','series2',hd_series');
    
    lim{6,1}(ii,:)=mlt_subplt(flipud(p_bin_t),t_series',6,row,col,ii,append('Heatmap@relative p of orientation vs TIME in each bin@bin=',num2str(bin)), ...
        'Time(s)','Heading Direction (deg)',color,append('p_o_r_i_e_n_t_a_t_i_o_n-Time bin',name(ii)),'heatmap','series2',hd_series');
   clear p_*
   %% 2-1) Compute the turning freq vs HD for different x/y/t bins
    c_x=bin_data_count2(vertcat(dat.pre_deg{:}),hdseries,vertcat(dat.turn_x{:}),xseries);
    pc_x=c_x./t_x;
    c_y=bin_data_count2(vertcat(dat.pre_deg{:}),hdseries,vertcat(dat.turn_y{:}),yseries);
    pc_y=c_y./t_y;
    t_ver=vertcat(dat.t0s{:});
    if isa(t_ver,'cell')
        t_ver=vertcat(t_ver{:});
    end 
    c_t=bin_data_count2(vertcat(dat.pre_deg{:}),hdseries,t_ver,tseries);
    pc_t=c_t./t_t;
     %% 2-2) Plot-->lim is saved for changing the caxis of heatma later
    lim{7,1}(ii,:)=mlt_subplt(flipud(pc_x),x_series',7,row,col,ii,append('Heatmap@turning freq vs HD vs X pos@bin=',num2str(bin)), ...
        'X(mm)','Heading Direction (deg)',color,append('Turn freq-',name(ii)),'heatmap','series2',hd_series');
    
    lim{8,1}(ii,:)=mlt_subplt(flipud(pc_y),y_series',8,row,col,ii,append('Heatmap@turning freq vs HD vs Y pos@bin=',num2str(bin)), ...
        'Y(mm)','Heading Direction (deg)',color,append('Turn freq-',name(ii)),'heatmap','series2',hd_series');
    
    lim{9,1}(ii,:)=mlt_subplt(flipud(pc_t),t_series',9,row,col,ii,append('Heatmap@turning freq vs HD vs TIME @bin=',num2str(bin)), ...
        'Time(s)','Heading Direction (deg)',color,append('Turn freq-',name(ii)),'heatmap','series2',hd_series');
    %% 2-3) Compute the larger turning freq vs HD in different x/y/y bins
    pre=vertcat(dat.pre_deg{:}); post=vertcat(dat.post_deg{:});
    turn_x=vertcat(dat.turn_x{:});turn_y=vertcat(dat.turn_y{:});
    ang=vertcat(dat.reorient_deg_abs{:});
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
    idx=find(ang>=deg_l);

    c_ang_x=bin_data_count2(pre(idx),hdseries,turn_x(idx),xseries);
    c_ang_y=bin_data_count2(pre(idx),hdseries,turn_y(idx),yseries);
    c_ang_t=bin_data_count2(pre(idx),hdseries,t0s(idx),tseries);
    
    pc_ang_x=c_ang_x./t_x;
    pc_ang_y=c_ang_y./t_y;
    pc_ang_t=c_ang_t./t_t;
    %% 2-4) plot the heatmap
    lim{13,1}(ii,:)=mlt_subplt(flipud(pc_ang_x),x_series',13,row,col,ii,append('Heatmap@Large turn freq vs HD vs X pos@bin=',num2str(bin)), ...
        'X(mm)','Heading Direction (deg)',color,append('Large turn freq-',name(ii)),'heatmap','series2',hd_series');
    
    lim{14,1}(ii,:)=mlt_subplt(flipud(pc_ang_y),y_series',14,row,col,ii,append('Heatmap@Large turn freq vs HD vs Y pos@bin=',num2str(bin)), ...
        'Y(mm)','Heading Direction (deg)',color,append('Large turn freq-',name(ii)),'heatmap','series2',hd_series');
    
    lim{15,1}(ii,:)=mlt_subplt(flipud(pc_ang_t),t_series',15,row,col,ii,append('Heatmap@Large turn freq vs HD vs TIME @bin=',num2str(bin)), ...
        'Time(s)','Heading Direction (deg)',color,append('Large turn freq-',name(ii)),'heatmap','series2',hd_series');
 
    
    %% 3-1) Compute the reorientation angle in different HD and y/x/t bin
    [~,deg_x]=bin_data_mean2(vertcat(dat.pre_deg{:}),hdseries,vertcat(dat.turn_x{:}),xseries,vertcat(dat.reorient_deg_abs{:}));
    [~,deg_y]=bin_data_mean2(vertcat(dat.pre_deg{:}),hdseries,vertcat(dat.turn_y{:}),yseries,vertcat(dat.reorient_deg_abs{:}));
    [~,deg_t]=bin_data_mean2(vertcat(dat.pre_deg{:}),hdseries,t_ver,tseries,vertcat(dat.reorient_deg_abs{:}));
    %% 3-2) Plot heatmap
    lim{10,1}(ii,:)=mlt_subplt(flipud(deg_x),x_series',10,row,col,ii,append('Heatmap@Reorientation angle vs HD vs X pos@bin=',num2str(bin)), ...
        'X(mm)','Heading Direction (deg)',color,append('Reorientation angle-',name(ii)),'heatmap','series2',hd_series');
    
    lim{11,1}(ii,:)=mlt_subplt(flipud(deg_y),y_series',11,row,col,ii,append('Heatmap@Reorientation angle vs HD vs Y pos@bin=',num2str(bin)), ...
        'Y(mm)','Heading Direction (deg)',color,append('Reorientation angle-',name(ii)),'heatmap','series2',hd_series');
    
    lim{12,1}(ii,:)=mlt_subplt(flipud(deg_t),t_series',12,row,col,ii,append('Heatmap@Reorientation angle vs HD vs TIME @bin=',num2str(bin)), ...
        'Time(s)','Heading Direction (deg)',color,append('Reorientation angle-',name(ii)),'heatmap','series2',hd_series');

   %% Clear variables
    clearvars -except row col name cond color *max *bin deg* outdir output_name *_series lim ii* fig_num *series
end
%% 1-3) & 2-5) Organized the heatmap caxis
for i=1:length(lim)
    axis_max=max(lim{i,1}(:,2));
    axis_min=min(lim{i,1}(:,1));

        for k=1:ii_num
            figure(i)
            subplot(row,col,k)
            caxis([axis_min axis_max]);
        end

end
%% save figures
save_all_figures(outdir);
close all
clear

end
