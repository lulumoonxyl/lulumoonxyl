function []=single_animal_analysis(tracker_num,genotype,condition,filename,output_name,varargin)
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

input_cond=[];%it will have a same length as name, telling us whether the input is 'ctr' or 'exp'
row_comp=2; col_comp=2; % these two will be used for the plot to show the comparison of ctr vs exp group, this number should be based on how many comparison you have
tmax=900;% tbin is used for binning the turn event's t0s
xmax=225; deg_list=60;

%cond_list is a arry of 1 and 0, 1 means attractive while 0 is
%aversive

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
    et=dat.et;
    x=dat.x;
    y=dat.y
    for i=1:length(et)
        et_len(i,1)=et{i,1}(end)-et{i,1}(1);
    end

    idx=find(et_len>=100);
    for j=1:length(idx)-10
        figure(1+j)
        hold on
        plot(x{idx(j),1},y{idx(j),1},'r.');
        text(x{idx(j),1}(1,1),y{idx(j),1}(1,1),'start')
        xlim([0 210]);
        ylim([0 210]);
        hold off
    end
end
end