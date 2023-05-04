function []=cal_PI_by_run(tracker_num,genotype,condition,filename,output_name,varargin)
%% get the directory for all the data
w=1;

for i=1:length(genotype)
    for j=1:length(condition)
        cond(w,1)=fullfile("/project/6010970/screen/olfactory_output/JB_JAABA",tracker_num,genotype{i,1},condition{j,1},filename);
        
%         cond(w,1)=fullfile("C:\Users\feihu\OneDrive - McGill University\matlab\Analysis based on larvae\JB_JAABA",tracker_num,genotype{i,1},condition{j,1},filename);
%       
        name(w,1)=append(genotype{i,1},"@",condition{j,1})
        w=w+1;
    end
end
%% set the output name
outdir=fullfile("/project/6010970/screen/olfactory_output/Preferential index/run events",output_name);
% outdir=fullfile("C:\Users\feihu\OneDrive - McGill University\matlab\Analysis based on larvae\JB_JAABA",output_name);

ii_num=length(cond);

row=3;col=3;
color=[0,0,0;0.659,0.659,0.659;0.627,0.82,1;0.439,0.733,1;0.255,0.643,1;0.871,0.416,0.451;1,0,0];
grad_min=0; %the min value of gradient (opto or odor)
grad_max=3;
xmax=225; xbin=3;
sti='odor';
vector=[];%vector indicates the pos of the highest stimulus, should be consistent with the vector input for the group_JB_JAABA_data() function

for i=1:2:length(varargin)
    if strcmp ('row',varargin{i})
        row=varargin{i+1};
    elseif strcmp ('col',varargin{i})
        col=varargin{i+1};
    elseif strcmp ('color',varargin{i})
        color=varargin{i+1};
    elseif strcmp(varargin{i},'grad_min')
        grad_min=varargin{i+1};
    elseif strcmp(varargin{i},'grad_max')
        grad_max=varargin{i+1};
    elseif strcmp(varargin{i},'xbin')
        xbin=varargin{i+1};
    elseif strcmp (varargin{i},'xmax')
        xmax=varargin{i+1};
    elseif strcmp(varargin{i},'stimulus')
        sti=varargin{i+1};
    elseif strcmp(varargin{i},'name')
        name=string(varargin{i+1});
    elseif strcmp(varargin{i},'vector')
        vector=varargin{i+1};

    end
end

for ii=1:ii_num
    %% load data
    load(cond(ii));
    ii
    %% 1) Compute the PI for each running event using run_x run_y run_et
    % 1-1) get the x,y,t of each running event for both the start and
    % stop using center of mass
    run_et=vertcat(dat.run_et{:});
    r0x=vertcat(dat.r0x{:});
    
    for i=1:height(dat)
        r0_idx=dat.r0_idx{i,1};
        r1_idx=dat.r1_idx{i,1};
        y=dat.y{i,1};
        x=dat.x{i,1};
        
        r1x{i,1}=x(r1_idx);
        r0y{i,1}=y(r0_idx);
        r1y{i,1}=y(r1_idx);
        clear y idx *_idx
    end
    
    r1x=vertcat(r1x{:});
    r0y=vertcat(r0y{:});
    r1y=vertcat(r1y{:});
    % 1-2)compute the PI for each running event
    x_diff=r1x-r0x;
    y_diff=r1y-r0y;
    dist=sqrt(x_diff.^2+y_diff.^2);
    sp=dist./run_et;
    vx=x_diff./run_et;
    PI=vx./sp;
    
    % 1-3)save the mean data
    if ii==1
        PI_mean=zeros(ii_num,2);
        PI_larva=[];group=[];
        PI_size=zeros(ii_num,1);
    end
    PI_mean(ii,1)=mean(PI,'omitnan');
    PI_mean(ii,2)=std(PI,'omitnan')/length(PI);
    
    % 1-4) Save the data for stat
    PI_larva=vertcat(PI_larva,PI);
    PI_size(ii,1)=length(PI);
    group=vertcat(group,repmat({name(ii)},PI_size(ii,1),1));
    
    
    %% 2) Bin the PI based on intensity
    if strcmp(sti,'opto')
        if ii==1
            PI_mean_int=cell(ii_num,1);
        end
        int=vertcat(dat.grad_r0{:});
        int_series=[grad_min:0.1:grad_max]';
        [PI_mean_int{ii,1}]=bin_data_mean(int,int_series,PI);
    end
    %% 3) Bin the PI based on regions
    
    xseries=[0:xmax/xbin:xmax];
    [a]=bin_data_mean(r0x,xseries,PI);
    if xbin<=6
        if ii==1
            PI_mean_x=zeros(ii,xbin);PI_sem_x=zeros(ii,xbin);
        end
        PI_mean_x(ii,:)=a(:,1)';PI_sem_x(ii,:)=a(:,2)';
    else
        if ii==1
            PI_mean_x=cell(ii_num,1);
        end
        PI_mean_x{ii,1}=a;
    end
    clear a PI
    %% 4) Compute the PI using theta
     for i=1:length(dat.x)
        x1=dat.x{i,1}(2:end);
        x2=dat.x{i,1}(1:end-1);

        y1=dat.y{i,1}(2:end);
        y2=dat.y{i,1}(1:end-1);
        deg=cal_ori_deg(x1,x2,y1,y2,vector);
        theta=deg2rad(deg);
        PI{i,1}=cos(theta);
        PI_theta(i,1)=mean(PI{i,1});
        clear x1 x2 y1 y2
    end 
    
    if ii==1
        PI_mean_theta=zeros(ii_num,2);
    end 
    PI_mean_theta(ii,1)=mean(PI_theta);
    PI_mean_theta(ii,2)=std(PI_theta)/length(PI_theta);
    clear PI PI_theta
    %%
    clearvars -except group PI* dat outdir row col name color cond ii* group grad* xbin xmax sti vector
end

% 1-6) Plot the PI
mlt_subplt(PI_mean,[],100,1,1,1,"Preferential Index","Groups","PI:v_x/sp",color,name,"bar",'ii_num',ii_num,'ebar',1,'ylim',[-0.1 0.1]);
mlt_subplt(PI_mean_theta,[],104,1,1,1,"Preferential Index,using theta","Groups","PI:cos(theta)",color,name,"bar",'ii_num',ii_num,'ebar',1,'ylim',[-0.1 0.1]);
if strcmp (sti,'opto')
    int_series=[grad_min+0.05:0.1:grad_max-0.05]';
    mlt_subplt(PI_mean_int,int_series,101,1,1,1,"Preferential Index across light intensity using running event","Light Intensity (uW/mm^2)","PI:v_x/sp",color,name,"line",'ii_num',ii_num,'ebar',1);
end
if xbin<=6
    
    for i=1:xbin
        xleg(1,i)=append("Region ", num2str(i));
    end
    data_x={PI_mean_x,PI_sem_x};
    mlt_subplt(data_x,[],103,1,1,1,append("Preferential Index for different regions_",num2str(xbin)," regions, using running event"),"Groups","PI:v_x/sp",color,name,"bar",'ii_num',ii_num,'ebar',1,'legends',xleg);
    
else
    a=xmax/xbin;
    x_series=[0+a/2:a:xmax-a/2]';
    mlt_subplt(PI_mean_x,x_series,102,1,1,1,"Preferential Index vs X position using running event","X(mm)","PI:v_x/sp",color,name,"line",'ii_num',ii_num,'ebar',1);
end
save_all_figures(outdir);
close all
% 1-5) Stat test for all the PI
group=vertcat(group{:});
[p_mwu,group_comp,p_k]=non_parametric_test(PI_larva,group,PI_size);
if length(PI_size)>2
    fig=gcf;
    n=fig.Number;
    outdir1=fullfile(outdir,"Stat_result");
    change_fig_prop(n,'title',"Kruskal Wallis test of the preferential index of different conditions",'xtitle',"Groups",'ytitle',"PI(v_x/sp)", ...
        'name',"Kruskal_Wallis_test_groups_boxplot_run_event");
    change_fig_prop(n-1,'name',"Kruskal_Wallis_test_groups_table_run_event");
    %save figures
    if ~isfolder(outdir1)
        mkdir(outdir1);
    end
    save_all_figures(outdir1);
end
%% save data

filename=fullfile(outdir1,'p_value.mat');

if isfile(filename)
    delete(filename);
end
save(filename,'p_k*','group','PI*','p_mwu*','group_comp*');
clear filename
end