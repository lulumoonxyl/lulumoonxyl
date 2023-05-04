function []=run_prop(tracker_num,genotype,condition,filename,output_name,varargin)
%% 1) Compute the frequency of running duration--bar chart--one is vtr vs exp, one is just barchart for each cond
%2) compute the distribution of absolute reorientation angle for all
%turning events--polar plot and bar chart --one for each cond or one for
%ctr vs exp

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
outdir=fullfile("/project/6010970/screen/olfactory_output/Properties_run_events",output_name);
ii_num=length(cond);

%% set the properties for plots
row=3;col=3;
color=[0.75,0.75,0.75;0.627,0.82,1;0.659,0.659,0.659;0.871,0.416,0.451];
rbin=2; % this bin is used for binning the running duration
tmax=900; % for binning data
input_cond=[];%it will have a same length as name, telling us whether the input is 'ctr' or 'exp'
row_comp=2; col_comp=2; % these two will be used for the plot to show the comparison of ctr vs exp group, this number should be based on how many comparison you have

for i=1:2:length(varargin)
    if strcmp ('row',varargin{i})
        row=varargin{i+1};
    elseif strcmp ('col',varargin{i})
        col=varargin{i+1};
    elseif strcmp ('rbin',varargin{i})
        rbin=varargin{i+1};
    elseif strcmp ('tmax',varargin{i})
        tmax=varargin{i+1};
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
    %% 1) bin the run duration for calculating the distribution later
    r_series=[0:rbin:tmax]';
    run_et=vertcat(dat.run_et{:});
    c=bin_data_count(run_et,r_series);
    f_rt(:,1)=c/sum(c);
    f_rt(:,2)=1.96.*(sqrt(c)./sum(c));% get the 95% CI using poisson approximation
    if ii==1
        f_rt_all=cell(ii_num,1);
    end 
    f_rt_all{ii,1}=f_rt;
    %% plot the distribution of running event duration
    r_series1=[0+rbin/2:rbin:tmax-rbin/2]';
    mlt_subplt(f_rt,r_series1,100,row,col,ii,append("Frequency of running duration ",num2str(output_name)),"Running duration (s)","Frequency",color,name(ii),"bar",'ii_num',ii_num,'xlim',[0 100],'ebar',1,'ylim',[0 0.35]);
    clear c

    %% 2) clear variables
    clearvars -except cond row* col* ii* *bin *max name color outdir f_deg* output_name f_rt_* input_cond
end


%% compare the ctr and exp in the same bar chart
if ~isempty(input_cond)
    idx=find(contains(input_cond,"ctr"));
    %% 3) plot the ctr and exp groups in the same subplot plot 
    %% running duration
    r_series1=[0+rbin/2:rbin:tmax-rbin/2]';
    [ctr,exp,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(f_rt_all,idx,name,color,'input_cond',input_cond);

    for i=1:length(exp)
        data={ctr{i,1};exp{i,1}};
        legends=[ctr_name(i,1);exp_name(i,1)];
        color1=[ctr_color(i,:);exp_color(i,:)];
        %% running duration
        mlt_subplt(data,r_series1,101,row_comp,col_comp,i,append("Frequency of running duration ",num2str(output_name),", ctr vs exp"),"Running duration (s)","Frequency",color1,append(ctr_name(i,1),"-",exp_name(i,1)),"overlapped bar",'legend',legends,'xlim',[0 100],'ylim',[0 0.3]);
    end
end

%% SAVE FIGURES
if ~isfolder(outdir)
    mkdir(outdir);
end
save_all_figures(outdir);
close all;
end