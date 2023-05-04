function []=reori_angle_4_quadrants(tracker_num,genotype,condition,filename,output_name,varargin)
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
hbin=45;% for binning data
input_cond=[];%it will have a same length as name, telling us whether the input is 'ctr' or 'exp'
row_comp=2; col_comp=2; % these two will be used for the plot to show the comparison of ctr vs exp group, this number should be based on how many comparison you have
tmax=900;% tbin is used for binning the turn event's t0s
xmax=225;
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
    %% get the data ready for binning
    pre=vertcat(dat.pre_deg{:});
    %     post=vertcat(dat.post_deg{:});
    turn_x=vertcat(dat.turn_x{:});
    %     turn_y=vertcat(dat.turn_y{:});
    reori_deg=vertcat(dat.reorient_deg_abs{:});
    t0s=vertcat(dat.t0s{:});
    if isa(t0s,'cell')
        t0s=vertcat(t0s{:});
    end 
    xbin2=xmax/2;
    hdseries=[-180:hbin:180]';
    % hd_series1=[-180+hbin/2:hbin:180-hbin/2]';
    %% 1-1) Compute the mean value of reorientation angle in each quadrant
 
    if ii==1
        deg=zeros(ii_num,l); deg_prox=zeros(ii_num,l); deg_dist=zeros(ii_num,l);
        sem=zeros(ii_num,l); sem_prox=zeros(ii_num,l); sem_dist=zeros(ii_num,l);
    end
    [~,deg_list_q3]=bin_data_mean(pre,hdseries,reori_deg);
    if l==4
        [deg_list]=assign_num_4quadrants(deg_list_q3);
    elseif l==3
        [deg_list]=assign_num_3quadrants(deg_list_q3);
    end
    %compute the mean reori angle in each quadrant
    for j=1:l
        deg(ii,j)=mean(deg_list{j,1},'omitnan');
        len=length(deg_list{j,1});
        sem(ii,j)=std(deg_list{j,1},'omitnan')/sqrt(len);
    end

    %% save the data for kruskal wallis test
    % 1) for testing if there is any difference between the 3q within each
    % condition 2) for tesring if there is any significant difference in
    % each q between different condition --str_3q,deg_3q
    if ii==1
        str_3q=cell(l,1);deg_3q=cell(l,1);size_3q=cell(l,1);
        str_cond=cell(ii_num,1);deg_cond=cell(ii_num,1);size_cond=cell(ii_num,1);
    end

    for j=1:l
        len=length(deg_list{j,1});
        deg_cond{ii,1}=vertcat(deg_list{j,1},deg_cond{ii,1});
        str_cond{ii,1}=vertcat(repmat(append(name(ii,1),"-",q_legends(j)),len,1),str_cond{ii,1});
        size_cond{ii,1}=vertcat(len,size_cond{ii,1});

        deg_3q{j,1}=vertcat(deg_list{j,1},deg_3q{j,1});
        str_3q{j,1}=vertcat(repmat(append(name(ii),"-",q_legends(j)),len,1),str_3q{j,1});
        size_3q{j,1}=vertcat(len,size_3q{j,1});
    end
    clear deg_list deg_list_q3
    %% 1-2)Compute the mean reorientation angle in different quadrants and postions (prox vs dist)
    [deg_list_q3_x2]=bin_data_mean2(pre,hdseries,turn_x,[0:xbin2:xmax]',reori_deg);
    if l==4
        [deg_list]=assign_num_4quadrants(deg_list_q3_x2);
    elseif l==3
        [deg_list]=assign_num_3quadrants(deg_list_q3_x2);
    end

    for j=1:l
        deg_prox(ii,j)=mean(deg_list{j,2},'omitnan');
        deg_dist(ii,j)=mean(deg_list{j,1},'omitnan');
        len_prox=length(deg_list{j,2});
        len_dist=length(deg_list{j,1});

        sem_prox(ii,j)=std(deg_list{j,2})/sqrt(len_prox);
        sem_dist(ii,j)=std(deg_list{j,1})/sqrt(len_dist);
    end
    %% save data for kw stat
    if ii==1
        str_prox_3q=cell(l,1);deg_prox_3q=cell(l,1);size_prox_3q=cell(l,1);
        str_prox_cond=cell(ii_num,1);deg_prox_cond=cell(ii_num,1);size_prox_cond=cell(ii_num,1);

        str_dist_3q=cell(l,1);deg_dist_3q=cell(l,1);size_dist_3q=cell(l,1);
        str_dist_cond=cell(ii_num,1);deg_dist_cond=cell(ii_num,1);size_dist_cond=cell(ii_num,1);
    end

    for j=1:l
        len=length(deg_list{j,2});

        deg_prox_cond{ii,1}=vertcat(deg_list{j,2},deg_prox_cond{ii,1});
        str_prox_cond{ii,1}=vertcat(repmat(append(name(ii),"-",q_legends(j)),len,1),str_prox_cond{ii,1});
        size_prox_cond{ii,1}=vertcat(len,size_prox_cond{ii,1});

        deg_prox_3q{j,1}=vertcat(deg_list{j,2},deg_prox_3q{j,1});
        str_prox_3q{j,1}=vertcat(repmat(append(name(ii),"-",q_legends(j)),len,1),str_prox_3q{j,1});
        size_prox_3q{j,1}=vertcat(len,size_prox_3q{j,1});

        len=length(deg_list{j,1});

        deg_dist_cond{ii,1}=vertcat(deg_list{j,1},deg_dist_cond{ii,1});
        str_dist_cond{ii,1}=vertcat(repmat(append(name(ii),"-",q_legends(j)),len,1),str_dist_cond{ii,1});
        size_dist_cond{ii,1}=vertcat(len,size_dist_cond{ii,1});

        deg_dist_3q{j,1}=vertcat(deg_list{j,1},deg_dist_3q{j,1});
        str_dist_3q{j,1}=vertcat(repmat(append(name(ii),"-",q_legends(j)),len,1),str_dist_3q{j,1});
        size_dist_3q{j,1}=vertcat(len,size_dist_3q{j,1});
    end
    clear deg_list deg_list_q3_x2
    %% 4) clear variables
    clearvars -except cond row* col* ii* *bin* *max name color outdir output_name input_cond *series* deg* t l sem* str* size* q_legends
end
%% 1-3) plot the mean reori angle vs different quadrants
data={deg',sem'};
data_prox={deg_prox',sem_prox'};
data_dist={deg_dist',sem_dist'};
mlt_subplt(data,[],100,1,1,1,append("Mean_reori_angle_bar_",num2str(output_name)),"Orientation ","Mean reorientation angle(deg)",color,q_legends,"bar",'ii_num',ii_num,'ebar',1,'legends',name,'ylim',[35 70],'bar_color',color)
mlt_subplt(data_prox,[],101,1,1,1,append("Mean_reori_angle_bar_",num2str(output_name),"_proximate"),"Orientation ","Mean reorientation angle(deg, prox)",color,q_legends,"bar",'ii_num',ii_num,'ebar',1,'legends',name,'ylim',[35 70],'bar_color',color)
mlt_subplt(data_dist,[],102,1,1,1,append("Mean_reori_angle_bar_",num2str(output_name),"_distant"),"Orientation ","Mean reorientation angle(deg, dist)",color,q_legends,"bar",'ii_num',ii_num,'ebar',1,'legends',name,'ylim',[35 70],'bar_color',color)
%% SAVE FIGURES
outdir1=fullfile(outdir,"reori angle",append(num2str(l)," quadrants"));
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;
%% 1-4) perform stat test on the data
p_mwu_cond=cell(ii_num,1); group_comp_cond=cell(ii_num,1);p_k_cond=cell(ii_num,1);
p_mwu_3q=cell(l,1); group_comp_3q=cell(l,1);p_k_3q=cell(l,1);

p_mwu_prox_cond=cell(ii_num,1); group_prox_comp_cond=cell(ii_num,1);p_k_prox_cond=cell(ii_num,1);
p_mwu_prox_3q=cell(l,1); group_prox_comp_3q=cell(l,1);p_k_prox_3q=cell(l,1);

p_mwu_dist_cond=cell(ii_num,1); group_dist_comp_cond=cell(ii_num,1);p_k_dist_cond=cell(ii_num,1);
p_mwu_dist_3q=cell(l,1); group_dist_comp_3q=cell(l,1);p_k_dist_3q=cell(l,1);
for i=1:ii_num
    [p_mwu_cond{i,1},group_comp_cond{i,1},p_k_cond{i,1}]=non_parametric_test(deg_cond{i,1},str_cond{i,1},size_cond{i,1},0);

%     fig=gcf;
%     n=fig.Number;
%     change_fig_prop(n,'title',append("Kruskal Wallis test of the reorientation angle-",name(i,1)),...
%         'xtitle',"Groups",'ytitle',"Reorientation angle (deg)",'name',append("KW_test_boxplot_",name(i,1)));
%     change_fig_prop(n-1,'name',append("KW_test_",name(i,1),"_table"));
% 
    %proximate region
    [p_mwu_prox_cond{i,1},group_prox_comp_cond{i,1},p_k_prox_cond{i,1}]=non_parametric_test(deg_prox_cond{i,1},str_prox_cond{i,1},size_prox_cond{i,1},0);

%     fig=gcf;
%     n=fig.Number;
%     change_fig_prop(n,'title',append("Kruskal Wallis test of the reorientation angle-Proximal-",name(i,1)),...
%         'xtitle',"Groups",'ytitle',"Reorientation angle (deg)",'name',append("KW_test_boxplot_proximate_",name(i,1)));
%     change_fig_prop(n-1,'name',append("KW_test_proximate_",name(i,1),"_table"));

    %distant region
    [p_mwu_dist_cond{i,1},group_dist_comp_cond{i,1},p_k_dist_cond{i,1}]=non_parametric_test(deg_dist_cond{i,1},str_dist_cond{i,1},size_dist_cond{i,1},0);

    fig=gcf;
%     n=fig.Number;
%     change_fig_prop(n,'title',append("Kruskal Wallis test of the reorientation angle-Distant-",name(i,1)),...
%         'xtitle',"Groups",'ytitle',"Reorientation angle (deg)",'name',append("KW_test_boxplot_distant_",name(i,1)));
%     change_fig_prop(n-1,'name',append("KW_test_distant_",name(i,1),"_table"));
end

for i=1:l
    [p_mwu_3q{i,1},group_comp_3q{i,1},p_k_3q{i,1}]=non_parametric_test(deg_3q{i,1},str_3q{i,1},size_3q{i,1},0);

%     fig=gcf;
%     n=fig.Number;
%     change_fig_prop(n,'title',append("Kruskal Wallis test of the reorientation angle-",q_legends(i,1)),...
%         'xtitle',"Groups",'ytitle',"Reorientation angle (deg)",'name',append("KW_test_boxplot_",q_legends(i,1)));
%     change_fig_prop(n-1,'name',append("KW_test_",q_legends(i,1),"_table"));
    %peximate region
    [p_mwu_prox_3q{i,1},group_prox_comp_3q{i,1},p_k_prox_3q{i,1}]=non_parametric_test(deg_prox_3q{i,1},str_prox_3q{i,1},size_prox_3q{i,1},0);

%     fig=gcf;
%     n=fig.Number;
%     change_fig_prop(n,'title',append("Kruskal Wallis test of the reorientation angle-Proximal-",q_legends(i,1)),...
%         'xtitle',"Groups",'ytitle',"Reorientation angle (deg)",'name',append("KW_test_boxplot_proximate_",q_legends(i,1)));
%     change_fig_prop(n-1,'name',append("KW_test_proximate_",q_legends(i,1),"_table"));
    %distant region
    [p_mwu_dist_3q{i,1},group_dist_comp_3q{i,1},p_k_dist_3q{i,1}]=non_parametric_test(deg_dist_3q{i,1},str_dist_3q{i,1},size_dist_3q{i,1},0);

%     fig=gcf;
%     n=fig.Number;
%     change_fig_prop(n,'title',append("Kruskal Wallis test of the reorientation angle-Distant-",q_legends(i,1)),...
%         'xtitle',"Groups",'ytitle',"Reorientation angle (deg)",'name',append("KW_test_boxplot_distant_",q_legends(i,1)));
%     change_fig_prop(n-1,'name',append("KW_test_distant_",q_legends(i,1),"_table"));
end

p_mwu_cond=vertcat(p_mwu_cond{:}); group_comp_cond=vertcat(group_comp_cond{:});
p_mwu_prox_cond=vertcat(p_mwu_prox_cond{:}); group_prox_comp_cond=vertcat(group_prox_comp_cond{:});
p_mwu_dist_cond=vertcat(p_mwu_dist_cond{:}); group_dist_comp_cond=vertcat(group_dist_comp_cond{:});

p_mwu_3q=vertcat(p_mwu_3q{:}); group_comp_3q=vertcat(group_comp_3q{:});
p_mwu_prox_3q=vertcat(p_mwu_prox_3q{:}); group_prox_comp_3q=vertcat(group_prox_comp_3q{:});
p_mwu_dist_3q=vertcat(p_mwu_dist_3q{:}); group_dist_comp_3q=vertcat(group_dist_comp_3q{:});
%% SAVE FIGURES
outdir1=fullfile(outdir,"reori angle",append(num2str(l)," quadrants"),"stat");
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;
%% save data


if length(p_mwu_cond)==length(p_mwu_prox_cond) && length(p_mwu_prox_cond)==length(p_mwu_dist_cond)
    filename4=fullfile(outdir1,'MWU_test_each_condition.csv');
    T=table(group_comp_cond,p_mwu_cond,group_prox_comp_cond,p_mwu_prox_cond,group_dist_comp_cond,p_mwu_dist_cond);
    writetable(T,filename4);
else
    filename4=fullfile(outdir1,'MWU_test_each_condition.mat');
    save(filename4,'p_k*_cond','p_mwu*_cond','group_comp*_cond');
end


filename1=fullfile(outdir1,'KW_test_each_condition.txt');
T_kw=table(p_k_cond,p_k_prox_cond,p_k_dist_cond);
if isfile(filename1)
    delete(filename1);
end
writetable(T_kw,filename1);


 

if length(p_mwu_3q)==length(p_mwu_prox_3q) && length(p_mwu_prox_3q)==length(p_mwu_dist_3q)
    filename2=fullfile(outdir1,'MWU_test_each_quadrant.csv');
    T_3q=table(group_comp_3q,p_mwu_3q,group_prox_comp_3q,p_mwu_prox_3q,group_dist_comp_3q,p_mwu_dist_3q);
    writetable(T_3q,filename2);
else 
    filename2=fullfile(outdir1,'MWU_test_each_quadrant.mat');
    save(filename2,'p_k*_3q','p_mwu*_3q','group_*_comp_3q');
end 

    
T_kw_3q=table(p_k_3q,p_k_prox_3q,p_k_dist_3q);

filename3=fullfile(outdir1,'KW_test_each_quadrant.txt');
if isfile(filename3)
    delete(filename3);
end
writetable(T_kw_3q,filename3);

clear
end