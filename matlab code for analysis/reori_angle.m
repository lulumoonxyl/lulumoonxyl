function []=reori_angle(tracker_num,genotype,condition,filename,output_name,varargin)
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
outdir=fullfile("/project/6010970/screen/olfactory_output/Properties_turn_events/",output_name);
ii_num=length(cond);
%% set the properties for plots
row=3;col=3;
color=[0.75,0.75,0.75;0.627,0.82,1;0.659,0.659,0.659;0.871,0.416,0.451];
hbin=5;% for binning data-->frequency of reorientation angle
input_cond=[];%it will have a same length as name, telling us whether the input is 'ctr' or 'exp'
row_comp=2; col_comp=2; % these two will be used for the plot to show the comparison of ctr vs exp group, this number should be based on how many comparison you have
hbin2=10;%use for the mean of reorientation angle
xmax=225;
for i=1:2:length(varargin)
    if strcmp ('row',varargin{i})
        row=varargin{i+1};
    elseif strcmp ('col',varargin{i})
        col=varargin{i+1};
    elseif strcmp ('hbin2',varargin{i})
        hbin2=varargin{i+1};
    elseif strcmp ('hbin',varargin{i})
        hbin=varargin{i+1};
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

    end
end

for ii=1:ii_num
    load(cond(ii));
    
    hd_series=[-180:hbin:180]';
    hd_series4=[-180:hbin2:180]';
    %% get the data ready for binning
    pre=vertcat(dat.pre_deg{:});
    %     post=vertcat(dat.post_deg{:});
    turn_x=vertcat(dat.turn_x{:});
    %     turn_y=vertcat(dat.turn_y{:});
    reori_deg=vertcat(dat.reorient_deg_abs{:});
    xbin2=xmax/2; %separate the x position into two regions
    %% 3-1) save this for violin  plot
    if ii==1
        size=zeros(ii_num,1); reori_vio=cell(ii_num,1); str=[];
        size_prox=zeros(ii_num,1); reori_vio_prox=cell(ii_num,1); str_prox=[];
        size_dist=zeros(ii_num,1); reori_vio_dist=cell(ii_num,1); str_dist=[];
    end
    reori_vio{ii,1}=reori_deg;
    size(ii,1)=length(reori_deg);
    str=vertcat( str,repmat({name(ii)},size(ii,1),1));

    idx_prox=find(turn_x>=0&turn_x<xmax/2);
    reori_vio_prox{ii,1}=reori_deg(idx_prox);
    size_prox(ii,1)=length(idx_prox);
    str_prox=vertcat( str_prox,repmat({name(ii)},size_prox(ii,1),1));

    idx_dist=find(turn_x>=xmax/2&turn_x<=xmax);
    reori_vio_dist{ii,1}=reori_deg(idx_dist);
    size_dist(ii,1)=length(idx_dist);
    str_dist=vertcat( str_dist,repmat({name(ii)},size_dist(ii,1),1));
    %% 1-1) Frequency of reorientation angle: Bin the reorientation angle from [-180 180]-->for polar plot later
    c=bin_data_count(reori_deg,hd_series);
    if ii==1
        f=cell(ii_num,1); f_all=cell(ii_num,1);
        f_prox=cell(ii_num,1); f_dist=cell(ii_num,1);
        f_prox_polar=cell(ii_num,1); f_dist_polar=cell(ii_num,1);
    end
    f{ii,1}(:,1)=c/sum(c);
    f{ii,1}(:,2)=1.96.*(sqrt(c)./sum(c));
    clear c

    %%  bin from [0 180]-->use for bar chart
    hd_series3=[0:hbin:180]';
    c=bin_data_count(reori_deg,hd_series3);
    f1(:,1)=c./sum(c);
    f1(:,2)=1.96.*(sqrt(c)./sum(c));
    hd_series1=[0+hbin/2:hbin:180-hbin/2]';
    f_all{ii,1}=f1;
    clear c

    %% Compute the frequency in proximate and distant for bar plot
    % polar plot
    c_x2=bin_data_count2(reori_deg,hd_series3,turn_x,[0:xbin2:xmax]');
    f_x2=c_x2./sum(c_x2);
    f_sem2=1.96.*(sqrt(c_x2)./sum(c_x2));
    f_prox{ii,1}(:,1)=f_x2(:,2); f_prox{ii,1}(:,2)=f_sem2(:,2);
    f_dist{ii,1}(:,1)=f_x2(:,1); f_dist{ii,1}(:,2)=f_sem2(:,1);
    clear f_x2 f_sem2 c_x2

    %% 1-3) Also compute the frequency of reorientation angle based on x for
    % polar plot
    c_x2=bin_data_count2(reori_deg,hd_series,turn_x,[0:xbin2:xmax]');
    f_x2=c_x2./sum(c_x2);
    f_sem2=1.96.*(sqrt(c_x2)./sum(c_x2));
    f_prox_polar{ii,1}(:,1)=f_x2(:,2); f_prox_polar{ii,1}(:,2)=f_sem2(:,2);
    f_dist_polar{ii,1}(:,1)=f_x2(:,1); f_dist_polar{ii,1}(:,2)=f_sem2(:,1);
    clear f_x2 f_sem2 c_x2
    %% 1-4) plot the frequency vs angle--bar plot
    mlt_subplt(f1,hd_series1,100,row,col,ii,append("Frequency of reorientation angle ",num2str(output_name)),"Reorientation (deg)","Frequency",color,name(ii),"bar",'ii_num',ii_num,'ebar',1);

    clear f1
    %% 2-1) Mean reorientation angle across HD and in different x positions (proximate vs distant and all positions)
    if ii==1
        deg=cell(ii_num,1); deg_prox=cell(ii_num,1);deg_dist=cell(ii_num,1);
    end

    [~,deg_x2,sem_x2]=bin_data_mean2(pre,hd_series4,turn_x,[0:xbin2:xmax]',reori_deg);
    deg_prox{ii,1}(:,1)=deg_x2(:,2); deg_prox{ii,1}(:,2)=sem_x2(:,2);
    deg_dist{ii,1}(:,1)=deg_x2(:,1); deg_dist{ii,1}(:,2)=sem_x2(:,1);

    [deg{ii,1}]=bin_data_mean(pre,hd_series4,reori_deg);

    %% 4) clear variables
    clearvars -except cond row* col* ii* *bin* *max name color outdir f* output_name input_cond *series* f* deg* str* reori_vio* size*
end

%% 1-5) plot the distribution of the reorientation angle --polar plot
mlt_subplt(f,[],101,row,col,ii,append("Frequency of reorientation angle_polar plot ",num2str(output_name))," "," ",color,name,"polar plot",'bins',hbin,'max',0.08,'tlim',[0 180],'ebar',1);
%% 2-2) plot the mean reorientation angle for different x positions
% mlt_subplt(deg,[],104,row,col,ii,append("Mean reorientation angle_polar plot ",num2str(output_name))," "," ",color,name,"polar plot",'bins',hbin2,'max',75,'min',30,'ebar',1);
% mlt_subplt(deg_prox,[],105,row,col,ii,append("Mean reorientation angle_Proximate region_polar plot ",num2str(output_name))," "," ",color,name,"polar plot",'bins',hbin2,'max',75,'min',30,'ebar',1);
% mlt_subplt(deg_dist,[],106,row,col,ii,append("Mean reorientation angle_Distant region_polar plot ",num2str(output_name))," "," ",color,name,"polar plot",'bins',hbin2,'max',75,'min',30,'ebar',1);

%% compare the ctr and exp in the same bar chart or polar plot
if ~isempty(input_cond)
    idx=find(contains(input_cond,"ctr"));
    [ctr_f,exp_f,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(f_all,idx,name,color,'input_cond',input_cond);
    [ctr_f_prox,exp_f_prox]=split_data_ctr_exp(f_prox,idx,name,color,'input_cond',input_cond);
    [ctr_f_dist,exp_f_dist]=split_data_ctr_exp(f_dist,idx,name,color,'input_cond',input_cond);
    [ctr_f_polar,exp_f_polar]=split_data_ctr_exp(f,idx,name,color,'input_cond',input_cond);

    [ctr_deg_polar1,exp_deg_polar1]=split_data_ctr_exp(deg,idx,name,color,'input_cond',input_cond);
    [ctr_deg_prox_polar,exp_deg_prox_polar]=split_data_ctr_exp(deg_prox,idx,name,color,'input_cond',input_cond);
    [ctr_deg_dist_polar,exp_deg_dist_polar]=split_data_ctr_exp(deg_dist,idx,name,color,'input_cond',input_cond);

    for i=1:length(exp_f)
        % get the color for the ctr and exp groups
        legends=[ctr_name(i,1);exp_name(i,1)];
        color1=[ctr_color(i,:);exp_color(i,:)];

        %% 1-6) frequency of reorientation deg for the ctr group--bar plot
        data_f={ctr_f{i,1};exp_f{i,1}};
        mlt_subplt(data_f,hd_series1,102,row_comp,col_comp,i,append("Frequency of reorientation angle_bar plot ",num2str(output_name),", ctr vs exp"),"Reorientation (deg)","Frequency",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1)),"overlapped bar",'legend',legends,'xlim',[0 180],'ylim',[0 0.08]);

        data_f_prox={ctr_f_prox{i,1};exp_f_prox{i,1}};
        mlt_subplt(data_f_prox,hd_series1,110,row_comp,col_comp,i,append("Frequency of reorientation angle_Proximate region_bar plot ",num2str(output_name),", ctr vs exp"),"Reorientation (deg)","Frequency",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),'_ _p_r_o_x_i'),"overlapped bar",'legend',legends,'xlim',[0 180],'ylim',[0 0.08]);

        data_f_dist={ctr_f_dist{i,1};exp_f_dist{i,1}};
        mlt_subplt(data_f_dist,hd_series1,111,row_comp,col_comp,i,append("Frequency of reorientation angle_Distant region_bar plot ",num2str(output_name),", ctr vs exp"),"Reorientation (deg)","Frequency",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),'_ _d_i_s_t'),"overlapped bar",'legend',legends,'xlim',[0 180],'ylim',[0 0.08]);

        %% 1-7) & 2-3) frequency and mean of reorientation deg for the ctr group--polar plot
        ctr_reori_f=ctr_f_polar{i,1};
        exp_reori_f=exp_f_polar{i,1};

        ctr_deg_prox_mean=ctr_deg_prox_polar{i,1};
        exp_deg_prox_mean=exp_deg_prox_polar{i,1};

        ctr_deg_dist_mean=ctr_deg_dist_polar{i,1};
        exp_deg_dist_mean=exp_deg_dist_polar{i,1};

        ctr_deg_mean=ctr_deg_polar1{i,1};
        exp_deg_mean=exp_deg_polar1{i,1};
        %% frequency of reorientaion angle
        hbin1=hbin/2;
        hd_series2=[-180+hbin1/2:hbin1:180-hbin1/2]';
        l1=length(hd_series2); l2=length(ctr_reori_f);
        em_arry=zeros(l1-l2,1);

        for j=1:2
            data_f1{1,1}(:,j)=reshape([ctr_reori_f(:,j)';em_arry'],1,[])';
            data_f1{2,1}(:,j)=reshape([em_arry';exp_reori_f(:,j)'],1,[])';
        end

        %% Mean of reorientation angle
        hbin3=hbin2/2;
        hd_series5=[-180+hbin3/2:hbin3:180-hbin3/2]';
        l1=length(hd_series5);
        l2=length(ctr_deg_mean);
        em_arry1(:,1)=ones(l1-l2,1).*20; % 45 is the 'min' input to the line 197-199
        em_arry1(:,2)=zeros(l1-l2,1); % 45 is the 'min' input to the line 197-199
        for j=1:2
            data_deg_mean{1,1}(:,j)=reshape([ctr_deg_mean(:,j)';em_arry1(:,j)'],1,[])';
            data_deg_mean{2,1}(:,j)=reshape([em_arry1(:,j)';exp_deg_mean(:,j)'],1,[])';

            data_deg_prox_mean{1,1}(:,j)=reshape([ctr_deg_prox_mean(:,j)';em_arry1(:,j)'],1,[])';
            data_deg_prox_mean{2,1}(:,j)=reshape([em_arry1(:,j)';exp_deg_prox_mean(:,j)'],1,[])';

            data_deg_dist_mean{1,1}(:,j)=reshape([ctr_deg_dist_mean(:,j)';em_arry1(:,j)'],1,[])';
            data_deg_dist_mean{2,1}(:,j)=reshape([em_arry1(:,j)';exp_deg_dist_mean(:,j)'],1,[])';
        end
        %% plot the frequency of reorientation angle
        mlt_subplt(data_f1,[],103,row_comp,col_comp,i,append("Frequency of reorientation angle_polar plot ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1)),"overlapped polar plot",'bins',hbin1,'tlim',[0 180],'max',0.08,'ebar',1);
        %% plot the mean of reorientation angle
        mlt_subplt(data_deg_mean,[],107,row_comp,col_comp,i,append("Mean reorientation angle_polar plot ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1)),"overlapped polar plot",'bins',hbin3,'max',85,'min',20,'ebar',1);
        mlt_subplt(data_deg_prox_mean,[],108,row_comp,col_comp,i,append("Mean reorientation angle_Proximate region_polar plot ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),'_ _p_r_o_x'),"overlapped polar plot",'bins',hbin3,'max',85,'min',20,'ebar',1);
        mlt_subplt(data_deg_dist_mean,[],109,row_comp,col_comp,i,append("Mean reorientation angl_Distant region_polar plot ",num2str(output_name),", ctr vs exp")," "," ",color1,append(ctr_name(i,1)," vs. ",exp_name(i,1),'_ _d_i_s_t'),"overlapped polar plot",'bins',hbin3,'max',85,'min',20,'ebar',1);

    end
end

%% 1-8) plot violin plot to show the distribution of reorientation angle
for i=1:length(name)
    name1{1,i}=name(i,:);
end 
figure(112)
[~,~,reori_mean,reori_med]=violin(reori_vio','xlabel',name1,'facecolor',color);
change_fig_prop(112,'title',"All positions",'xtitle',"Groups",'ytitle',"Reorientation angle(deg)", ...
    'name',"Violinplot reorientation angle for different group_all positions","ylim",[-15 250]);
figure(113)
[~,~,reori_mean_prox,reori_med_prox]=violin(reori_vio_prox','xlabel',name1,'facecolor',color);
change_fig_prop(113,'title',"Proximate region",'xtitle',"Groups",'ytitle',"Reorientation angle(deg)", ...
    'name',"Violinplot reorientation angle for different group_proximate region","ylim",[-15 250]);
figure(114)
[~,~,reori_mean_dist,reori_med_dist]=violin(reori_vio_dist','xlabel',name1,'facecolor',color);
change_fig_prop(114,'title',"Distant region",'xtitle',"Groups",'ytitle',"Reorientation angle(deg)", ...
    'name',"Violinplot for reorientation angle for different group_distant region","ylim",[-15 250]);
%% SAVE FIGURES
outdir_fig=fullfile(outdir,"reori angle");
if ~isfolder(outdir_fig)
    mkdir(outdir_fig);
end
save_all_figures(outdir_fig);
close all;
%% stat test
str=vertcat(str{:});
reori=vertcat(reori_vio{:});
[p_mwu,group_comp,p_k]=non_parametric_test(reori,str,size,0);
% fig=gcf;
% n=fig.Number;
outdir1=fullfile(outdir_fig,"Reori_angle_result");
% change_fig_prop(n,'title',"Kruskal Wallis test of the reorientation angle for different conditions",'xtitle',"Groups",'ytitle',"Reorientation angle(deg)", ...
%     'name',"Reori_angle_Kruskal_Wallis_test_groups_boxplot");
% change_fig_prop(n-1,'name',"Reori_angle_Kruskal_Wallis_test_groups_table");
% %save figures
% if ~isfolder(outdir1)
%     mkdir(outdir1);
% end
% save_all_figures(outdir1);
% close all
% for proximate region
str_prox=vertcat(str_prox{:});
reori_prox=vertcat(reori_vio_prox{:});
[p_mwu_prox,group_comp_prox,p_k_prox]=non_parametric_test(reori_prox,str_prox,size_prox,0);
% fig=gcf;
% n=fig.Number;
% 
% change_fig_prop(n,'title',"Proximate: Kruskal Wallis test of the reorientation angle for different conditions",'xtitle',"Groups",'ytitle',"Reorientation angle(deg)", ...
%     'name',"Reori_angle_Proximate_Kruskal_Wallis_test_groups_boxplot");
% change_fig_prop(n-1,'name',"Reori_angle_Proximate_region_Kruskal_Wallis_test_groups_table");
% %save figures
% if ~isfolder(outdir1)
%     mkdir(outdir1);
% end
% save_all_figures(outdir1);
% close all
% for distant region
str_dist=vertcat(str_dist{:});
reori_dist=vertcat(reori_vio_dist{:});
[p_mwu_dist,group_comp_dist,p_k_dist]=non_parametric_test(reori_dist,str_dist,size_dist,0);
% fig=gcf;
% n=fig.Number;
% 
% change_fig_prop(n,'title',"Distant: Kruskal Wallis test of the reorientation angle for different conditions",'xtitle',"Groups",'ytitle',"Reorientation angle(deg)", ...
%     'name',"Reori_angle_Distant_Kruskal_Wallis_test_groups_boxplot");
% change_fig_prop(n-1,'name',"Reori_angle_Distant_region_Kruskal_Wallis_test_groups_table");
%save figures
if ~isfolder(outdir1)
    mkdir(outdir1);
end
% save_all_figures(outdir1);
% close all
%save p value 
filename=fullfile(outdir1,'reori_ang.mat');

if isfile(filename)
    delete(filename);
end

save(filename,'p_k*','p_mwu*','group_comp*','reori_*');
clear filename

data=table(group_comp,p_mwu,group_comp_prox,p_mwu_prox,group_comp_dist,p_mwu_dist);
filename=fullfile(outdir1,'MWU_test.xlsx');
if isfile(filename)
    delete(filename);
end
writetable(data,filename);
close all
clear
end