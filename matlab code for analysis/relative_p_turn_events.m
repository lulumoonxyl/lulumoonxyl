function relative_p_turn_events(tracker_num,genotype,condition,filename,output_name,varargin)
%condition should be a string array, genotype can also be a string array
%output_name is the output folder name 
w=1;
%     name=["10E-1GA","10E-2GA","10E-3GA","10^-^5EA","H2O"];
for i=1:length(genotype)
    for j=1:length(condition)
        cond(w,1)=fullfile("/project/6010970/screen/olfactory_output/JB_JAABA",tracker_num,genotype{i,1},condition{j,1},filename);
        name(w,1)=append(genotype{i,1},"@",condition{j,1});
        w=w+1;
    end
end

cond

% outdir="C:\Users\feihu\OneDrive - McGill University\matlab\Analysis based on larvae\odor\relative probability\updated analysis\turning frequency";
outdir=fullfile("/project/6010970/screen/olfactory_output/turn_freq",output_name);
ii_num=length(cond);
%set the properties for plots
row=3;col=3;
color=[0,0,0;0.659,0.659,0.659;0.255,0.643,1;0.439,0.733,1;0.627,0.82,1;0.871,0.416,0.451;1,0,0];
xmax=230;ymax=250;tmax=900;bin=10;tbin=30;hdbin=20;deg_l=60;
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
    elseif strcmp(varargin{i},'tbin')
        tbin=varargin{i+1};
    elseif strcmp(varargin{i},'bin')
        bin=varargin{i+1};
    elseif strcmp(varargin{i},'hdbin')
        hdbin=varargin{i+1};
    elseif strcmp(varargin{i},'name')
        clear name;
        name=string(varargin{i+1})
    end
end
name
xseries=[0:bin:xmax]'; x_series=[0+bin/2:bin:xmax-bin/2]';
yseries=[0:bin:ymax]'; y_series=[0+bin/2:bin:ymax-bin/2]';
tseries=[0:tbin:tmax]';  t_series=[0+tbin/2:tbin:tmax-tbin/2]';
hdseries=[-180:hdbin:180]'; hd_series=[-180+hdbin/2:hdbin:180-hdbin/2]';

for ii=1:ii_num
    load(cond(ii));
 
          % next reduce the turning frequency plot from 3D to 2D, therefore plot
    % only turning frquency vs HD only, x only, y only, and t only
    if ii==1
        pc_hd_all=cell(ii_num,1);pc_x_all=cell(ii_num,1);pc_y_all=cell(ii_num,1);pc_t_all=cell(ii_num,1);
    end
    c_hd=sum(c_x,2,'omitnan'); t_hd2=sum(t_x1,2,'omitnan');pc_hd_all{ii,1}=c_hd./t_hd2;
    c_x1=sum(c_x,1,'omitnan'); t_x2=sum(t_x1,1,'omitnan');pc_x_all{ii,1}=(c_x1./t_x2)';
    c_y1=sum(c_y,1,'omitnan'); t_y2=sum(t_y1,1,'omitnan');pc_y_all{ii,1}=(c_y1./t_y2)';
    c_t1=sum(c_t,1,'omitnan'); t_t2=sum(t_t1,1,'omitnan');pc_t_all{ii,1}=(c_t1./t_t2)';

    % next we will bin the x into three regions and check the turning
    % frequency vs HD
    xbin3=xmax/3;ybin3=ymax/3;tbin3=tmax/3;
    t_x3=bin_data_sum_t2(dat.orientation,dat.x,hdseries,[0:xbin3:xmax]',dat.et);
    t_y3=bin_data_sum_t2(dat.orientation,dat.y,hdseries,[0:ybin3:ymax]',dat.et);
    t_t3=bin_data_sum_t2(dat.orientation,dat.et,hdseries,[0:tbin3:tmax]',dat.et);

    c_x3=bin_data_count2(vertcat(dat.pre_deg{:}),hdseries,vertcat(dat.turn_x{:}),[0:xbin3:xmax]');
    pc_x3=c_x3./t_x3;
    c_y3=bin_data_count2(vertcat(dat.pre_deg{:}),hdseries,vertcat(dat.turn_y{:}),[0:ybin3:ymax]');
    pc_y3=c_y3./t_y3;
    c_t3=bin_data_count2(vertcat(dat.pre_deg{:}),hdseries,vertcat(t_ver{:}),[0:tbin3:tmax]');
    pc_t3=c_t3./t_t3;

    if ii==1
        pc_x_dis=cell(ii_num,1);pc_x_cent=cell(ii_num,1);pc_x_prox=cell(ii_num,1);
        pc_y_up=cell(ii_num,1);pc_y_cent=cell(ii_num,1);pc_y_low=cell(ii_num,1);
        pc_t300=cell(ii_num,1);pc_t600=cell(ii_num,1);pc_t900=cell(ii_num,1);
    end
    pc_x_dis{ii,1}=pc_x3(:,1);pc_x_cent{ii,1}=pc_x3(:,2);pc_x_prox{ii,1}=pc_x3(:,3);
    pc_y_low{ii,1}=pc_y3(:,1);pc_y_cent{ii,1}=pc_y3(:,2);pc_y_up{ii,1}=pc_y3(:,3);
    pc_t300{ii,1}=pc_t3(:,1);pc_t600{ii,1}=pc_t3(:,2);pc_t900{ii,1}=pc_t3(:,3);

    %% 6) Compute the reorientation angle in different HD and y/x/t bin
    [deg_alx,deg_x]=bin_data_mean2(vertcat(dat.pre_deg{:}),hdseries,vertcat(dat.turn_x{:}),xseries,vertcat(dat.reorient_deg_abs{:}));
    [deg_aly,deg_y]=bin_data_mean2(vertcat(dat.pre_deg{:}),hdseries,vertcat(dat.turn_y{:}),yseries,vertcat(dat.reorient_deg_abs{:}));
    [deg_alt,deg_t]=bin_data_mean2(vertcat(dat.pre_deg{:}),hdseries,vertcat(t_ver{:}),tseries,vertcat(dat.reorient_deg_abs{:}));

  
    % similar to the turn freq, we can reduce it from 3D to 2D; therefore
    % plot the reorient angle vs HD/X/Y/T
    if ii==1
        deg_x_all=cell(ii_num,1); deg_y_all=cell(ii_num,1); deg_t_all=cell(ii_num,1);deg_hd_all=cell(ii_num,1);
    end
    for i=1:height(deg_alx)
        deg={deg_alx{i,:}}';
        deg=vertcat(deg{:});
        deg_hd1(i,1)=mean(deg);
        deg_hd1(i,2)=std(deg)/sqrt(length(deg));
        clear deg
    end

    for i=1:width(deg_alx)
        deg={deg_alx{:,i}}';
        deg=vertcat(deg{:});
        deg_x1(i,1)=mean(deg);
        deg_x1(i,2)=std(deg)/sqrt(length(deg));
        clear deg
    end
    for i=1:width(deg_aly)
        deg={deg_aly{:,i}}';
        deg=vertcat(deg{:});
        deg_y1(i,1)=mean(deg);
        deg_y1(i,2)=std(deg)/sqrt(length(deg));
        clear deg
    end
    for i=1:width(deg_alt)
        deg={deg_alt{:,i}}';
        deg=vertcat(deg{:});
        deg_t1(i,1)=mean(deg);
        deg_t1(i,2)=std(deg)/sqrt(length(deg));
        clear deg
    end
    deg_x_all{ii,1}=deg_x1;deg_y_all{ii,1}=deg_y1;deg_t_all{ii,1}=deg_t1;deg_hd_all{ii,1}=deg_hd1;
    % we can also separate x,y,t into three regions
    xbin3=xmax/3;ybin3=ymax/3;tbin3=tmax/3;
    t_ver=vertcat(dat.t0s{:});
    [~,deg_x3,sem_x3]=bin_data_mean2(vertcat(dat.pre_deg{:}),hdseries,vertcat(dat.turn_x{:}),[0:xbin3:xmax]',vertcat(dat.reorient_deg_abs{:}));
    [~,deg_y3,sem_y3]=bin_data_mean2(vertcat(dat.pre_deg{:}),hdseries,vertcat(dat.turn_y{:}),[0:ybin3:ymax]',vertcat(dat.reorient_deg_abs{:}));
    [~,deg_t3,sem_t3]=bin_data_mean2(vertcat(dat.pre_deg{:}),hdseries,vertcat(t_ver{:}),[0:tbin3:tmax]',vertcat(dat.reorient_deg_abs{:}));
    if ii==1
        deg_x_dis=cell(ii_num,1);deg_x_cent=cell(ii_num,1);deg_x_prox=cell(ii_num,1);
        deg_y_up=cell(ii_num,1);deg_y_cent=cell(ii_num,1);deg_y_low=cell(ii_num,1);
        deg_t300=cell(ii_num,1);deg_t600=cell(ii_num,1);deg_t900=cell(ii_num,1);
    end
    %save the results in cell arrys
    deg_x_dis{ii,1}(:,1)=deg_x3(:,1);deg_x_cent{ii,1}(:,1)=deg_x3(:,2);deg_x_prox{ii,1}(:,1)=deg_x3(:,3);
    deg_y_low{ii,1}(:,1)=deg_y3(:,1);deg_y_cent{ii,1}(:,1)=deg_y3(:,2);deg_y_up{ii,1}(:,1)=deg_y3(:,3);
    deg_t300{ii,1}(:,1)=deg_t3(:,1);deg_t600{ii,1}(:,1)=deg_t3(:,2);deg_t900{ii,1}(:,1)=deg_t3(:,3);
    deg_x_dis{ii,1}(:,2)=sem_x3(:,1);deg_x_cent{ii,1}(:,2)=sem_x3(:,2);deg_x_prox{ii,1}(:,2)=sem_x3(:,3);
    deg_y_low{ii,1}(:,2)=sem_y3(:,1);deg_y_cent{ii,1}(:,2)=sem_y3(:,2);deg_y_up{ii,1}(:,2)=sem_y3(:,3);
    deg_t300{ii,1}(:,2)=sem_t3(:,1);deg_t600{ii,1}(:,2)=sem_t3(:,2);deg_t900{ii,1}(:,2)=sem_t3(:,3);

    clear deg_*1 deg_al* deg_x deg_y deg_t
    % 7) compute the acceptance rate, toward odor counts as accept, avoid is non-acceptance
    pre=vertcat(dat.pre_deg{:}); post=vertcat(dat.post_deg{:}); turn_x=vertcat(dat.turn_x{:});turn_y=vertcat(dat.turn_y{:});
    ang=vertcat(dat.reorient_deg_abs{:});t0s1=vertcat(dat.t0s{:});t0s=vertcat(t0s1{:});
    q_post=check_qua_turn(post);
    q_pre=check_qua_turn(pre);
    q_ver=find(q_pre==2|q_pre==4);

    idx=find((q_pre==2|q_pre==4)&abs(pre)>abs(post));
    ac(ii,1)=length(idx)/length(q_ver);
    ac(ii,2)=1.96*sqrt(ac(ii,1)*(1-ac(ii,1))/length(q_ver));
    ver_turn_all(ii,1)=length(q_ver);
    ver_turn_ac(ii,1)=length(idx);
    %get the acceptance rate for three region
    xbin3=xmax/3;ybin3=ymax/3;tbin3=tmax/3;
    ac_x=bin_data_count(turn_x(idx),[0:xbin3:xmax]);
    ac_x_tt=bin_data_count(turn_x(q_ver),[0:xbin3:xmax]);
    ac_x3(ii,:)=ac_x./ac_x_tt;
    ac_x3_sem(ii,:)=1.96*sqrt(ac_x3(ii,:).*(1-ac_x3(ii,:))./ac_x_tt');
    %get the acceptance rate for large reorientation event (>60) and three
    %regions
    q_verl=find((q_pre==2|q_pre==4)&ang>deg_l);
    idxl=find((q_pre==2|q_pre==4)&abs(pre)>abs(post)&ang>deg_l);
    acl(ii,1)=length(idxl)/length(q_verl);
    acl(ii,2)=1.96*sqrt(acl(ii,1)*(1-acl(ii,1))/length(q_verl));

    ac_xl=bin_data_count(turn_x(idxl),[0:xbin3:xmax]);
    ac_x_ttl=bin_data_count(turn_x(q_verl),[0:xbin3:xmax]);
    ac_x3l(ii,:)=ac_xl./ac_x_ttl;
    ac_x3_seml(ii,:)=1.96*sqrt(ac_x3l(ii,:).*(1-ac_x3l(ii,:))./ac_x_ttl');
    ver_turn_all_l(ii,1)=length(q_verl);
    ver_turn_ac_l(ii,1)=length(idxl);
    clear ac_x ac_x_tt idx idxl ac_xl ac_x_ttl
    %% 8) calculate count of turning event,turning frequency
    idx=find(ang>=deg_l);
    c_ang_x=bin_data_count2(pre(idx),hdseries,turn_x(idx),xseries);
    c_ang_y=bin_data_count2(pre(idx),hdseries,turn_y(idx),yseries);
    c_ang_t=bin_data_count2(pre(idx),hdseries,t0s(idx),tseries);

    pc_ang_x=c_ang_x./t_x1;
    pc_ang_y=c_ang_y./t_y1;
    pc_ang_t=c_ang_t./t_t1;

    %% get the turning freq vs H,X,Y,T for reorientation event larger than
    %60 deg
    c_ang_hd1=bin_data_count(pre(idx),hdseries);
    c_ang_x1=bin_data_count(turn_x(idx),xseries);
    c_ang_y1=bin_data_count(turn_y(idx),yseries);
    c_ang_t1=bin_data_count(t0s(idx),tseries);
    if ii==1
        pc_ang_t_all=cell(ii_num,1);pc_ang_hd_all=cell(ii_num,1);pc_ang_x_all=cell(ii_num,1);pc_ang_y_all=cell(ii_num,1);
    end
    pc_ang_t_all{ii,1}=c_ang_t1./t_t2;
    pc_ang_hd_all{ii,1}=c_ang_hd1./t_hd2;
    pc_ang_x_all{ii,1}=c_ang_x1./t_x2;
    pc_ang_y_all{ii,1}=c_ang_y1./t_y2;

    c_ang_x3=bin_data_count2(pre(idx),hdseries,turn_x(idx),[0:xbin3:xmax]');
    c_ang_y3=bin_data_count2(pre(idx),hdseries,turn_y(idx),[0:ybin3:ymax]');
    c_ang_t3=bin_data_count2(pre(idx),hdseries,t0s(idx),[0:tbin3:tmax]');

    pc_ang_x3=c_ang_x3./t_x3;pc_ang_y3=c_ang_y3./t_y3;pc_ang_t3=c_ang_t3./t_t3;
    if ii==1
        pc_ang_x_dis=cell(ii_num,1);pc_ang_x_cent=cell(ii_num,1);pc_ang_x_prox=cell(ii_num,1);
        pc_ang_y_up=cell(ii_num,1);pc_ang_y_cent=cell(ii_num,1);pc_ang_y_low=cell(ii_num,1);
        pc_ang_t300=cell(ii_num,1);pc_ang_t600=cell(ii_num,1);pc_ang_t900=cell(ii_num,1);
    end

    pc_ang_x_dis{ii,1}=pc_ang_x3(:,1);pc_ang_x_cent{ii,1}=pc_ang_x3(:,2);pc_ang_x_prox{ii,1}=pc_ang_x3(:,3);
    pc_ang_y_up{ii,1}=pc_ang_y3(:,1);pc_ang_y_cent{ii,1}=pc_ang_y3(:,2);pc_ang_y_low{ii,1}=pc_ang_y3(:,3);
    pc_ang_t300{ii,1}=pc_ang_t3(:,1);pc_ang_t600{ii,1}=pc_ang_t3(:,2);pc_ang_t900{ii,1}=pc_ang_t3(:,3);


    % 9) plot the distribution of pre_turn deg and post turn deg
    c_pre=bin_data_count(pre,hdseries);
    p_pre=c_pre./length(pre);
    c_post=bin_data_count(post,hdseries);
    p_post=c_post./length(post);
    row1=3;col1=ii_num;
    c_reori=bin_data_count(ang,[0:10:180]);
    p_reori=c_reori./length(ang);
    color1=[color;color;color];
    mlt_subplt(p_pre,hd_series,110,row1,col1,ii,'','Pre turn (deg)','Frequency',color1,name(ii),'bar','ii_num',ii_num);
    mlt_subplt(p_post,hd_series,110,row1,col1,ii+ii_num,'','Post turn (deg)','Frequency',color1,'','bar','ii_num',ii_num);
    mlt_subplt(p_reori,[5:10:175]',110,row1,col1,ii+ii_num+ii_num,append('Distribution of pre_turn, post_turn degree, reorientatoin angle@bin=',num2str(bin),',tbin=',num2str(tbin),',hdbin=',num2str(hdbin), ',histograms'),'Reorientation angle (deg)','Frequency',color1,'','bar','ii_num',ii_num);

    clearvars -except ii* p_mean p_all row col name color p_t *max *bin *series p_*_all cond pc_x* pc_y* pc_*_all pc_t* deg_*_all outdir ac ver_turn_all* deg_x* deg_y* deg_t* ver_turn_ac* ac_x3* acl deg_l pc_ang_* ac_x3l ac_x3_seml
end


%% plot the turning f vs HD/x/y/t
mlt_subplt(pc_hd_all,hd_series,102,2,2,1,'','Heading Direction (deg)','Turn freq (s^-^1)',color,name,'line','ii_num',ii_num');
mlt_subplt(pc_x_all,x_series,102,2,2,2,'','X(mm)','Turn freq (s^-^1)',color,name,'line','ii_num',ii_num);
mlt_subplt(pc_y_all,y_series,102,2,2,3,'','Y(mm)','Turn freq (s^-^1)',color,name,'line','ii_num',ii_num);
mlt_subplt(pc_t_all,t_series,102,2,2,4,append('Turn freq (with turn events) vs HD,X,Y,T@bin=',num2str(bin),',tbin=',num2str(tbin),',hdbin=',num2str(hdbin)),'Time(s)','Turn freq (s^-^1)',color,name,'line','ii_num',ii_num);
%% plot the turn f vs HD in different region
mlt_subplt(pc_x_dis,hd_series,103,3,3,1,'','Heading Direction (deg)','Turn freq (s^-^1) vs x',color,name,'line','ii_num',ii_num,'title','Distant');
mlt_subplt(pc_x_cent,hd_series,103,3,3,2,'','Heading Direction (deg)','Turn freq (s^-^1) vs x',color,name,'line','ii_num',ii_num,'title','Central');
mlt_subplt(pc_x_prox,hd_series,103,3,3,3,'','Heading Direction (deg)','Turn freq (s^-^1) vs x',color,name,'line','ii_num',ii_num,'title','Proximate');
mlt_subplt(pc_y_low,hd_series,103,3,3,4,'','Heading Direction (deg)','Turn freq (s^-^1) vs y',color,name,'line','ii_num',ii_num,'title','Lower');
mlt_subplt(pc_y_cent,hd_series,103,3,3,5,'','Heading Direction (deg)','Turn freq (s^-^1) vs y',color,name,'line','ii_num',ii_num,'title','Central');
mlt_subplt(pc_y_up,hd_series,103,3,3,6,'','Heading Direction (deg)','Turn freq (s^-^1) vs y',color,name,'line','ii_num',ii_num,'title','Upper');
mlt_subplt(pc_t300,hd_series,103,3,3,7,'','Heading Direction (deg)','Turn freq (s^-^1) vs t',color,name,'line','ii_num',ii_num,'title','0-300s');
mlt_subplt(pc_t600,hd_series,103,3,3,8,'','Heading Direction (deg)','Turn freq (s^-^1) vs t',color,name,'line','ii_num',ii_num,'title','300-600s');
mlt_subplt(pc_t900,hd_series,103,3,3,9,append('Turn freq (with turn events) vs HD@bin=',num2str(bin),',tbin=',num2str(tbin),',hdbin=',num2str(hdbin)),'Heading Direction (deg)','Turn freq (s^-^1) vs t',color,name,'line','ii_num',ii_num,'title','600-900s');
% plot the turning f vs HD/x/y/t (>60)
mlt_subplt(pc_ang_hd_all,hd_series,111,2,2,1,'','Heading Direction (deg)',append('Turn freq (>',num2str(deg_l),',s^-^1)'),color,name,'line','ii_num',ii_num');
mlt_subplt(pc_ang_x_all,x_series,111,2,2,2,'','X(mm)',append('Turn freq (>',num2str(deg_l),',s^-^1)'),color,name,'line','ii_num',ii_num');
mlt_subplt(pc_ang_y_all,y_series,111,2,2,3,'','Y(mm)',append('Turn freq (>',num2str(deg_l),',s^-^1)'),color,name,'line','ii_num',ii_num');
mlt_subplt(pc_ang_t_all,t_series,111,2,2,4,append('Large turning freq vs HD,X,Y,T@bin=',num2str(bin),',tbin=',num2str(tbin),',hdbin=',num2str(hdbin)),'Time(s)',append('Turn freq (>',num2str(deg_l),',s^-^1)'),color,name,'line','ii_num',ii_num');
% plot the turning freq vs HD for three regions of X,Y,T
mlt_subplt(pc_ang_x_dis,hd_series,112,3,3,1,'','Heading Direction (deg)',append('Turn freq (>',num2str(deg_l),',s^-^1) vs x'),color,name,'line','ii_num',ii_num,'title','Distant');
mlt_subplt(pc_ang_x_cent,hd_series,112,3,3,2,'','Heading Direction (deg)',append('Turn freq (>',num2str(deg_l),',s^-^1) vs x'),color,name,'line','ii_num',ii_num,'title','Central');
mlt_subplt(pc_ang_x_prox,hd_series,112,3,3,3,'','Heading Direction (deg)',append('Turn freq (>',num2str(deg_l),',s^-^1) vs x'),color,name,'line','ii_num',ii_num,'title','Proximate');
mlt_subplt(pc_ang_y_low,hd_series,112,3,3,4,'','Heading Direction (deg)',append('Turn freq (>',num2str(deg_l),',s^-^1) vs y'),color,name,'line','ii_num',ii_num,'title','Lower');
mlt_subplt(pc_ang_y_cent,hd_series,112,3,3,5,'','Heading Direction (deg)',append('Turn freq (>',num2str(deg_l),',s^-^1) vs y'),color,name,'line','ii_num',ii_num,'title','Central');
mlt_subplt(pc_ang_y_up,hd_series,112,3,3,6,'','Heading Direction (deg)',append('Turn freq (>',num2str(deg_l),',s^-^1) vs y'),color,name,'line','ii_num',ii_num,'title','Upper');
mlt_subplt(pc_ang_t300,hd_series,112,3,3,7,'','Heading Direction (deg)',append('Turn freq (>',num2str(deg_l),',s^-^1) vs t'),color,name,'line','ii_num',ii_num,'title','0-300s');
mlt_subplt(pc_ang_t600,hd_series,112,3,3,8,'','Heading Direction (deg)',append('Turn freq (>',num2str(deg_l),',s^-^1) vs t'),color,name,'line','ii_num',ii_num,'title','300-600s');
mlt_subplt(pc_ang_t900,hd_series,112,3,3,9,append('Large turning freq vs HD@bin=',num2str(bin),',tbin=',num2str(tbin),',hdbin=',num2str(hdbin)),'Heading Direction (deg)',append('Turn freq (>',num2str(deg_l),',s^-^1) vs t'),color,name,'line','ii_num',ii_num,'title','600-900s');


%% plot reorient deg vs HD/X/Y/T
mlt_subplt(deg_hd_all,hd_series,104,2,2,1,'','Heading Direction (deg)','Reorientation angle (deg)',color,name,'line','ii_num',ii_num,'ebar',1);
mlt_subplt(deg_x_all,x_series,104,2,2,2,'','X(mm)','Reorientation angle (deg)',color,name,'line','ii_num',ii_num,'ebar',1);
mlt_subplt(deg_y_all,y_series,104,2,2,3,'','Y(mm)','Reorientation angle (deg)',color,name,'line','ii_num',ii_num,'ebar',1);
mlt_subplt(deg_t_all,t_series,104,2,2,4,append('Reori angle vs HD,X,Y,T@bin=',num2str(bin),',tbin=',num2str(tbin),',hdbin=',num2str(hdbin)),'Time(s)','Reorientation angle (deg)',color,name,'line','ii_num',ii_num,'ebar',1);
%% plot reorient deg vs HD in different region of x,y,t
mlt_subplt(deg_x_dis,hd_series,105,3,3,1,'','Heading Direction (deg)','Reori angle (deg) vs x',color,name,'line','ii_num',ii_num,'title','Distant','ebar',1);
mlt_subplt(deg_x_cent,hd_series,105,3,3,2,'','Heading Direction (deg)','Reori angle (deg) vs x',color,name,'line','ii_num',ii_num,'title','Central','ebar',1);
mlt_subplt(deg_x_prox,hd_series,105,3,3,3,'','Heading Direction (deg)','Reori angle (deg) vs x',color,name,'line','ii_num',ii_num,'title','Proximate','ebar',1);

mlt_subplt(deg_y_low,hd_series,105,3,3,4,'','Heading Direction (deg)','Reori angle (deg) vs y',color,name,'line','ii_num',ii_num,'title','Lower','ebar',1);
mlt_subplt(deg_y_cent,hd_series,105,3,3,5,'','Heading Direction (deg)','Reori angle (deg) vs y',color,name,'line','ii_num',ii_num,'title','Central','ebar',1);
mlt_subplt(deg_y_up,hd_series,105,3,3,6,'','Heading Direction (deg)','Reori angle (deg) vs y',color,name,'line','ii_num',ii_num,'title','Upper','ebar',1);

mlt_subplt(deg_t300,hd_series,105,3,3,7,'','Heading Direction (deg)','Reori angle (deg) vs t',color,name,'line','ii_num',ii_num,'title','0-300s','ebar',1);
mlt_subplt(deg_t600,hd_series,105,3,3,8,'','Heading Direction (deg)','Reori angle (deg) vs t',color,name,'line','ii_num',ii_num,'title','300-600s','ebar',1);
mlt_subplt(deg_t900,hd_series,105,3,3,9,append('Reori angle vs HD@bin=',num2str(bin),',tbin=',num2str(tbin),',hdbin=',num2str(hdbin)),'Heading Direction (deg)','Reori angle (deg) vs t',color,name,'line','ii_num',ii_num,'title','600-0-900s','ebar',1);

%% acceptance rate
mlt_subplt(ac,[],106,2,3,1,'',"Groups","Acceptance rate",color,name,"bar",'ii_num',ii_num,'ebar',1);
ac_x3_all{1,1}=ac_x3;ac_x3_all{1,2}=ac_x3_sem;
ac_x3_all_l{1,1}=ac_x3l;ac_x3_all_l{1,2}=ac_x3_seml;
mlt_subplt(ac_x3_all,[],106,2,3,2,'',"Groups","Acceptance rate",color,name,"bar",'ii_num',ii_num,'ebar',1,'legends',["Distant";"Central";"Proximate"]);
mlt_subplt(acl,[],106,2,3,4,'',"Groups",append("Acceptance rate (>",num2str(deg_l),")"),color,name,"bar",'ii_num',ii_num,'ebar',1,'title',append('Acceptance rate for reorientation >', num2str(deg_l)));
mlt_subplt(ac_x3_all_l,[],106,2,3,5,'Acceptance rate',"Groups",append("Acceptance rate (>",num2str(deg_l),")"),color,name,"bar",'ii_num',ii_num,'ebar',1,'title',append('Acceptance rate for reorientation >', num2str(deg_l)),'legends',["Distant";"Central";"Proximate"]);

group=['Conditions';name];

str1=["Number of total turn";string(ver_turn_all)];
str2=["Number of turn toward odor";string(ver_turn_ac)];
str3=[append("Number of total turn >",num2str(deg_l));string(ver_turn_all_l)];
str4=[append("Number of turn toward odor >",num2str(deg_l));string(ver_turn_ac_l)];
figure(106)
subplot(2,3,2)
hold on
text(ii_num+2,max(ac(:,1)),group)
text(ii_num+3,max(ac(:,1)),str1);
text(ii_num+4,max(ac(:,1)),str2);
text(ii_num+5,max(ac(:,1)),str3);
text(ii_num+6,max(ac(:,1)),str4);
hold off
ac
save_all_figures(outdir);
close all
clear
end