function []=opto_history_properties(tracker_num,genotype,condition,filename,output_name,varargin)
%% get the directory for all the data
w=1;

for i=1:length(genotype)
    for j=1:length(condition)
        cond(w,1)=fullfile("/project/6010970/screen/olfactory_output/JB_JAABA",tracker_num,genotype{i,1},condition{j,1},filename);
        name(w,1)=append(genotype{i,1},"@",condition{j,1});
        w=w+1;
    end
end
%make sure there is a file exists
a=1;
del=[];
for i=1:length(cond)
    if ~isfile(cond(i))
        del(a,1)=i;
        a=a+1;
    end
end

if ~isempty(del)
    cond(del)=[];
end

cond

clear a del
%% set the output name
outdir=fullfile("/project/6010970/screen/olfactory_output/properties_opto_turn",output_name);
ii_num=length(cond);
%set the properties for plots
row=3;col=3;
color=[0,0,0;0.627,0.82,1;0.659,0.659,0.659;0.871,0.416,0.451];
t_pre=10;%how many second we want to go back to calculate the change in odor gradient
t_post=10; %same as t_pre
deg_list=60;% this is used for the cutoff for the turn size -->can be an arry (size ii_num*1) or a number
input_cond=["ctr1";"exp1";"exp1";"ctr2";"exp2";"exp2"];
thred=0.2;%turning event with a starting point below this thredshold will be eliminated
row_comp=1; col_comp=2;
tbin=0.5; % binning for the time, make sure it is not smaller than the timeframe of the tracking (around 0.2s)
t_el=600; % if turning events start after 600s, it will be eliminated
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
    elseif strcmp ('t_bef',varargin{i})
        t_pre=varargin{i+1};
    elseif strcmp ('t_aft',varargin{i})
        t_post=varargin{i+1};
    elseif strcmp ('thred',varargin{i})
        thred=varargin{i+1};
    elseif strcmp(varargin{i},'deg_l')
        deg_list=varargin{i+1};
    elseif strcmp(varargin{i},'input_cond')
        input_cond=string(varargin{i+1});
    elseif strcmp(varargin{i},'row_comp')
        row_comp=varargin{i+1};
    elseif strcmp(varargin{i},'col_comp')
        col_comp=varargin{i+1};
    elseif strcmp(varargin{i},'t_el')
        t_el=varargin{i+1};
    end
end
%% load data
for ii=1:ii_num
    load(cond(ii));
    %% 1) plot the x vs opto intensity just to confirm the shape of the fitted curve
    mlt_subplt(vertcat(dat.grad_x{:}),vertcat(dat.x{:}),1,row,col,ii,'x vs opto gradient', ...
        'X position(mm)','Light Intensity (uW/mm^2)',color,name,'scatter','marker_size',2,'alpha',1);
    %the results will be a single figure holding curve for all conditions

    %% 2) plot the distribution of opto intensity for turning events
    % a) use t0&cen of mass b) use midpoint of turning &cen of mass
    % c) use t0&head of turning

    series=0:0.05:3; % set bins for light intensity

    [c_a]=bin_data_count(vertcat(dat.grad_t0{:}),series);
    [c_b]=bin_data_count(vertcat(dat.grad_turn_x{:}),series);
    [c_c]=bin_data_count(vertcat(dat.grad_turn_t0_head{:}),series);

    f_a=c_a./sum(c_a);
    f_b=c_b./sum(c_b);
    f_c=c_c./sum(c_c);
    %plot the count vs the opto series for each different methods
    for i=1:length(series)-1
        series1(i,1)=(series(i)+series(i+1))/2;
    end

    mlt_subplt([f_a,f_b,f_c],series1,2,row,col,ii,'Light intensity of turns','Light Intensity of turns (uW/mm^2)','Frequency',[0 0 1; 1 0 1;0 0 0],name(ii),'multiple lines','ii_num',ii_num,'legends',["t0&cen of mass";"mid&cen of mass";"t0&head"])

    clear f_* c_* series*
    %% 3 i) get the data
    grad_t0=dat.grad_t0;
    grad_t1=dat.grad_t1;
    t0s=dat.t0s;
    t1s=dat.t1s;
    reori_deg=dat.reorient_deg_abs;
    grad_x=dat.grad_x;
    t=dat.et;

    %% 3 ii) eliminate the turning with a opto gradient lower than the given threshold  eliminate the turning events that occur after t_el
    len=length(grad_t0);
    grad_t0_el=cell(len,1);
    grad_t1_el=cell(len,1);
    t0s_el=cell(len,1);
    t1s_el=cell(len,1);
    grad_t0_el=cell(len,1);
    deg_el=cell(len,1);
    for j=1:len
        if ~isempty(grad_t0{j,1})
            idx_a=find(grad_t0{j,1}>thred&t0s{j,1}<t_el);

            grad_t0_el{j,1}=grad_t0{j,1}(idx_a);
            grad_t1_el{j,1}=grad_t1{j,1}(idx_a);
            t0s_el{j,1}=t0s{j,1}(idx_a);
            t1s_el{j,1}=t1s{j,1}(idx_a);
            deg_el{j,1}=reori_deg{j,1}(idx_a);

            clear idx*

        end
    end


    %% 4) plot the the opto intensity <<t_bef>> seconds before and <<t_aft>> seconds after the turning
    % use t0&cen of mass
    tseries_pre=[-t_pre:tbin:0]';
    tseries_post=[0:tbin:t_post]';
    % A) data with elimination B) data without elimination
    a=1;;% a and b is used for the idx to save data of pre and post-turn data
    % for data with elimination
    for j=1:len
        t0s_larva_el=t0s_el{j,1};
        t1s_larva_el=t1s_el{j,1};
        t_larva=t{j,1};
        int_x_larva=grad_x{j,1};

        if ~isempty(t0s_larva_el) %a) use t0&center of mass
            a_pre=t0s_larva_el-t_pre;
            a_post=t1s_larva_el+t_post;
            %if the turning events do not have enough data to go back to
            %its history, set the start point the start of its tracking
            %time
            a_pre(a_pre<=t_larva(1))=t_larva(1);
            a_post(a_post>t_larva(end))=t_larva(end);
            for i=1:length(a_pre)
                % get the index for the <<t_bef>> and the start of
                % turns
                [~,I(i,1)]=min(abs(t_larva-a_pre(i)));
                [~,I1(i,1)]=min(abs(t_larva-t0s_larva_el(i)));
                list=I(i,1):I1(i,1);
                clear I I1
                % get the light intensity and time
                int_pre=int_x_larva(list);
                t_pre1=t_larva(list)-t0s_larva_el(i);
                % compute the light gradient
                deri=diff(int_pre)./diff(t_pre1);
                % compute the light gradient according to the time bin
                for k=1:length(tseries_pre)
                    [~,I(k,1)]=min(abs(t_pre1-tseries_pre(k)));
                end
                % two ways to compute: get two points or use the average
                % value in that each timebin
                deri_pre_arry{a,1}=(diff(int_pre(I))./diff(t_pre1(I)))';
                t_pre1(1)=[]; % make sure the deri and t_pre1 have the same length
                for k=1:length(tseries_pre)-1
                    idx=find(t_pre1>tseries_pre(k)&t_pre1<=tseries_pre(k+1));
                    deri_pre_arry2{a,1}(1,k)=mean(deri(idx));
                    clear idx
                end

                clear list t_pre1 int_pre I deri

                [~,I2(i,1)]=min(abs(t_larva-a_post(i)));
                [~,I3(i,1)]=min(abs(t_larva-t1s_larva_el(i)));
                list=I3(i,1):I2(i,1);
                int_post=int_x_larva(list);
                t_post1=t_larva(list)-t1s_larva_el(i);

                % compute the light gradient
                deri=diff(int_post)./diff(t_post1);
                % similar to the code in 166-177
                for k=1:length(tseries_pre)
                    [~,I(k,1)]=min(abs(t_post1-tseries_post(k)));
                end
                deri_post_arry{a,1}=(diff(int_post(I))./diff(t_post1(I)))';

                t_post1(end)=[];
                for k=1:length(tseries_post)-1
                    idx=find(t_post1>tseries_post(k)&t_post1<=tseries_post(k+1));
                    deri_post_arry2{a,1}(1,k)=mean(deri(idx));
                    clear idx
                end
                a=a+1;
                clear I* idx t_post1 int_post deri
            end
        end
    end

    % concatenate all the cell arrys
    deri_pre=cell2mat(deri_pre_arry);
    deri_post=cell2mat(deri_post_arry);


    %compute the mean light intensity using the data
    if ii==1
        deri_pre_all=cell(ii_num,1); deri_post_all=cell(ii_num,1);
    end
    deri_pre_all{ii,1}(:,1)=mean(deri_pre,'omitnan');
    deri_post_all{ii,1}(:,1)=mean(deri_post,'omitnan');
    
   
    % Compute the sem
%     deri_pre_all{ii,1}(:,2)=std(deri_pre,'omitnan')./sqrt(height(deri_pre));
%     deri_post_all{ii,1}(:,2)=std(deri_post,'omitnan')./sqrt(height(deri_post));

    % concatenate all the cell arrys
    deri_pre2=cell2mat(deri_pre_arry2);
    deri_post2=cell2mat(deri_post_arry2);


    %compute the mean light intensity using the data
    if ii==1
        deri_pre_all2=cell(ii_num,1); deri_post_all2=cell(ii_num,1);
    end
    deri_pre_all2{ii,1}(:,1)=mean(deri_pre2,'omitnan');
    deri_post_all2{ii,1}(:,1)=mean(deri_post2,'omitnan');

%     deri_pre_all2{ii,1}(:,2)=std(deri_pre2,'omitnan')./sqrt(height(deri_pre2));
%     deri_post_all2{ii,1}(:,2)=std(deri_post2,'omitnan')./sqrt(height(deri_post2));

    
     % Compute the 95% CI
    for k=1:width(deri_pre)
        pd_pre=fitdist(deri_pre(:,k),'Normal');
        pd_post=fitdist(deri_post(:,k),'Normal');
        
        ci=paramci(pd_pre);
        ci2=paramci(pd_post);
        
        ci_pre(k,1)=ci(2,1)-pd_pre.mu; 
        ci_post(k,1)=ci2(2,1)-pd_post.mu; 
        
        clear ci2 ci pd_pre pd_post
        pd_pre=fitdist(deri_pre2(:,k),'Normal');
        pd_post=fitdist(deri_post2(:,k),'Normal');
        
        ci=paramci(pd_pre);
        ci2=paramci(pd_post);
        
        ci_pre2(k,1)=ci(2,1)-pd_pre.mu; 
        ci_post2(k,1)=ci2(2,1)-pd_post.mu; 
        clear ci ci2
    end 
    
     deri_pre_all{ii,1}(:,2)=ci_pre;
     deri_post_all{ii,1}(:,2)=ci_post;
     
     deri_pre_all2{ii,1}(:,2)=ci_pre2;
     deri_post_all2{ii,1}(:,2)=ci_post2;
     
    %% 5) plot the light intensity vs t for large and small turning events


    if length(deg_list)==1
        deg_l=deg_list;
    else
        deg_l=deg_list(ii);
    end
    disp(append("The cutoff for large turning event of the group ", name(ii)," is ",num2str(deg_l)," deg"));

    % Split the data into large and small turning events
    deg_turn_el=vertcat(deg_el{:});
    idx_L=find(deg_turn_el>=deg_l);
    idx_s=find(deg_turn_el<deg_l);

    deri_pre_L=deri_pre(idx_L,:);
    deri_post_L=deri_post(idx_L,:);

    deri_pre_s=deri_pre(idx_s,:);
    deri_post_s=deri_post(idx_s,:);

    % Compute the mean light intensity of large and small turning events
    if ii==1
        deri_pre_L_all=cell(ii_num,1); deri_post_L_all=cell(ii_num,1);
        deri_pre_s_all=cell(ii_num,1); deri_post_s_all=cell(ii_num,1);
    end

    deri_pre_L_all{ii,1}(:,1)=mean(deri_pre_L,'omitnan');
    deri_post_L_all{ii,1}(:,1)=mean(deri_post_L,'omitnan');

%     deri_pre_L_all{ii,1}(:,2)=std(deri_pre_L,'omitnan')./sqrt(height(deri_pre_L));
%     deri_post_L_all{ii,1}(:,2)=std(deri_post_L,'omitnan')./sqrt(height(deri_post_L));

    deri_pre_s_all{ii,1}(:,1)=mean(deri_pre_s,'omitnan');
    deri_post_s_all{ii,1}(:,1)=mean(deri_post_s,'omitnan');

%     deri_pre_s_all{ii,1}(:,2)=std(deri_pre_s,'omitnan')./sqrt(height(deri_pre_s));
%     deri_post_s_all{ii,1}(:,2)=std(deri_post_s,'omitnan')./sqrt(height(deri_post_s));


    % Split the data into large and small turning events

    deri_pre_L2=deri_pre2(idx_L,:);
    deri_post_L2=deri_post2(idx_L,:);

    deri_pre_s2=deri_pre2(idx_s,:);
    deri_post_s2=deri_post2(idx_s,:);

    % Compute the mean light intensity of large and small turning events
    if ii==1
        deri_pre_L_all2=cell(ii_num,1); deri_post_L_all2=cell(ii_num,1);
        deri_pre_s_all2=cell(ii_num,1); deri_post_s_all2=cell(ii_num,1);
    end

    deri_pre_L_all2{ii,1}(:,1)=mean(deri_pre_L2,'omitnan');
    deri_post_L_all2{ii,1}(:,1)=mean(deri_post_L2,'omitnan');

%     deri_pre_L_all2{ii,1}(:,2)=std(deri_pre_L2,'omitnan')./sqrt(height(deri_pre_L2));
%     deri_post_L_all2{ii,1}(:,2)=std(deri_post_L2,'omitnan')./sqrt(height(deri_post_L2));

    deri_pre_s_all2{ii,1}(:,1)=mean(deri_pre_s2,'omitnan');
    deri_post_s_all2{ii,1}(:,1)=mean(deri_post_s2,'omitnan');

%     deri_pre_s_all2{ii,1}(:,2)=std(deri_pre_s2,'omitnan')./sqrt(height(deri_pre_s2));
%     deri_post_s_all2{ii,1}(:,2)=std(deri_post_s2,'omitnan')./sqrt(height(deri_post_s2));

             % Compute the 95% CI
    for k=1:width(deri_pre_s)
        pd_pre_L=fitdist(deri_pre_L(:,k),'Normal');
        pd_post_L=fitdist(deri_post_L(:,k),'Normal');
        
        ci_L=paramci(pd_pre_L);
        ci_L2=paramci(pd_post_L);
        
        ci_pre_L(k,1)=ci_L(2,1)-pd_pre_L.mu; 
        ci_post_L(k,1)=ci_L2(2,1)-pd_post_L.mu; 
        clear ci_L2 ci_L pd_pre_L pd_post_L
        
        pd_pre_s=fitdist(deri_pre_s(:,k),'Normal');
        pd_post_s=fitdist(deri_post_s(:,k),'Normal');
        
        ci_s=paramci(pd_pre_s);
        ci_s2=paramci(pd_post_s);
        
        ci_pre_s(k,1)=ci_s(2,1)-pd_pre_s.mu; 
        ci_post_s(k,1)=ci_s2(2,1)-pd_post_s.mu; 
        clear ci_s2 ci_s pd_pre_s pd_post_s
         
        pd_pre_L=fitdist(deri_pre_L2(:,k),'Normal');
        pd_post_L=fitdist(deri_post_L2(:,k),'Normal');
        
        ci_L=paramci(pd_pre_L);
        ci_L2=paramci(pd_post_L);
        
        ci_pre_L2(k,1)=ci_L(2,1)-pd_pre_L.mu; 
        ci_post_L2(k,1)=ci_L2(2,1)-pd_post_L.mu; 
        clear ci_L2 ci_L pd_pre_L pd_post_L
        
        pd_pre_s=fitdist(deri_pre_s2(:,k),'Normal');
        pd_post_s=fitdist(deri_post_s2(:,k),'Normal');
        
        ci_s=paramci(pd_pre_s);
        ci_s2=paramci(pd_post_s);
        
        ci_pre_s2(k,1)=ci_s(2,1)-pd_pre_s.mu; 
        ci_post_s2(k,1)=ci_s2(2,1)-pd_post_s.mu; 
        clear ci_s2 ci_s pd_pre_s pd_post_s
         
    end 
    %save the 95% CI
    deri_pre_s_all2{ii,1}(:,2)=ci_pre_s2;
    deri_pre_L_all2{ii,1}(:,2)=ci_pre_L2;
    deri_pre_s_all{ii,1}(:,2)=ci_pre_s;
    deri_pre_L_all{ii,1}(:,2)=ci_pre_L;
    
    
    deri_post_s_all2{ii,1}(:,2)=ci_post_s2;
    deri_post_L_all2{ii,1}(:,2)=ci_post_L2;
    deri_post_s_all{ii,1}(:,2)=ci_post_s;
    deri_post_L_all{ii,1}(:,2)=ci_post_L;
    clearvars -except dat len row* col* ii* name color*  cond thred outdir deg_list deri_*_all t_pre t_post input_cond tbin tseries* t_el deri_*_all2


end
%% plot the light intensity <<t_bef>> and <<t_aft>> relative to t0s and t1s (turns occur in the zone with light intensity less than <<thred>> are deleted)

% Plot the derivative of light intensity for each ctr and exp group
tseries=[linspace(-t_pre,0,length(tseries_pre)-1)]';
tseries1=[linspace(0,t_post,length(tseries_post)-1)]';
if ~isempty(input_cond)
    idx=find(contains(input_cond,"ctr"));
    [ctr_deri_pre,exp_deri_pre,~,~,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(deri_pre_all,idx,name,color,'input_cond',input_cond);
    [ctr_deri_pre_s,exp_deri_pre_s]=split_data_ctr_exp(deri_pre_s_all,idx,name,color,'input_cond',input_cond);
    [ctr_deri_pre_L,exp_deri_pre_L]=split_data_ctr_exp(deri_pre_L_all,idx,name,color,'input_cond',input_cond);

    [ctr_deri_post,exp_deri_post]=split_data_ctr_exp(deri_post_all,idx,name,color,'input_cond',input_cond);
    [ctr_deri_post_s,exp_deri_post_s]=split_data_ctr_exp(deri_post_s_all,idx,name,color,'input_cond',input_cond);
    [ctr_deri_post_L,exp_deri_post_L]=split_data_ctr_exp(deri_post_L_all,idx,name,color,'input_cond',input_cond);

    [ctr_deri_pre2,exp_deri_pre2]=split_data_ctr_exp(deri_pre_all2,idx,name,color,'input_cond',input_cond);
    [ctr_deri_pre_s2,exp_deri_pre_s2]=split_data_ctr_exp(deri_pre_s_all2,idx,name,color,'input_cond',input_cond);
    [ctr_deri_pre_L2,exp_deri_pre_L2]=split_data_ctr_exp(deri_pre_L_all2,idx,name,color,'input_cond',input_cond);

    [ctr_deri_post2,exp_deri_post2]=split_data_ctr_exp(deri_post_all2,idx,name,color,'input_cond',input_cond);
    [ctr_deri_post_s2,exp_deri_post_s2]=split_data_ctr_exp(deri_post_s_all2,idx,name,color,'input_cond',input_cond);
    [ctr_deri_post_L2,exp_deri_post_L2]=split_data_ctr_exp(deri_post_L_all2,idx,name,color,'input_cond',input_cond);



    [~,ind,ind0]=unique(ctr_name);
    for k=1:length(ind)
        legends=[ctr_name(ind(k),1);exp_name((ind0==k),1)];
        color1=[ctr_color(ind(k),:);exp_color((ind0==k),:)];

        data_pre={ctr_deri_pre{ind(k),1};exp_deri_pre(ind0==k)};
        data_pre=vertcat(data_pre{:});

        data_pre_s={ctr_deri_pre_s{ind(k),1};exp_deri_pre_s(ind0==k)};
        data_pre_s=vertcat(data_pre_s{:});

        data_pre_L={ctr_deri_pre_L{ind(k),1};exp_deri_pre_L(ind0==k)};
        data_pre_L=vertcat(data_pre_L{:});

        data_post={ctr_deri_post{ind(k),1};exp_deri_post(ind0==k)};
        data_post=vertcat(data_post{:});

        data_post_s={ctr_deri_post_s{ind(k),1};exp_deri_post_s(ind0==k)};
        data_post_s=vertcat(data_post_s{:});

        data_post_L={ctr_deri_post_L{ind(k),1};exp_deri_post_L(ind0==k)};
        data_post_L=vertcat(data_post_L{:});

        data_pre2={ctr_deri_pre2{ind(k),1};exp_deri_pre2(ind0==k)};
        data_pre2=vertcat(data_pre2{:});

        data_pre_s2={ctr_deri_pre_s2{ind(k),1};exp_deri_pre_s2(ind0==k)};
        data_pre_s2=vertcat(data_pre_s2{:});

        data_pre_L2={ctr_deri_pre_L2{ind(k),1};exp_deri_pre_L2(ind0==k)};
        data_pre_L2=vertcat(data_pre_L2{:});

        data_post2={ctr_deri_post2{ind(k),1};exp_deri_post2(ind0==k)};
        data_post2=vertcat(data_post2{:});

        data_post_s2={ctr_deri_post_s2{ind(k),1};exp_deri_post_s2(ind0==k)};
        data_post_s2=vertcat(data_post_s2{:});

        data_post_L2={ctr_deri_post_L2{ind(k),1};exp_deri_post_L2(ind0==k)};
        data_post_L2=vertcat(data_post_L2{:});

        mlt_subplt(data_pre,tseries,100,row_comp,col_comp,k,'','','',color1,legends,'line','ii_num',length(data_pre),'ebar',1);
        mlt_subplt(data_post,tseries1,100,row_comp,col_comp,k,append('Light intensity derivative',num2str(t_pre),' s before and after turn and threshold is ', num2str(thred),'_eliminate last ', num2str(900-t_el),'s'),'Time (s)','Derivative of light intensity (uW/mm^2*s)',color1,legends,'line','ii_num',length(data_post),'ebar',1);
        mlt_subplt(data_pre_L,tseries,101,row_comp,col_comp,k,'','','',color1,legends,'line','ii_num',length(data_pre_L),'ebar',1);
        mlt_subplt(data_post_L,tseries1,101,row_comp,col_comp,k,append('Large turn_Light intensity derivative',num2str(t_pre),' s before and after turn and threshold is ', num2str(thred),'_eliminate last ', num2str(900-t_el),'s'),'Time (s)','Derivative of light intensity (uW/mm^2*s,>rms turn size)',color1,legends,'line','ii_num',length(data_post_L),'ebar',1);
        mlt_subplt(data_pre_s,tseries,102,row_comp,col_comp,k,'','','',color1,legends,'line','ii_num',length(data_pre_s),'ebar',1);
        mlt_subplt(data_post_s,tseries1,102,row_comp,col_comp,k,append('Small turn_Light intensity derivative',num2str(t_pre),' s before and after turn and threshold is ', num2str(thred),'_eliminate last ', num2str(900-t_el),'s'),'Time (s)','Derivative of light intensity (uW/mm^2*s,<rms turn size)',color1,legends,'line','ii_num',length(data_post_s),'ebar',1);
        % use the mean value to compute the deri
        mlt_subplt(data_pre2,tseries,103,row_comp,col_comp,k,'','','',color1,legends,'line','ii_num',length(data_pre),'ebar',1);
        mlt_subplt(data_post2,tseries1,103,row_comp,col_comp,k,append('Using mean value_Light intensity derivative',num2str(t_pre),' s before and after turn and threshold is ', num2str(thred),'_eliminate last ', num2str(900-t_el),'s'),'Time (s)','Derivative of light intensity (uW/mm^2*s)',color1,legends,'line','ii_num',length(data_post),'ebar',1);
        mlt_subplt(data_pre_L2,tseries,104,row_comp,col_comp,k,'','','',color1,legends,'line','ii_num',length(data_pre_L),'ebar',1);
        mlt_subplt(data_post_L2,tseries1,104,row_comp,col_comp,k,append('Using mean value_Large turn_Light intensity derivative',num2str(t_pre),' s before and after turn and threshold is ', num2str(thred),'_eliminate last ', num2str(900-t_el),'s'),'Time (s)','Derivative of light intensity (uW/mm^2*s,>rms turn size)',color1,legends,'line','ii_num',length(data_post_L),'ebar',1);
        mlt_subplt(data_pre_s2,tseries,105,row_comp,col_comp,k,'','','',color1,legends,'line','ii_num',length(data_pre_s),'ebar',1);
        mlt_subplt(data_post_s2,tseries1,105,row_comp,col_comp,k,append('Using mean value_Small turn_Light intensity derivative',num2str(t_pre),' s before and after turn and threshold is ', num2str(thred),'_eliminate last ', num2str(900-t_el),'s'),'Time (s)','Derivative of light intensity (uW/mm^2*s,<rms turn size)',color1,legends,'line','ii_num',length(data_post_s),'ebar',1);
    end

end



%% save figures
outdir1=fullfile(outdir,append(num2str(t_pre),"Light gradient before and after turn"));
if ~isfolder(outdir1)
    mkdir(outdir1);
end
save_all_figures(outdir1);
close all;
end