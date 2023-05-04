function group_JB_JAABA_data_small(tracker_num,genotype,condition,odor_pos)
% this code is for exp using the 9cm agar plate
% eliminate the edge data using the 'r', radius
% odor_pos is the location of the odor
%set output folder
outdir=fullfile("/project/6010970/screen/olfactory_output/JB_JAABA",tracker_num,genotype,condition);
%% 1) load JB data

dir_JB=fullfile("/project/6010970/screen/jb-results",tracker_num,genotype,condition);
[dat_JB,uname_JB]=get_JB_data(dir_JB);
uname_JB
height(dat_JB)
%get x and y from the xspine and yspine
for i=1:length(dat_JB.xspine)
    dat_JB.x{i,1}=dat_JB.xspine{i,1}(:,6);
    dat_JB.y{i,1}=dat_JB.yspine{i,1}(:,6);
end
%% 2) load JAABA data

dir_JAABA=fullfile("/project/6010970/screen/JAABA_processed",tracker_num,genotype,condition);
[dat_JAABA,uname_JAABA]=get_JAABA_data(dir_JAABA);
uname_JAABA
height(dat_JAABA)
%% 3) rearrange the data array for the jaaba, calculate the total
% tracking time and body length for each object
% dat=align_JB_JAABA_data(dat_JB,dat_JAABA);
idx_JB=find(diff(dat_JB.AN)<0);
idx_JB=[0;idx_JB;height(dat_JB)];
idx_JAABA=find(diff(dat_JAABA.AN)<0);
idx_JAABA=[0;idx_JAABA;height(dat_JAABA)];

if length(idx_JB)~=length(idx_JAABA)
    warning('The number of timestamps in the jb_results and JAABA_processed are different!');
    quit
end
%find the animal number that is in JAABA but not in JB and delete it
for j=length(idx_JB):-1:2
    arry_JB=[idx_JB(j-1)+1:idx_JB(j)]';
    arry_JAABA=[idx_JAABA(j-1)+1:idx_JAABA(j)]';
    AN_JB=dat_JB.AN(arry_JB,1);
    AN_JAABA=dat_JAABA.AN(arry_JAABA,1);
    ind=find(~ismember(AN_JB,AN_JAABA));
    dat_JB(arry_JB(ind),:)=[];

    ind1=find(~ismember(AN_JAABA,AN_JB));
    dat_JAABA(arry_JAABA(ind1),:)=[];

    clear AN_* arry* ind*
end

if length(dat_JB.AN)==length(dat_JAABA.AN)
    deleted_JB_length=length(dat_JB.AN)
else
    warning('The length of JB and JAABA array are different!');
    quit
end

dat_JAABA=removevars(dat_JAABA,"AN");
dat=[dat_JB dat_JAABA];

% rotate the x,y to make the odor pos on the right side-->
% counter-clockwise 90 deg


for i=1:length(dat.x)
    x=-dat.y{i,1};
    y=dat.x{i,1};
    xspine=-dat.yspine{i,1};
    yspine=dat.xspine{i,1};
    dat.y{i,1}=y;
    dat.yspine{i,1}=yspine;
    dat.x{i,1}=x;
    dat.xspine{i,1}=xspine;
    clear x y xspine yspine
end

ymin=min(vertcat(dat.y{:}));
xmin=min(vertcat(dat.x{:}));

for i=1:length(dat.x)
    dat.y{i,1}=dat.y{i,1}-ymin;
    dat.yspine{i,1}=dat.yspine{i,1}-ymin;
    dat.x{i,1}=dat.x{i,1}-xmin;
    dat.xspine{i,1}=dat.xspine{i,1}-xmin;
end
%% calculate tsum,midline,and change t0s from horizontal arry to vertical

for i=length(dat.x):-1:1
    if length(dat.x{i,1})<3
        dat(i,:)=[];
    else

        dat.tsum(i,1)=dat.et{i,1}(end)-dat.et{i,1}(1);
        dat.mid(i,1)=sqrt((dat.xspine{i,1}(1,11)-dat.xspine{i,1}(1,1))^2+(dat.yspine{i,1}(1,11)-dat.yspine{i,1}(1,1))^2);
        dat.t0s{i,1}=(dat.t0s{i,1})';
        dat.t1s{i,1}=(dat.t1s{i,1})';
    end
end

%% 4) find the index for t0s and t1s

for i=1:length(dat.t0s)
    dat0=dat.t0s{i,1};
    dat1=dat.t1s{i,1};
    time=dat.et{i,1};

    %find the t0s and t1s that is out of the time range and delete them
    if ~isempty(dat0)&&~isempty(dat1)
        ind0=find(~ismember(dat0,time))';
        ind1=find(~ismember(dat1,time))';
        ind2=[ind0,ind1];
        ind=unique(ind2);
        dat0(ind)=[];
        dat1(ind)=[];

        dat.t0s{i,1}(ind)=[];
        dat.t1s{i,1}(ind)=[];
        %get the index for t0s and t1s based on et
        if isempty(dat0)||isempty(dat1)
            dat.t0_idx{i,1}=[];
            dat.t1_idx{i,1}=[];
        else
            idx0=find(ismember(time,dat0));
            idx1=find(ismember(time,dat1));
            dat.t0_idx{i,1}=idx0;
            dat.t1_idx{i,1}=idx1;
        end
    else
        dat.t0_idx{i,1}=[];
        dat.t1_idx{i,1}=[];
    end
    clear idx* dat1 dat0
end
%% 5) calculate the orientation for each datapoints using the spine data


for i=1:length(dat.xspine)
    x1=dat.xspine{i,1}(:,1);
    x2=dat.xspine{i,1}(:,6);
    y1=dat.yspine{i,1}(:,1);
    y2=dat.yspine{i,1}(:,6);

    vector(:,1)=odor_pos(1)-x2;
    vector(:,2)=odor_pos(2)-y2;

    deg{i,1}=cal_ori_deg(x1,x2,y1,y2,vector);
    clear x1 x2 y1 y2 vector
end

dat.orientation=deg;
clear deg;
%% 6) combine and delete some short turning events
for i=1:length(dat.t0_idx)
    if length(dat.t0_idx{i,1})>1&&length(dat.t1_idx{i,1})>1

        turn_interval=dat.t0_idx{i,1}(2:end)-dat.t1_idx{i,1}(1:end-1);
        %if the difference between two turning events is less than 5
        %sections, we combine them
        ind=find(turn_interval<5);
        dat.t0_idx{i,1}(ind+1)=[];
        dat.t1_idx{i,1}(ind)=[];
        dat.t0s{i,1}(ind+1)=[];
        dat.t1s{i,1}(ind)=[];

        turn_diff=dat.t1_idx{i,1}-dat.t0_idx{i,1};
        %after combination if the turning event takes less than 4
        %timepoints to turn, delete it
        idx=find(turn_diff<5);
        dat.t1_idx{i,1}(idx)=[];
        dat.t0_idx{i,1}(idx)=[];
        dat.t1s{i,1}(idx)=[];
        dat.t0s{i,1}(idx)=[];
    end
end

%% 7) get the index for running event and also delete turning envet happens
%right at the beginning or end of the trajectory
for i=1:length(dat.t0_idx)
    if isempty(dat.t0_idx{i,1})||isempty(dat.t1_idx{i,1})
        dat.r0_idx{i,1}=1;dat.r0s{i,1}=dat.et{i,1}(1);
        dat.r1_idx{i,1}=length(dat.et{i,1});dat.r1s{i,1}=dat.et{i,1}(dat.r1_idx{i,1});
        continue
    end
    idx0=dat.t0_idx{i,1};
    idx1=dat.t1_idx{i,1};
    len=length(dat.et{i,1})-2;
    dat.r0_idx{i,1}=idx1+1;
    dat.r1_idx{i,1}=idx0-1;

    if idx0(1)<3 && idx1(end)>len
        %if a turning event occur at the beginning or the end of a
        %trajectory, delete it
        dat.r1_idx{i,1}(1)=[];dat.r0_idx{i,1}(end)=[];
        dat.t0_idx{i,1}(1)=[];dat.t1_idx{i,1}(1)=[];dat.t0s{i,1}(1)=[];dat.t1s{i,1}(1)=[];
        if ~isempty(dat.t0_idx{i,1})
            dat.t0_idx{i,1}(end)=[];dat.t1_idx{i,1}(end)=[];dat.t0s{i,1}(end)=[];dat.t1s{i,1}(end)=[];
        end
    elseif idx0(1)<3 && idx1(end)<=len
        %if there is a turning event in the beginning but not in the end of
        %the trajecoty
        dat.r1_idx{i,1}(1)=[];dat.r1_idx{i,1}=[dat.r1_idx{i,1};len+2];
        dat.t0_idx{i,1}(1)=[];dat.t1_idx{i,1}(1)=[];dat.t0s{i,1}(1)=[];dat.t1s{i,1}(1)=[];
    elseif idx0(1)>=3 && idx1(end)>len
        %if there is a turning event in not the beginning but in the end of
        %the trajecoty
        dat.r0_idx{i,1}(end)=[];dat.r0_idx{i,1}=[1;dat.r0_idx{i,1}];
        dat.t0_idx{i,1}(end)=[];dat.t1_idx{i,1}(end)=[];dat.t0s{i,1}(end)=[];dat.t1s{i,1}(end)=[];
    else
        dat.r1_idx{i,1}=[dat.r1_idx{i,1};len+2];dat.r0_idx{i,1}=[1;dat.r0_idx{i,1}];
    end
    dat.r0s{i,1}=dat.et{i,1}(dat.r0_idx{i,1});
    dat.r1s{i,1}=dat.et{i,1}(dat.r1_idx{i,1});
end
%% 8) get reorientation angle, x position where turning occurs, pre/post-turning angles
dis=5;
%dis is used for pre/post-turning angles (we use the point 'dis' before the turning event to calculate the vector)
for i=1:length(dat.t0_idx)
    if isempty(dat.t0_idx{i,1})||isempty(dat.t1_idx{i,1})
        turn_x{i,1}=[];
        turn_y{i,1}=[];
        pre{i,1}=[];post{i,1}=[];
        continue;
    end
    ind=find(dat.t0_idx{i,1}==1|dat.t1_idx{i,1}==length(dat.x{i,1}));%if a turning event happens in the first and the last point of the trajectory, delete it
    dat.t0_idx{i,1}(ind)=[];dat.t1_idx{i,1}(ind)=[];
    dat.t0s{i,1}(ind)=[];dat.t1s{i,1}(ind)=[];
    clear ind idx0 idx1

    if isempty(dat.t0_idx{i,1})||isempty(dat.t1_idx{i,1})
        turn_x{i,1}=[];
        turn_y{i,1}=[];
        pre{i,1}=[];post{i,1}=[];
        continue;
    end

    idx0=dat.t0_idx{i,1};
    idx1=dat.t1_idx{i,1};
    %get the x position where turning occurs
    x0=dat.x{i,1}(idx0);
    x1=dat.x{i,1}(idx1);
    turn_x{i,1}=x0+(x1-x0)./2;
    y0=dat.y{i,1}(idx0);
    y1=dat.y{i,1}(idx1);
    turn_y{i,1}=y0+(y1-y0)./2;

    clear x0 x1 y0 y1
    %% calculate the pre_deg
    idx0_pre=idx0-dis;
    neg_val=find(idx0_pre<1); %if there is a negative value of index, change it to 1, so that we can calculate the pre-deg
    idx0_pre(neg_val)=1;

    x1=dat.x{i,1}(idx0);% then use the x and y value to get a vector
    y1=dat.y{i,1}(idx0);

    x2=dat.x{i,1}(idx0_pre);
    y2=dat.y{i,1}(idx0_pre);

    vector(:,1)=odor_pos(1)-x2;
    vector(:,2)=odor_pos(2)-y2;

    deg=cal_ori_deg(x1,x2,y1,y2,vector);
    pre{i,1}=deg;
    clear x1 y1 x2 y2 neg_val deg vector
    %% calculate the post_deg
    idx1_post=idx1+dis;
    pos_val=find(idx1_post>length(dat.x{i,1})); %if the index is longer than the total length, change it
    idx1_post(pos_val)=length(dat.x{i,1});

    x2=dat.x{i,1}(idx1);% then use the x and y value to get a vector
    y2=dat.y{i,1}(idx1);

    x1=dat.x{i,1}(idx1_post);
    y1=dat.y{i,1}(idx1_post);

    vector(:,1)=odor_pos(1)-x2;
    vector(:,2)=odor_pos(2)-y2;

    deg=cal_ori_deg(x1,x2,y1,y2,vector);
    post{i,1}=deg;
    clear x1 y1 x2 y2 deg pos_val idx* vector
end
dat.pre_deg=pre; dat.post_deg=post;dat.turn_x=turn_x;dat.turn_y=turn_y;
clear pre post turn_x turn_y
%% 9) get the change in reorientation for each turning event

for i=1:length(dat.pre_deg)
    if isempty(dat.pre_deg{i,1})||isempty(dat.post_deg{i,1})
        reori_deg{i,1}=[];
        continue
    end
    %if the pre/post angles are in second or third quadrant ([90 180] and [-180 90])
    %we the reorientation angle will be 360-abs(pre)-abs(post)
    %the rest will just be post-pre
    reori_deg{i,1}=abs(dat.post_deg{i,1}-dat.pre_deg{i,1});
    idx=find(reori_deg{i,1}>180);
    reori_deg{i,1}(idx)=360-reori_deg{i,1}(idx);

    clear idx
end
dat.reorient_deg_abs=reori_deg;
clear reori_deg

%% 10) Calculate properties for running event (length, position, and degree)

for i=1:length(dat.x)
    %calculate the length of run
    run_et{i,1}=dat.r1s{i,1}-dat.r0s{i,1};
    %get the where run starts
    r0x{i,1}=dat.x{i,1}(dat.r0_idx{i,1});

    x2=dat.x{i,1}(dat.r0_idx{i,1});
    x1=dat.x{i,1}(dat.r1_idx{i,1});

    y2=dat.y{i,1}(dat.r0_idx{i,1});
    y1=dat.y{i,1}(dat.r1_idx{i,1});

    vector(:,1)=odor_pos(1)-x2;
    vector(:,2)=odor_pos(2)-y2;

    deg=cal_ori_deg(x1,x2,y1,y2,vector);
    run_deg{i,1}=deg;
    clear deg x1 x2 x3 y1 y2 y3 vector

end
dat.run_deg=run_deg; dat.run_et=run_et; dat.r0x=r0x;
clear run* r0x

%% 11) eliminate data that has 99% of datapoints inside 3l* midline
for i=1:length(dat.x)
    x=dat.x{i,1}(1,1);
    y=dat.y{i,1}(1,1);
    L=dat.mid(i,1);
    l=length(dat.x{i,1});
    dis=sqrt((dat.x{i,1}-x).^2+(dat.y{i,1}-y).^2);
    ind=find(dis>=3*L);
    p_out_L(i,1)=length(ind)/l;
    clear x y ind dis L l
end
idx=find(p_out_L<0.01);

dat(idx,:)=[];
clear idx* p_out_L

%% save the data file for each condition

filename=fullfile(outdir,'data.mat');
if ~isfolder(outdir)
    mkdir(outdir)
end

if isfile(filename)
    delete(filename);
end
save(filename,'dat','-v7.3');
disp('data has been sucessfully processed')

clear
end