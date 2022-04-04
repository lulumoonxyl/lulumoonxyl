ii=1;%ii indicates which condition you want to load
%
outdir="C:\Users\feihu\OneDrive - McGill University\matlab\Analysis based on larvae\Fall_2021_Replication\All\example";
%% Load data from JB and JAABA
choredir_JB='C:\Users\feihu\OneDrive\Desktop\JB_data';
[dat_JB,uname_JB]=get_JB_data(choredir_JB,ii);
%get x and y
dat_JB=get_xy_center(dat_JB);

%load JAABA data
choredir_JAABA='C:\Users\feihu\Downloads\JAABA_data\jaaba_result';
% field t0 is when turning starts in second, t1 is when turning events end
% in second
[dat_JAABA,uname_JAABA]=get_JAABA_data(choredir_JAABA,ii);
%% rearrange the data array for the jaaba
dat_JAABA=arrange_JAABA_data(dat_JAABA);%combine the cells in each fields
dat_JAABA=align_JB_JAABA_data(dat_JB,dat_JAABA);%delete the data that is in JAABA but not in JB
%find the index for t0s and t0s
[dat_JAABA.t0_idx,dat_JAABA.t1_idx,dat_JAABA.t0s,dat_JAABA.t1s]=find_ind_t0s_t1s(dat_JB.et,dat_JAABA.t0s,dat_JAABA.t1s);
%% calculate the heading direction at each timepoint using xspine and
%%yspine
vector=[-1 0];
dat_JB.heading_direction=calculate_heading_direction(dat_JB.xspine,dat_JB.yspine,vector);
%% calculate the pre/post turning degree using x and y
dat_JAABA=cal_deg_pre_post_turn(vector,dat_JB,dat_JAABA);
%get the x, y position for turning event
[dat_JAABA.turn_x,dat_JAABA.turn_y]=get_tunring_xy_pos(dat_JB.x,dat_JB.y,dat_JAABA.t0_idx,dat_JAABA.t1_idx);
%% calculate the absolute value of reorientation angle
dat_JAABA.reorient_deg=cal_reorient_deg(dat_JAABA.pre_deg,dat_JAABA.post_deg);
%% calculate the turning frequency(tf_mean), 
% turning frequency larger than deg(60,tfl_mean), abs reorientation angle (deg_mean) based on x position
deg=60; 
bins=20;
[tf_mean,tfl_mean,deg_mean]=cal_reorient_f(dat_JAABA.turn_x, ...
    dat_JB.x,dat_JB.et,dat_JAABA.reorient_deg,bins,deg,'x'); %keep in mind that odor is at x=0

%% same as above but based on heading direction 
bins=10;
[tf_mean_hd,tfl_mean_hd,deg_mean_hd]=cal_reorient_f(dat_JAABA.pre_deg, ...
    dat_JB.heading_direction,dat_JB.et,dat_JAABA.reorient_deg,bins,deg,'t');
%% plot trajectory
list=[1 2];% just which animal you want to plot 
plot_specific_trajectory(list,dat_JB.x,dat_JB.y,dat_JAABA.turn_x,dat_JAABA.turn_y,dat_JAABA.t0_idx,dat_JAABA.t1_idx);
%%  save figures
save_all_figures(outdir);
%% save data
file_dir=fullfile(outdir,'Data');

if ~isfolder(file_dir)
    mkdir(file_dir);
end
filename=fullfile(file_dir,'dat.mat');
if isfile(filename)
delete(filename);
end
save( filename, 'dat*')
